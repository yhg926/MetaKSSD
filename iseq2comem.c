#include "iseq2comem.h"

#include "command_dist.h"
#include "global_basic.h"
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <zlib.h>
#include <math.h>
#include <string.h>
#include <err.h>
#include <errno.h>
#include <sys/mman.h>
#include <fcntl.h> // open function
#include <unistd.h> //close function



#define HIBITSET1 0x8000000000000000LLU
#define _64MASK 0xffffffffffffffffLLU
/* mv to global_basic.h 20220714
#define H1(K,HASH_SZ) ((K)%(HASH_SZ))
#define H2(K,HASH_SZ) ( 1 + (K) % ( (HASH_SZ) - 1 ) )
#define HASH(K,I,HASH_SZ) ( ( H1(K,HASH_SZ) + I * H2(K,HASH_SZ) ) % HASH_SZ )
*/
static int rand_id;
static int half_ctx_len;
static int half_subctx_len;
static int half_outctx_len;
static int drlevel;

static int comp_bittl;  //64-BITTL ;
static int crvsaddmove ; //bits the new base should right move before add to crvstuple
static llong tupmask; //TUPMASK;
static int TL ;  //TUPLEN ;
//MASK for Dimension Reduction
static  llong domask ;
static  llong undomask ;
//dim shuffled array arg
dim_shuffle_t* dim_shuffle; //dim stat + dim arr
static int *dim_shuf_arr ;
static int dim_shuf_arr_len;
static int dim_start; //reduction dim start
static int dim_end;  // reduction dim end

unsigned int hashsize;
unsigned int hashlimit; //limit context count in hashtable
int component_num;

void seq2co_global_var_initial(void)
{
	rand_id = dim_shuffle->dim_shuffle_stat.id ;  
	half_ctx_len = dim_shuffle->dim_shuffle_stat.k ;
	half_subctx_len = dim_shuffle->dim_shuffle_stat.subk ;
	half_outctx_len = half_ctx_len - half_subctx_len;
  drlevel = dim_shuffle->dim_shuffle_stat.drlevel;
	hashlimit = hashsize * LD_FCTR ;
	printf("rand_id=%d\thalf_ctx_len=%d\thashsize=%d\thashlimit=%d\n",rand_id,half_ctx_len,hashsize,hashlimit);
//	printf("%d\t%d\t%d\n",COMPONENT_SZ,half_ctx_len,drlevel );
	component_num = half_ctx_len - drlevel > COMPONENT_SZ ? 
								 1LU << 4*(half_ctx_len - drlevel - COMPONENT_SZ ) : 1  ; 
//	printf("compnum=%d\t compsize=%d\n",component_num,1<<4*COMPONENT_SZ) ;	
	comp_bittl = 64-4*half_ctx_len;  //64-BITTL ;
	crvsaddmove = 4*half_ctx_len-2; //4*half_ctx_len, bits the new base should right move before add to crvstuple
	tupmask = _64MASK >> comp_bittl; //TUPMASK;
	TL = 2*half_ctx_len ;  //TUPLEN ;

//MASK for Dimension Reduction
	//( (1LLU << ( half_subctx_len*2 ) ) - 1 ) << (2*half_ctx_len + 2);
	domask = ( (1LLU << ( half_subctx_len*4 ) ) - 1 ) << (2*half_outctx_len); 
	undomask =  ( (1LLU << (half_outctx_len*2)) - 1 ) 
									<< (2*(half_ctx_len + half_subctx_len)); //<< (2*(half_ctx_len + half_subctx_len + 1));	

	dim_shuf_arr = dim_shuffle->shuffled_dim;
	dim_shuf_arr_len = 1LLU << (4*half_subctx_len) ;
	dim_start = 0;
	dim_end = MIN_SUBCTX_DIM_SMP_SZ ;
};

// reads2mco, per read per co, combined into mco, repeat k-mer allowed 
int reads2mco(char* seqfname,const char *co_dir, char * pipecmd){
#define unit_incrs 1000000 //realloc increase unit 
int comp_code_bits = half_ctx_len - drlevel > COMPONENT_SZ ? 4*(half_ctx_len - drlevel - COMPONENT_SZ ) : 0  ;

	size_t **cof_count = malloc( component_num * sizeof (size_t*) ); 
	FILE **outindf = malloc(component_num * sizeof(FILE *));
	FILE **outf  = malloc(component_num * sizeof(FILE *));
	char indexfname[PATHLEN]; char combined_cof[PATHLEN];
	
	size_t cof_count_sz = unit_incrs;
	for(int i=0;i<component_num ;i++)
  {
		cof_count[i] = (size_t *)calloc( cof_count_sz , sizeof(size_t) );
		sprintf(combined_cof,"%s/combco.%d",co_dir,i);
		sprintf(indexfname,"%s/combco.index.%d",co_dir,i);
  	if( (outf[i] = fopen(combined_cof,"wb")) == NULL) err(errno,"%s",combined_cof);
    if( (outindf[i] = fopen(indexfname,"wb")) == NULL) err(errno,"%s",indexfname);				
  };	


	FILE *infp;
	char fas_fname[PATHLEN];
	if(pipecmd[0] != '\0'){
		sprintf(fas_fname,"%s %s",pipecmd,seqfname);
		if( (infp=popen(fas_fname,"r")) == NULL ) err(errno,"reads2mco():%s",fas_fname);
	}
	else
		if( (infp=fopen(seqfname,"r")) == NULL ) err(errno,"reads2mco():%s",fas_fname);;
	
	
	char seqin_buff[ READSEQ_BUFFSZ + 1 ];
	int newLen =  fread(seqin_buff, sizeof(char),READSEQ_BUFFSZ,infp);	
	if(! (newLen >0) )  err(errno,"reads2mco():eof or fread error file=%s",seqfname);

	llong base = 1; char ch; int basenum;
	llong tuple = 0LLU, crvstuple = 0LLU,
  unituple, drtuple, pfilter;
	llong readn = 0; //initial readn;
	for(int pos = 0; pos <= newLen; pos++)
  {
    if(pos == newLen){
        newLen =  fread(seqin_buff, sizeof(char),READSEQ_BUFFSZ,infp);
      if(newLen > 0)
        pos = 0;
      else break;
    };
    ch = seqin_buff[pos];
    basenum = Basemap[(int)ch];

    if(basenum != DEFAULT) //make sure basenum is not negative
    {
      tuple = ( ( tuple<< 2 ) | (llong)basenum ) & tupmask ;
      crvstuple = ( crvstuple >> 2 ) + (((llong)basenum^3LLU) << crvsaddmove);
      base++;
    }
    else if ( (ch == '\n') || (ch == '\r') ) { continue;}
    else if (isalpha(ch)){ base=1; continue; }
    else if ( ch == '>' )
    {
			if( readn > cof_count_sz ){ //realloc
				for(int i=0;i<component_num ;i++){
					size_t *newtmp = (size_t *)realloc(cof_count[i], (cof_count_sz + unit_incrs) * sizeof(size_t));
					if (newtmp != NULL) {
 						cof_count[i] = newtmp;
						memset(cof_count[i] + cof_count_sz  ,0, unit_incrs * sizeof(size_t) );	
					}
					else err(errno,"cof_count[%d] realloc failed",i);
				}	
				cof_count_sz += unit_incrs;
			}	
			readn++;
 
     while( (pos < newLen ) && ( seqin_buff[pos] != '\n' ) )
      {
        if (pos < newLen - 1)
          pos++;
        else
        {
          newLen = fread(seqin_buff, sizeof(char),READSEQ_BUFFSZ,infp);
          if(newLen > 0) pos = -1; // pos = 0 is bug ? should be pos = -1 , since continue and pos++
          else err(errno,"fasta2co(): can not find seqences head start from '>' %d",newLen);
        };
      };
      base = 1;
      continue;
    }
    else {
    //  warnx("ignorn illegal charactor%c \n",ch);
      base=1;
      continue;
    };

    if( base > TL ) // if(!(base < TL))
    {
      //make sure unituple == min(tuple,crvstuple);
      unituple = tuple < crvstuple ? tuple:crvstuple;
      //only for 64bits Kmer storage ,
      //important !!!! make sure is right
      int dim_tup = ((unituple & domask) >> ( (half_outctx_len)*2 ) ) ;
      pfilter = dim_shuf_arr[dim_tup];
      if( ( pfilter >= dim_end) || (pfilter < dim_start ) ) continue;
      pfilter = pfilter - dim_start; //add 190307: prevent pfilter > MIN_SUBCTX_DIM_SMP_SZ when dim_start >0
      drtuple = ( ( (unituple & undomask) //left half outer ctx
              + ( ( unituple & ( ( 1LLU<< ( half_outctx_len*2) ) - 1)) << (TL*2 - half_outctx_len*4) ) )
              >> ( drlevel*4 ) ) // subctx dim reduced
              +  pfilter ; //  subctx dim after reduction
	
			cof_count[drtuple % component_num][readn]++;	//ctx_count for each components each read				
			unsigned int newid = (unsigned int)( drtuple >> comp_code_bits ) ;
			fwrite( &newid, sizeof(unsigned int),1,outf[(int)( drtuple % component_num )] );				 			
   };
  };//end file hashing
	
	for(int i=0;i<component_num ;i++){ //write index
		llong cumval = 0 ; // cumulated value 
		for(llong n=0; n<=readn;n++ ){
			cumval += cof_count[i][n]; 
			fwrite(&cumval, sizeof(llong),1,outindf[i]);
		}

		fclose(outf[i]);
		fclose(outindf[i]);
	}
	printf("decomposing %s by reads is complete!\n",seqfname);  	
	return 1;
}

const char gzpipe_cmd[]= "zcat -fc"; //const char gzpipe_cmd[]= "unpigz -fc";

llong * fasta2co(char* seqfname, llong *co, char * pipecmd) //20190910, enhancement: pipecmd 
{
	llong tuple = 0LLU, crvstuple = 0LLU,
    unituple, drtuple, pfilter;

	memset(co,0LLU,hashsize*sizeof(llong));		
	char seqin_buff[ READSEQ_BUFFSZ + 1 ]; // seq infile buffer
	FILE *infp;
	char fas_fname[PATHLEN];
	if(pipecmd[0] != '\0')
		sprintf(fas_fname,"%s %s",pipecmd,seqfname);//other pipecmd except decompress cmd
	else
		sprintf(fas_fname,"%s %s",gzpipe_cmd,seqfname);

	if( (infp=popen(fas_fname,"r")) == NULL ) err(errno,"fasta2co():%s",fas_fname);

	int newLen =  fread(seqin_buff, sizeof(char),READSEQ_BUFFSZ,infp);
	if(! (newLen >0) ) 	err(errno,"fastco():eof or fread error file=%s",seqfname);

	llong base = 1; char ch; int basenum;
	unsigned int keycount = 0;
	// core function begin
	for(int pos = 0; pos <= newLen; pos++) 
	{	
		if(pos == newLen){
				newLen =  fread(seqin_buff, sizeof(char),READSEQ_BUFFSZ,infp);	
			if(newLen > 0)
				pos = 0;
			else break;
		};
		ch = seqin_buff[pos];
		basenum = Basemap[(int)ch];

		if(basenum != DEFAULT) //make sure basenum is not negative
		{
			tuple = ( ( tuple<< 2 ) | (llong)basenum ) & tupmask ;
			crvstuple = ( crvstuple >> 2 ) + (((llong)basenum^3LLU) << crvsaddmove);
			base++;
		}
		else if ( (ch == '\n') || (ch == '\r') ) { continue;}
		else if (isalpha(ch)){ base=1; continue; }
		else if ( ch == '>' )
		{
			while( (pos < newLen ) && ( seqin_buff[pos] != '\n' ) )
			{			
				if (pos < newLen - 1)
					pos++;
				else
				{
					newLen = fread(seqin_buff, sizeof(char),READSEQ_BUFFSZ,infp);
					if(newLen > 0) pos = -1; // pos = 0 is bug ? should be pos = -1 , since continue and pos++
					else err(errno,"fasta2co(): can not find seqences head start from '>' %d",newLen);
				};
			};
      base = 1;
      continue;
		}
		else {
    //  warnx("ignorn illegal charactor%c \n",ch);
      base=1;
      continue;
    };

		if( base > TL ) // if(!(base < TL))
		{
			//make sure unituple == min(tuple,crvstuple);
			unituple = tuple < crvstuple ? tuple:crvstuple;
			//only for 64bits Kmer storage ,
    	//important !!!! make sure is right
			int dim_tup = ((unituple & domask) >> ( (half_outctx_len)*2 ) ) ;
			pfilter = dim_shuf_arr[dim_tup];
			if( ( pfilter >= dim_end) || (pfilter < dim_start ) ) continue;
			pfilter = pfilter - dim_start; //add 190307: prevent pfilter > MIN_SUBCTX_DIM_SMP_SZ when dim_start >0	
			drtuple = ( ( (unituple & undomask) //left half outer ctx
							+ ( ( unituple & ( ( 1LLU<< ( half_outctx_len*2) ) - 1)) << (TL*2 - half_outctx_len*4) ) )										
            	>> ( drlevel*4 ) ) // subctx dim reduced
            	+  pfilter ; //  subctx dim after reduction

			unsigned int i,n ;
    	for(i=0;i<hashsize;i++)
    	{
      	n = HASH(drtuple,i,hashsize);
      	if (co[n] == 0)
      	{
        	co[n] = drtuple;
        	keycount++;
        	if( keycount > hashlimit)
          	err(errno,"the context space is too crowd, try rerun the program using -k%d", half_ctx_len + 1);
        	break;
      	}
      	else if ( co[n] == drtuple )
          break;
			};//end kmer hashing
		};
	};//end file hashing
		pclose(infp);
  return co;
};//end func
// macro for fastq2co
// Q:quality score, M: required least occurence of kmer
// make sure LEN is enough for long read
#define LEN 20000 //4096, 20220720 update for pacbio reads
#define CT_BIT 4 //bits for Kmer count
#define CT_MAX 0xfLLU //make sure smaller than 1LLU<<CT_BTI

llong * fastq2co(char* seqfname, llong *co, char *pipecmd, int Q, int M ) //20190910 enhanced pipecmd
{
	if(M >= CT_MAX) err(errno,"fastq2co(): Occurence num should smaller than %d", (int)CT_MAX);

	llong tuple = 0LLU, crvstuple = 0LLU, unituple, drtuple, pfilter;	
	memset(co,0LLU,hashsize*sizeof(llong));

	FILE *infp;
	char fq_fname[PATHLEN];
	if(pipecmd[0] != '\0')
    sprintf(fq_fname,"%s %s",pipecmd,seqfname);//other pipecmd except decompress cmd
  else
    sprintf(fq_fname,"%s %s",gzpipe_cmd,seqfname);

	if( (infp=popen(fq_fname,"r")) == NULL ) err(errno,"fastq2co():%s",fq_fname);
	//char seq[LEN];
	//char qual[LEN];
	char *seq = malloc(LEN+10);
	char *qual = malloc(LEN+10);
	fgets(seq,LEN,infp); fgets(seq,LEN,infp); 
	fgets(qual,LEN,infp);  fgets(qual,LEN,infp);
//	if(Q != -1) { 		};
	llong base = 1; char ch ; int basenum,line_num = 0 ;
	unsigned int keycount =0 ;
	int sl = strlen(seq); 
	for(int pos = 0; pos < sl; pos++){
		if(seq[pos] == '\n' ){
			fgets(seq,LEN,infp); fgets(seq,LEN,infp);
			fgets(qual,LEN,infp); fgets(qual,LEN,infp);
			sl = strlen(seq);
//			printf("seq=%s\nqual=%s\n==\n",seq,qual); 
			line_num+=4;
			if( !feof(infp) ) {
				base = 1; 
				pos = -1;
				continue ;									
			}
			else break;	
		}		
		else{
			ch = seq[pos];
			basenum = Basemap[(int)ch];
			if( (basenum != DEFAULT) && ( qual[pos] >= Q ) ){
				tuple = ( ( tuple<< 2 ) | (llong)basenum ) & tupmask ;
				crvstuple = ( crvstuple >> 2 ) + (((llong)basenum^3LLU) << crvsaddmove);
				base++;
			} 

			else {
			//	if (qual[pos] < Q)
			//		warnx("in file %s line %d skip low quality charactor '%c[%c]'\n",seqfname, line_num+2 ,ch,qual[pos]);
		 //	else warnx("in file %s line %d ignorn illegal charactor '%c'\n",seqfname,line_num+2,ch);
				base = 1;
				continue;
			};
		};
	
		if( base > TL ){ // serious bug ?!!!:  if( base >= TL ) 

			unituple = tuple < crvstuple ? tuple:crvstuple;
			int dim_tup = ((unituple & domask) >> ( (half_outctx_len)*2 ) ) ;
      pfilter = dim_shuf_arr[dim_tup];
      if( ( pfilter >= dim_end) || (pfilter < dim_start ) ) continue;
			pfilter = pfilter - dim_start;
      drtuple = ( ( (unituple & undomask) //left half outer ctx
              + ( ( unituple & ( ( 1LLU<< ( half_outctx_len*2) ) - 1)) << (TL*2 - half_outctx_len*4) ) )
              >> ( drlevel*4 ) ) // subctx dim reduced
              +  pfilter ;
			
			unsigned int i,n ;
			for(i=0;i<hashsize;i++)	{
				n = HASH(drtuple, i, hashsize);
				if (co[n] == 0LLU){
					//20190910:bug fixed, otherwise occrence==1 k-mer would be skipped
					if( M == 1) co[n] = (drtuple << CT_BIT) | CT_MAX; //--
					else co[n] = (drtuple << CT_BIT) + 1LLU;					
          if( keycount > hashlimit)
            err(errno,"the context space is too crowd, try rerun the program using -k%d", half_ctx_len + 1);
					break;
				} 
				else if ( ( co[n] >> CT_BIT ) == drtuple ) {
					//20190910:bug fixed, test if already to CT_MAX
	          if( (co[n] & CT_MAX) == CT_MAX ) break;
						co[n] += 1LLU;
						if( !((co[n] & CT_MAX) <  M) ) co[n]|= CT_MAX ;
						break ;					
        };
			}; //end kmer hashing
		};
	}// end file hashing 
	printf("%d reads detected\n",line_num);

	pclose(infp);
	free(seq);
  free(qual);
	return co;
}; // end fastq2co()

// fastq2co with occurence
// higher 48 bits for Kmer, lower 16 bits for count
#define OCCRC_BIT 16 // make sure it is 64 bits machine
#define OCCRC_MAX 0xffffLLU //must be 1<<OCCRC_BIT - 1
//20190910 enhanced pipecmd
llong * fastq2koc (char* seqfname, llong *co, char *pipecmd, int Q) 
{

  llong tuple = 0LLU, crvstuple = 0LLU, unituple, drtuple, pfilter;
  memset(co,0LLU,hashsize*sizeof(llong));

  FILE *infp;
  char fq_fname[PATHLEN];
	if(pipecmd[0] != '\0')
    sprintf(fq_fname,"%s %s",pipecmd,seqfname);//other pipecmd except decompress cmd
  else
    sprintf(fq_fname,"%s %s",gzpipe_cmd,seqfname);

  if( (infp=popen(fq_fname,"r")) == NULL ) err(errno,"fastq2koc():%s",fq_fname);

  //char seq[LEN];
  //char qual[LEN];
	char *seq = malloc(LEN+10);
  char *qual = malloc(LEN+10);

  fgets(seq,LEN,infp); fgets(seq,LEN,infp);
  fgets(qual,LEN,infp);  fgets(qual,LEN,infp);

  llong base = 1; char ch ; int basenum,line_num = 0 ;
  unsigned int keycount =0 ;
	int sl = strlen(seq); //20220720: old code put strlen(seq) in loop largely reduce the speed 
  for(int pos = 0; pos < sl; pos++){
    if(seq[pos] == '\n' ){
      fgets(seq,LEN,infp); fgets(seq,LEN,infp);
      fgets(qual,LEN,infp); fgets(qual,LEN,infp);
			sl = strlen(seq);
      line_num+=4;
      if( !feof(infp) ) {
        base = 1;
        pos = -1;
        continue ;
      }
      else break;
    }
    else{
      ch = seq[pos];
      basenum = Basemap[(int)ch];
      if( (basenum != DEFAULT) && ( qual[pos] >= Q ) ){
        tuple = ( ( tuple<< 2 ) | (llong)basenum ) & tupmask ;
        crvstuple = ( crvstuple >> 2 ) + (((llong)basenum^3LLU) << crvsaddmove);
        base++;
      }

      else {
        base = 1;
        continue;
      };
    };

    if( base > TL ){ 
      unituple = tuple < crvstuple ? tuple:crvstuple;
      unsigned int dim_tup = ((unituple & domask) >> ( (half_outctx_len)*2 ) ) ;
      pfilter = dim_shuf_arr[dim_tup];
      if( ( pfilter >= dim_end) || (pfilter < dim_start ) ) continue;
      pfilter = pfilter - dim_start;
      drtuple = ( ( (unituple & undomask) //left half outer ctx
              + ( ( unituple & ( ( 1LLU<< ( half_outctx_len*2) ) - 1)) << (TL*2 - half_outctx_len*4) ) )
              >> ( drlevel*4 ) ) // subctx dim reduced
              +  pfilter ;

      unsigned int i,n ;
      for(i=0;i<hashsize;i++) {
        n = HASH(drtuple, i, hashsize);
        if (co[n] == 0LLU){
          co[n] = (drtuple << OCCRC_BIT) + 1LLU ;
          keycount++;
          if( keycount > hashlimit )
            err(errno,"the context space is too crowd, try rerun the program using -k%d", half_ctx_len + 1);
          break;
        } else if ( ( co[n] >> OCCRC_BIT ) == drtuple ) {

            if( (co[n] & OCCRC_MAX) < OCCRC_MAX )
              co[n]+=1LLU;
            //else co[n]|= OCCRC_MAX ;
          break ;
        };
      }; //end kmer hashing
    };
  }// end file hashing
  pclose(infp);
	free(seq);
	free(qual);
  return co;
}; // end fastq2koc();

unsigned int write_fqkoc2files(char* cofilename, llong *co)
{
  int comp_code_bits = half_ctx_len - drlevel > COMPONENT_SZ ? 4*(half_ctx_len - drlevel - COMPONENT_SZ ) : 0  ;
  FILE **outf,**abdf; //sketch and occurence files
  outf = malloc(component_num * sizeof(FILE *));
	abdf = malloc(component_num * sizeof(FILE *));
  char cofilename_with_component[PATHLEN];

  for(int i=0;i<component_num ;i++)
  {
    sprintf(cofilename_with_component,"%s.%d",cofilename,i);
    if ( (outf[i] = fopen(cofilename_with_component,"wb") ) == NULL )
      err(errno,"write_fqkoc2files()") ;
	
		sprintf(cofilename_with_component,"%s.%d.a",cofilename,i);
		if ( (abdf[i] = fopen(cofilename_with_component,"wb") ) == NULL )
      err(errno,"write_fqkoc2files()") ;
  };

  unsigned int count, wr = 0, newid,compi; //compi: component index
	unsigned short abdc; //abundance 
		

  for(count=0;count < hashsize; count++)
  {
    if( co[count] > 0 ) {
			
			compi = (co[count] >> OCCRC_BIT ) % component_num;

      newid = (unsigned int)(co[count] >> (comp_code_bits + OCCRC_BIT));
      fwrite( &newid, sizeof(newid),1,outf[compi] );
			
			abdc = co[count] & OCCRC_MAX;
			fwrite( &abdc, sizeof(abdc),1,abdf[compi] );
		
      wr++;
    }
  }

  for(int i=0;i<component_num ;i++){
    fclose(outf[i]);
		fclose(abdf[i]);
	}
  free(outf);
	free(abdf);
  return wr;
};

llong write_fqkoc2file(char* cofilename, llong *co)
{
  int comp_code_bits = half_ctx_len - drlevel > COMPONENT_SZ ? 4*(half_ctx_len - drlevel - COMPONENT_SZ ) : 0  ;
  FILE **outf;
  outf = malloc(component_num * sizeof(FILE *));
  char cofilename_with_component[PATHLEN];
  for(int i=0;i<component_num ;i++)
  {
    sprintf(cofilename_with_component,"%s.%d",cofilename,i);
    if ( (outf[i] = fopen(cofilename_with_component,"wb") ) == NULL )
      err(errno,"write_fqkoc2file()") ;
  };

  unsigned int count, wr = 0;
	llong newid;
  for(count=0;count < hashsize; count++)
  {
    if( co[count] > 0 ) {
			newid = ((co[count] >> (comp_code_bits + OCCRC_BIT)) << OCCRC_BIT ) | (co[count] & OCCRC_MAX) ;	
      fwrite( &newid, sizeof(llong),1,outf[ (co[count] >> OCCRC_BIT ) % component_num ] );
      wr++;
    }
  }

  for(int i=0;i<component_num ;i++)
    fclose(outf[i]);

	free(outf);
  return wr;
};

// need special write fastq 2 co 
llong write_fqco2file(char* cofilename, llong *co)
{
	int comp_code_bits = half_ctx_len - drlevel > COMPONENT_SZ ? 4*(half_ctx_len - drlevel - COMPONENT_SZ ) : 0  ;
	FILE **outf;
	outf = malloc(component_num * sizeof(FILE *));
	char cofilename_with_component[PATHLEN];
	for(int i=0;i<component_num ;i++)
  {
    sprintf(cofilename_with_component,"%s.%d",cofilename,i);
    if ( (outf[i] = fopen(cofilename_with_component,"wb") ) == NULL )
      err(errno,"write_fqco2file()") ;
  };
	unsigned int count, wr = 0, newid;
	for(count=0;count < hashsize; count++)
	{
		if( (co[count] & CT_MAX) == CT_MAX ) { 		
			newid =  (unsigned int)( co[count] >> (comp_code_bits + CT_BIT) ) ;					
			fwrite( &newid, sizeof(unsigned int),1,outf[(int)( (co[count] >> CT_BIT ) % component_num )] );
      wr++;				
		}	
	}
	for(int i=0;i<component_num ;i++)
  	fclose(outf[i]);
	free(outf);
	return wr;
};

//write co to compnent files
//component determeined by modula context,namely use inner subctx
llong wrt_co2cmpn_use_inn_subctx(char* cofilename, llong *co)
{	
	int comp_code_bits = half_ctx_len - drlevel > COMPONENT_SZ ? 4*(half_ctx_len - drlevel - COMPONENT_SZ ) : 0  ;  
	FILE **outf;
	outf = malloc(component_num * sizeof(FILE *));
	char cofilename_with_component[PATHLEN];
  for(int i=0;i<component_num ;i++)
  {
    sprintf(cofilename_with_component,"%s.%d",cofilename,i);
    if ( (outf[i] = fopen(cofilename_with_component,"wb") ) == NULL )
      err(errno,"wrt_co2cmpn_use_inn_subctx()") ;
  };
  unsigned int count, wr = 0, newid;
	for(count=0;count < hashsize; count++)
	{
		if( co[count] != 0 )
		{	
			newid = (unsigned int)( co[count] >> comp_code_bits ) ;
			fwrite( &newid, sizeof(unsigned int),1,outf[(int)( co[count] % component_num )] );
			wr++;
		}
	}
	for(int i=0;i<component_num ;i++)
   fclose(outf[i]);

	free(outf);
  return wr;
};


