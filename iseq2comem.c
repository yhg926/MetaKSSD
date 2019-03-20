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
#define H1(K,HASH_SZ) ((K)%(HASH_SZ))
#define H2(K,HASH_SZ) ( 1 + (K) % ( (HASH_SZ) - 1 ) )
#define HASH(K,I,HASH_SZ) ( ( H1(K,HASH_SZ) + I * H2(K,HASH_SZ) ) % HASH_SZ )

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


const char gzpipe_cmd[]= "zcat -fc"; //const char gzpipe_cmd[]= "unpigz -fc";

llong * fasta2co(char* seqfname, llong *co)
{
	llong tuple = 0LLU, crvstuple = 0LLU,
    unituple, drtuple, pfilter;

	memset(co,0LLU,hashsize*sizeof(llong));		
	char seqin_buff[ READSEQ_BUFFSZ + 1 ]; // seq infile buffer
	FILE *infp;
	char fas_fname[PATHLEN];
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
#define LEN 4096 
#define CT_BIT 4 //bits for Kmer count
#define CT_MAX 0xfLLU //make sure smaller than 1LLU<<CT_BTI

llong * fastq2co(char* seqfname, llong *co, int Q, int M )
{
	if(M >= CT_MAX) err(errno,"fastq2co(): Occurence num should smaller than %d", (int)CT_MAX);

	llong tuple = 0LLU, crvstuple = 0LLU, unituple, drtuple, pfilter;	
	memset(co,0LLU,hashsize*sizeof(llong));

	FILE *infp;
	char fq_fname[PATHLEN];
	sprintf(fq_fname,"%s %s",gzpipe_cmd,seqfname);
	if( (infp=popen(fq_fname,"r")) == NULL ) err(errno,"fastq2co():%s",fq_fname);
	
	char seq[LEN];
	char qual[LEN];

	fgets(seq,LEN,infp); fgets(seq,LEN,infp); 
	fgets(qual,LEN,infp);  fgets(qual,LEN,infp);
//	if(Q != -1) { 		};
	llong base = 1; char ch ; int basenum,line_num = 0 ;
	unsigned int keycount =0 ;
	for(int pos = 0; pos < strlen(seq); pos++){
		if(seq[pos] == '\n' ){
			fgets(seq,LEN,infp); fgets(seq,LEN,infp);
			fgets(qual,LEN,infp); fgets(qual,LEN,infp); 
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
					co[n] = (drtuple << CT_BIT) + 1LLU ;					
					keycount++;
          if( keycount > hashlimit)
            err(errno,"the context space is too crowd, try rerun the program using -k%d", half_ctx_len + 1);
					break;
				} else if ( ( co[n] >> CT_BIT ) == drtuple ) {

						if( (co[n] & CT_MAX) < M )			
							co[n]+=1LLU; 	
						else 
							co[n]|= CT_MAX ; 
					break ;
        };
			}; //end kmer hashing
		};
	}// end file hashing 
	pclose(infp);
	return co;
}; // end fastq2co()
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
  return wr;
};


