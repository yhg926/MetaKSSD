#include "command_dist.h"
#include "command_shuffle.h"
#include "iseq2comem.h"
#include "co2mco.h"
#include "mytime.h"
#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <dirent.h>
#include <err.h>
#include <errno.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <sys/sysinfo.h>
#include <math.h>
#include <tgmath.h>
#include <time.h>
#include <unistd.h>
#include <zlib.h>
#include <stdbool.h>
#include <malloc.h>
#include <libgen.h>
#ifdef _OPENMP
	#include <omp.h>
#endif

mem_dispatch_t mem_dispatch = {0,0,0,0,0,0,0, OBJ_TRADITION , STD , MAKE_REF ,0};
mem_usage_stat_t mem_usage_stat = { 0, 0, 1e7 };
static llint real_time_mem;

const char co_dstat[] = "cofiles.stat"; //stat file name in co dir 
const char mco_dstat[] = "mcofiles.stat";
const char distoutdir[] = "dist";
const char logfpath[] = "dist_log.out";
FILE *logfp;

//distance compute type and core fun
typedef int mco_co_dist_t[BIN_SZ];
typedef unsigned int ctx_obj_ct_t; 
static inline unsigned int * mco_co_mmpdist_core(gidobj_t** unit_arrmco, char *co_fcode_in, unsigned int *ctx_obj_ct_in );
static inline void mco_co_dist_core(gidobj_t** unit_arrmco, char *co_fcode_in, int bin_sz,
            mco_co_dist_t shared_ctx_num_in);

//global distance threshold parameter, need initialize before use
static int neighbor_num ;
static double mutdistance_thre ;
//==============functions================//

/* dispatch functions  */
int dist_dispatch(dist_opt_val_t *opt_val) 
{
	logfp = fpathopen(opt_val->outdir, logfpath,"a") ;	
	// get available mem. from pre-set or system(if no pre-set)
	if( opt_val->mmry == 0 )
		opt_val->mmry = get_sys_mmry();
	llint base_mem = (llint)( (opt_val->mmry - get_sys_mmry()) * BBILLION ) - mem_usage_stat.others ;	
	real_time_mem = base_mem + (llint)(get_sys_mmry()*BBILLION) ; 
	fprintf(logfp,"availbe mem.=%lf\treal time mem.=%llu\n", get_sys_mmry(), real_time_mem);	
	
	int p_fit_mem = 1;	
	/***** invoke database search mode: *********/
	if( ( opt_val->refpath[0] != '\0' ) ) // refpath set
  {	
		const char *refco_dstat_fpath = test_get_fullpath(opt_val->refpath,co_dstat);	
		const char *refmco_dstat_fpath = test_get_fullpath(opt_val->refpath,mco_dstat);

 		//if not mco or co then prepare for stage I mode: convert .fas/.fq to .mco 
		if( (refco_dstat_fpath == NULL ) && ( refmco_dstat_fpath == NULL) ){

			infile_tab_t *ref_stat = dist_organize_refpath(opt_val);
			if(ref_stat->infile_num == 0)
				err(errno,"no valid input .fas/.fq file or absent %s | %s in %s \n", co_dstat, mco_dstat, opt_val->refpath ) ;	
			
			infile_fmt_count_t* ref_fmt_count;
      ref_fmt_count = infile_fmt_count(ref_stat);
			if( (ref_stat->infile_num == 0) || (ref_fmt_count->fasta + ref_fmt_count->fastq != ref_stat->infile_num ))
				err(errno,"not a valid input files: make sure input files are .fas/.fq format" 
					"or .co.num |.mco.num files with %s | %s in %s \n", co_dstat, mco_dstat, opt_val->refpath );

			const char * dist_rslt_dir = mk_dist_rslt_dir(opt_val->outdir,"kssd_dist_rslt");
      const char *dist_refco_dir = mk_dist_rslt_dir(dist_rslt_dir,"ref_co");

			int *shuffled_refname_ind = shuffleN( ref_stat->infile_num, 0 );
			mem_usage_stat.input_file_name_sz = ref_stat->infile_num * ( sizeof(llong) + PATHLEN*sizeof(char) );
real_time_mem -=  mem_usage_stat.input_file_name_sz;			

			dim_shuffle = get_dim_shuffle(opt_val);
//			printf("shufid=%d,k=%d,L=%d\n",dim_shuffle->dim_shuffle_stat.id,dim_shuffle->dim_shuffle_stat.k,
//				dim_shuffle->dim_shuffle_stat.drlevel);
  		hashsize = get_hashsz(dim_shuffle);
  		seq2co_global_var_initial();
  		mem_usage_stat.shuffled_subctx_arr_sz = (1LLU << 4*( dim_shuffle->dim_shuffle_stat.subk) )*sizeof(int);
real_time_mem -= mem_usage_stat.shuffled_subctx_arr_sz;
			//estimate threads num for stageI
			p_fit_mem = real_time_mem / ( (hashsize + 1) * sizeof(llong));
			if( opt_val->p < p_fit_mem )
				p_fit_mem = opt_val->p; //opt_val->p == 1 if !def _OPENMP
			else if(p_fit_mem < 1)
				err(errno,"dist_dispatch():\n"
				" Kmer hashing need mem.(%lf G) exceed the mem. system or user provide (%lf G)\n"
				" user can either consider specify more mem.(-m ) or use smaller k value ( -k)\n" 
				" or increase dimension reduction level ( -L)\n", 
				(double)hashsize * sizeof(llong)/1e9, opt_val->mmry);
		
			size_t est_pthreadcoMem4stageII = (llong)p_fit_mem * (hashsize + 1) * sizeof(llong);
real_time_mem -= est_pthreadcoMem4stageII ;			
      fprintf(logfp, "COMPONENT_SZ=%d\tcompnum=%d\t.fas:%d\t.fq:%d\tall:%d\nthreadnum:%d\tp hash mem. %lu\n",
        COMPONENT_SZ,component_num,ref_fmt_count->fasta,ref_fmt_count->fastq,ref_stat->infile_num,p_fit_mem,est_pthreadcoMem4stageII);
			if(ref_fmt_count->fastq >0 )
     		fprintf(logfp,"quality filter=%d\tKmer Occrence filter=%d\n",opt_val->kmerqlty,opt_val->kmerocrs);

			//if (real_time_mem > est_unitllmco_mem()) then build .mco in mem.	

			//run stage I
			const char *refcostat = run_stageI(opt_val, ref_stat, shuffled_refname_ind, dist_refco_dir, p_fit_mem);
			//estimate threads num for stageII
			size_t est_mem4unitllmco = precise_est_unitllmco_mem(refcostat); 
			p_fit_mem = real_time_mem / est_mem4unitllmco ;		
			fprintf(logfp,"\nest unit llmco mem.= %lu\tavailable mem.= %lld\t p_fit_mem = %d\n\n",
        est_mem4unitllmco, real_time_mem, p_fit_mem);
			
			if( opt_val->p < p_fit_mem ) p_fit_mem = opt_val->p;	
			else if(p_fit_mem < 1) err(errno,"dist_dispatch():\n" 
			" creat .mco db file need mem.(%lf G) exceed the mem. system or user provide (%lf G)\n"  
			" user might specify more mem. ( use -m ) or increase dimension reduction level ( use -L)\n",
			 est_mem4unitllmco/1e9, opt_val->mmry);

			run_stageII(refcostat,p_fit_mem);			
			free(dim_shuffle->shuffled_dim);
real_time_mem += mem_usage_stat.shuffled_subctx_arr_sz;
		}
		else if ( refco_dstat_fpath != NULL ){

			size_t est_mem4unitllmco = precise_est_unitllmco_mem(refco_dstat_fpath);
			p_fit_mem = real_time_mem / est_mem4unitllmco ;
			if( opt_val->p < p_fit_mem ) p_fit_mem = opt_val->p;
			else if(p_fit_mem < 1) err(errno,"dist_dispatch():creat .mco db file need mem.(%lf G)\n"
      " exceed the mem. system or user provide (%lf G) user might specify more mem. ( use -m )\n",
       est_mem4unitllmco/1e9, opt_val->mmry);

			run_stageII(refco_dstat_fpath,p_fit_mem);
		}
		else if(refmco_dstat_fpath != NULL) {;} //read mco..

		free((char*)refco_dstat_fpath);
		free((char*)refmco_dstat_fpath);
	} 
	// if has query
	if ( (opt_val->num_remaining_args >0) || (opt_val->fpath[0] != '\0' ) ) 
	{		//mem_dispatch.has_arg = 1;
		const char *qryco_dstat_fpath = NULL;
		const char *qrymco_dstat_fpath	= NULL;
		if(opt_val->num_remaining_args >0){
			qryco_dstat_fpath = test_get_fullpath(opt_val->remaining_args[0],co_dstat);
			qrymco_dstat_fpath = test_get_fullpath(opt_val->remaining_args[0],mco_dstat);
		}			
	
		// database serach mode; 
		if(opt_val->refpath[0] != '\0'){
			const char *ref_db = test_get_fullpath(opt_val->refpath, mco_dstat);
			const char *distdir = mk_dist_rslt_dir(opt_val->outdir, distoutdir);
			
			if(ref_db==NULL) err(errno,"need speficy the mco path for -r to run the query-ref search model");
			FILE * ref_mco_stat_fp;
		  if ((ref_mco_stat_fp = fopen(ref_db,"rb")) == NULL) err(errno,"mco stat file:%s",ref_db );
			mco_dstat_t mco_ref_dstat;
			fread( &mco_ref_dstat, sizeof(mco_dstat_t),1,ref_mco_stat_fp );
			fclose(ref_mco_stat_fp);				
			//search by .co file in opt_val->remaining_args[0]
			if( qryco_dstat_fpath != NULL ){
				FILE *qry_co_stat_fp;
				if (( qry_co_stat_fp = fopen(qryco_dstat_fpath,"rb")) == NULL) err(errno,"qry co stat file:%s",qryco_dstat_fpath); 					
				co_dstat_t co_qry_dstat;
				fread( &co_qry_dstat, sizeof(co_dstat_t), 1, qry_co_stat_fp);	
				if(co_qry_dstat.shuf_id != mco_ref_dstat.shuf_id ) 
					err(errno, "qry shuf_id: %d not match ref shuf_id: %d\ntry regenerate .co dir and feed -s the .shuf"
											"file used to generated ref database",co_qry_dstat.shuf_id,mco_ref_dstat.shuf_id);
				else if(co_qry_dstat.comp_num != mco_ref_dstat.comp_num)
					err(errno, "qry comp_num: %d not match ref comp_num: %d",co_qry_dstat.comp_num, mco_ref_dstat.comp_num);

				//initiital global distance output parameter before distance compute
				neighbor_num = opt_val->num_neigb ;
				mutdistance_thre = opt_val->mut_dist_max ;			
				mco_co_dist(opt_val->refpath, opt_val->remaining_args[0], distdir, opt_val->p);
		  }
			else if (qrymco_dstat_fpath != NULL)
				err(errno,"when -r specified, the query sould not be .mco format, the valid query format shoulde be .fas/.fq file or .co");
			else{ // test seq file format and search by .fas .fq 
				infile_tab_t *infile_stat = dist_organize_infiles(opt_val);
				infile_fmt_count_t* qry_fmt_count;
    		qry_fmt_count = infile_fmt_count(infile_stat);
				bool is_valid_fas_fq_in = (infile_stat->infile_num != 0) &&
    			(qry_fmt_count->fasta + qry_fmt_count->fastq == infile_stat->infile_num );
				if(is_valid_fas_fq_in){	
				/* fas/fq search db code */
				}
				else err(errno,"please specify valid query genomes seq or .co file for database search");		
			}
		}
		// if only .co query 
		else if( qryco_dstat_fpath != NULL) { // convert .co folder opt_val->remaining_args[0]  to .mco
     	size_t est_mem4unitllmco = precise_est_unitllmco_mem(qryco_dstat_fpath);
			p_fit_mem = real_time_mem / est_mem4unitllmco ;

			if( opt_val->p < p_fit_mem ) p_fit_mem = opt_val->p;
      else if(p_fit_mem < 1) err(errno,"dist_dispatch():\n"
      " creat .mco db file need mem.(%lf G) exceed the mem. system or user provide (%lf G)\n"
      " user might specify more mem. ( use -m )\n",
       est_mem4unitllmco/1e9, opt_val->mmry);

      fprintf(logfp,"\nest unit llmco mem.= %lu\tavailable mem.= %lld\t p_fit_mem = %d\n",
        est_mem4unitllmco, real_time_mem, p_fit_mem);

      run_stageII(qryco_dstat_fpath, p_fit_mem);
    }
		else if ( (qrymco_dstat_fpath != NULL) && (opt_val->num_remaining_args > 1 )){
			//run .mco combine mode
		}
		else {	//if only query => test if is valid raw seq format				
				infile_tab_t *infile_stat = dist_organize_infiles(opt_val);
        infile_fmt_count_t* qry_fmt_count;
        qry_fmt_count = infile_fmt_count(infile_stat);
        bool is_valid_fas_fq_in = (infile_stat->infile_num != 0) &&
          (qry_fmt_count->fasta + qry_fmt_count->fastq == infile_stat->infile_num );
				// if is valid raw seq format		
        if(is_valid_fas_fq_in){ //convert .fas .fq to co

      		const char * dist_rslt_dir = mk_dist_rslt_dir(opt_val->outdir,"kssd_dist_rslt");
      		const char *dist_co_dir = mk_dist_rslt_dir(dist_rslt_dir,"co");
					
					int *shuffled_refname_ind = shuffleN( infile_stat->infile_num, 0 );
					mem_usage_stat.input_file_name_sz = infile_stat->infile_num * ( sizeof(llong) + PATHLEN*sizeof(char) );
real_time_mem -=  mem_usage_stat.input_file_name_sz;

					dim_shuffle = get_dim_shuffle(opt_val);
      		hashsize = get_hashsz(dim_shuffle);
      		seq2co_global_var_initial();
      		mem_usage_stat.shuffled_subctx_arr_sz = (1LLU << 4*( dim_shuffle->dim_shuffle_stat.subk) )*sizeof(int);
real_time_mem -= mem_usage_stat.shuffled_subctx_arr_sz;

					p_fit_mem = real_time_mem / ( (hashsize + 1) * sizeof(llong));    
					if( opt_val->p < p_fit_mem ) p_fit_mem = opt_val->p;
      		else if(p_fit_mem < 1)
        		err(errno,"dist_dispatch():\n"
        		" Kmer hashing need mem.(%lf G) exceed the mem. system or user provide (%lf G)\n"
        		" user can either consider specify more mem.(-m ) or use smaller k value ( -k)\n"
        		" or increase dimension reduction level ( -L)\n",
        		(double)hashsize * sizeof(llong)/1e9, opt_val->mmry);

					size_t est_pthreadcoMem4stageII = (llong)p_fit_mem * (hashsize + 1) * sizeof(llong);
real_time_mem -= est_pthreadcoMem4stageII ;

					fprintf(logfp,".fas:%d\t.fq:%d\tall:%d\nthreadnum:%d\tp hash mem. %lu\n",
        		qry_fmt_count->fasta,qry_fmt_count->fastq, infile_stat->infile_num,p_fit_mem,est_pthreadcoMem4stageII);
					if(qry_fmt_count->fastq >0 )
						fprintf(logfp,"quality filter=%d\tKmer Occrence filter=%d\n",opt_val->kmerqlty,opt_val->kmerocrs);	
					
      		run_stageI(opt_val,infile_stat,shuffled_refname_ind,dist_co_dir,p_fit_mem);
        }
				else err(errno,"not valid raw seq format");
		}// only query and is raw seq format	
	}// if has query end

	fclose(logfp);
	return 0;
}; // dispatch() end

dim_shuffle_t *get_dim_shuffle( dist_opt_val_t *opt_val_in )
{
	char shuf_infile_name_prefix[PATHLEN+9];
  char shuf_infile_name[PATHLEN];
  strcpy(shuf_infile_name, opt_val_in->dr_file);

  if( strcmp( shuf_infile_name, "" ) == 0 )
  {
    fprintf(logfp,"addLen=%d\tsubctxlen=%d\n",add_len_drlevel2subk(),add_len_drlevel2subk() + opt_val_in->dr_level);
    srand ( time(NULL) );
    dim_shuffle_stat_t dim_shuffle_stat =
    {
      rand(), //id
      opt_val_in->k, //k
      opt_val_in->dr_level + add_len_drlevel2subk(), //subk, add len compute according to MIN_SUBCTX_DIM_SMP_SZ
      opt_val_in->dr_level,  //drlevel
    };

		struct stat outd;
		if( (stat(opt_val_in->outdir, &outd) != 0 ) || (! S_ISDIR(outd.st_mode) ) ) 
			mkdir(opt_val_in->outdir,0777);

		sprintf(shuf_infile_name_prefix, "%s/default",opt_val_in->outdir);	
    write_dim_shuffle_file( &dim_shuffle_stat,shuf_infile_name_prefix); //"default"
		fprintf(logfp,"subcontext shuffled dimension file: %s.shuf created\n",shuf_infile_name_prefix);
    sprintf(shuf_infile_name,"%s.shuf",shuf_infile_name_prefix);
  };
  return read_dim_shuffle_file(shuf_infile_name);
};

int get_hashsz(dim_shuffle_t *dim_shuffle_in )
{
	int dim_reduce_rate = 1 << 4*dim_shuffle_in->dim_shuffle_stat.drlevel;
	llong ctx_space_sz = 1LLU << 4*( dim_shuffle_in->dim_shuffle_stat.k - dim_shuffle_in->dim_shuffle_stat.drlevel );
  int primer_ind = 4*( dim_shuffle_in->dim_shuffle_stat.k - dim_shuffle_in->dim_shuffle_stat.drlevel ) - CTX_SPC_USE_L - 7;
  // check primer_ind
  if(primer_ind < 0 || primer_ind > 24 ){
    int k_add = primer_ind < 0 ?  (1 + (0-primer_ind)/4)  :  - (1 + ( primer_ind - 24 )/4)  ;
    err(errno,"get_hashsz(): primer_ind: %d out of range(0 ~ 24), by formula:\n"
    "int primer_ind = 4*(opt_val->k - dim_shuffle->dim_shuffle_stat.drlevel) - CTX_SPC_USE_L - 7\n"
    "this might caused by too small or too large k\n"
    "kmer length = %d\n"
    "dim reduction level = %d\n"
    "ctx_space size = %llu\n"
    "CTX space usage limit = %lf\n\n"
    "try rerun the program with option -k = %d",
    primer_ind, dim_shuffle_in->dim_shuffle_stat.k, dim_shuffle_in->dim_shuffle_stat.drlevel,ctx_space_sz,
       (double)1/(1 << CTX_SPC_USE_L),dim_shuffle_in->dim_shuffle_stat.k + k_add );
  };
  int hashsize_get = primer[primer_ind];
  fprintf(logfp,"dimension reduced %d\n"
          "ctx_space size=%llu\n"
          "k=%d\n"
          "drlevel=%d\n"
          "primer_ind=%d\n"
          "hashsize=%u\n",
    dim_reduce_rate,ctx_space_sz,dim_shuffle_in->dim_shuffle_stat.k, dim_shuffle_in->dim_shuffle_stat.drlevel, primer_ind, hashsize_get);

	return hashsize_get ;
}

//test if parent_path contain file dstat_f, return full path if has otherwise return NULL
//orginally develop fo test dstat_f , but can generally test other file 
const char* test_get_fullpath(const char *parent_path, const char *dstat_f)
{
	struct stat path_stat;
  if( stat(parent_path, &path_stat) < 0 ) 
		err(errno,"test_get_fullpath()::%s",parent_path);
  if( S_ISDIR(path_stat.st_mode) ){
		char* fullpath = malloc(PATHLEN+1);
		sprintf((char *)fullpath,"%s/%s", parent_path, dstat_f);
		FILE *fp;
		if ( (fp = fopen(fullpath,"rb")) != NULL ){
			fclose(fp);
			return fullpath;
		}
		else{
			free((char*)fullpath);
			return NULL;
		}
	}
	else
		return NULL;
};

const char * run_stageI (dist_opt_val_t *opt_val, infile_tab_t *seqfile_stat,
					 int* shuffled_seqfname_ind, const char *co_dir, int p_fit_mem)
{
      // get bin num
      int ref_full_bin_num = seqfile_stat->infile_num / BIN_SZ ;
      int binsz ;
			
			// seq2co by bin 
      llong **CO = malloc( p_fit_mem * sizeof(llong *) );
      for(int i = 0; i< p_fit_mem; i++ )
        CO[i] = (llong *)malloc(hashsize * sizeof(llong) );

			llong all_ctx_ct = 0 ;
			ctx_obj_ct_t *ctx_ct_list = malloc(sizeof(ctx_obj_ct_t) * seqfile_stat->infile_num);			

      for(int i = 0; i<= ref_full_bin_num; i++)
      {
        if( i == ref_full_bin_num )
        {
          binsz =  seqfile_stat->infile_num % BIN_SZ;
          if(binsz == 0) break; //if no remainder, for loop should stop when i == ref_full_bin_num
        }
        else binsz = BIN_SZ;
#pragma omp parallel for num_threads(p_fit_mem) reduction(+:all_ctx_ct) schedule(guided)
        for(int j=0;j< binsz; j++)
        {
          int tid = 0;
#ifdef _OPENMP
          tid = omp_get_thread_num();
#endif
          char* seqfname = seqfile_stat->organized_infile_tab[ shuffled_seqfname_ind[i*BIN_SZ + j ] ].fpath;
					char cofname[PATHLEN];
					//path/binnum.gid.co.componentid
					sprintf(cofname,"%s/%d.%d.co", co_dir,i, j );
					printf("decomposing %s\n",seqfname) ;
					llong *co;
					if(isOK_fmt_infile(seqfname,fastq_fmt,FQ_FMT_SZ)){
						co = fastq2co(seqfname,CO[tid],opt_val->kmerqlty,opt_val->kmerocrs);
				  	ctx_ct_list[i*BIN_SZ + j] = write_fqco2file(cofname,co);
						all_ctx_ct += ctx_ct_list[i*BIN_SZ + j] ;
					}else{
						co = fasta2co(seqfname,CO[tid]);
						ctx_ct_list[i*BIN_SZ + j] = wrt_co2cmpn_use_inn_subctx(cofname,co);			
						all_ctx_ct += ctx_ct_list[i*BIN_SZ + j] ;
					}						
        };
      }// end all bin loop
      // free **CO
      for(int i = 0; i< p_fit_mem; i++ )
        free(CO[i]);
      free(CO);

      //write co_stat file
      co_dstat_t co_dstat_wrout;
			co_dstat_wrout.shuf_id = dim_shuffle->dim_shuffle_stat.id ; 
			co_dstat_wrout.kmerlen = dim_shuffle->dim_shuffle_stat.k * 2;
			co_dstat_wrout.dim_rd_len = dim_shuffle->dim_shuffle_stat.drlevel * 2 ;
      co_dstat_wrout.comp_num = component_num ;
      co_dstat_wrout.infile_num = seqfile_stat->infile_num;
			co_dstat_wrout.all_ctx_ct = all_ctx_ct; 

      char *co_dstat_fullname = malloc(PATHLEN*sizeof(char) );
      sprintf(co_dstat_fullname, "%s/%s",co_dir,co_dstat);

      FILE *coutfp;
      if ( ( coutfp = fopen(co_dstat_fullname,"wb")) == NULL ) err(errno,"%s",co_dstat_fullname);
      fwrite(&co_dstat_wrout,sizeof(co_dstat_wrout),1,coutfp);
			//write file ctx_ct and names			
			fwrite(ctx_ct_list,sizeof(ctx_obj_ct_t),co_dstat_wrout.infile_num,coutfp );
			free(ctx_ct_list);

			for(int i = 0; i< co_dstat_wrout.infile_num; i++)
				fwrite(seqfile_stat->organized_infile_tab[ shuffled_seqfname_ind[i] ].fpath,PATHLEN,1,coutfp);				
      fclose(coutfp);			
			return (const char *)co_dstat_fullname;
}

void run_stageII(const char * co_dstat_fpath, int p_fit_mem)
{
			//get full path of co stat and mco stat file						
		  const char* dist_co_dir = get_pathname(co_dstat_fpath,co_dstat);
		  const char* dist_rslt_dir = malloc(PATHLEN*sizeof(char));		
			sprintf((char*)dist_rslt_dir,"%s/..",dist_co_dir);
      const char* dist_mco_dir = mk_dist_rslt_dir(dist_rslt_dir,"mco");
      const char* mco_dstat_fpath = malloc(PATHLEN*sizeof(char));
			sprintf((char*)mco_dstat_fpath,"%s/%s",dist_mco_dir,mco_dstat);			

			FILE *co_stat_fp,*mco_stat_fp;
			//read stat file from co dir
			if( ( co_stat_fp = fopen(co_dstat_fpath,"rb")) == NULL ) err(errno,"run_stageII(():%s",co_dstat_fpath);
			co_dstat_t co_dstat_readin;	
			//subtract last pointer size if co_dstat_t last member is pointer
			fread( &co_dstat_readin, sizeof(co_dstat_t),1,co_stat_fp );
			//char (*seqfilebasename)[PATHLEN] = malloc(co_dstat_readin.infile_num * PATHLEN);			
			// write stat file to mco dir
			if( ( mco_stat_fp = fopen(mco_dstat_fpath,"wb")) == NULL ) err(errno,"run_stageII(():%s",mco_dstat_fpath);
			mco_dstat_t mco_dstat_writeout;	
			mco_dstat_writeout.shuf_id = co_dstat_readin.shuf_id;		
			mco_dstat_writeout.kmerlen = co_dstat_readin.kmerlen ;
      mco_dstat_writeout.dim_rd_len = co_dstat_readin.dim_rd_len ;
			mco_dstat_writeout.infile_num = co_dstat_readin.infile_num;
			mco_dstat_writeout.comp_num = co_dstat_readin.comp_num ;
			fwrite(&mco_dstat_writeout,sizeof(mco_dstat_writeout),1, mco_stat_fp );
			//cp ctx_ct list and filname to mco dir			
			unsigned int *tmp_ctx_ct = malloc(sizeof(unsigned int)*co_dstat_readin.infile_num);
			fread(tmp_ctx_ct,sizeof(unsigned int),co_dstat_readin.infile_num,co_stat_fp);		
			fwrite(tmp_ctx_ct,sizeof(unsigned int),co_dstat_readin.infile_num,mco_stat_fp);
			char (*tmpname)[PATHLEN] = malloc( PATHLEN * co_dstat_readin.infile_num );			
			fread(tmpname,PATHLEN,co_dstat_readin.infile_num,co_stat_fp);
			fwrite(tmpname,PATHLEN,co_dstat_readin.infile_num,mco_stat_fp);

			free(tmp_ctx_ct);
			free(tmpname);
			fclose(co_stat_fp);
			fclose(mco_stat_fp);
			// unit .mco binning parametors 	
			int last_binsz = co_dstat_readin.infile_num % BIN_SZ ;			
			int full_bin_num = co_dstat_readin.infile_num / BIN_SZ ;
			int bin_num = last_binsz == 0 ? full_bin_num : full_bin_num + 1 ;

			int unit_arrmco_ct = bin_num * co_dstat_readin.comp_num ;				
    	printf("p_fit_mem=%d\tBIN_SZ=%d\treadin component_num=%d\torigin component_num=%d\tbinnum=%d\tlast_binsz=%d\tunit_arrmco_ct=%d\tfnum=%d\n", p_fit_mem,BIN_SZ,co_dstat_readin.comp_num,component_num,bin_num,last_binsz,unit_arrmco_ct,co_dstat_readin.infile_num );
			fprintf(logfp,"p_fit_mem=%d\tBIN_SZ=%d\treadin component_num=%d\torigin component_num=%d\tbinnum=%d\tlast_binsz=%d\tunit_arrmco_ct=%d\tfnum=%d\n", p_fit_mem,BIN_SZ,co_dstat_readin.comp_num,component_num,bin_num,last_binsz,unit_arrmco_ct,co_dstat_readin.infile_num );	
			//make unit arrmco and write to file
#pragma omp parallel for  num_threads(p_fit_mem) schedule(guided) 
      for ( int n = 0; n < unit_arrmco_ct ; n++ )
      {
				int binsz;
				if( n >= full_bin_num*co_dstat_readin.comp_num ){	
					binsz = last_binsz;
				}else binsz = BIN_SZ;

				int tidm  = omp_get_thread_num();
				printf("co2unitllmco tid= %d\n",tidm);
				int bin_code = n / co_dstat_readin.comp_num ;
				int comp_code = n % co_dstat_readin.comp_num ;				

        mco_entry_stat_t** llmco = co2unitllmco(dist_co_dir, binsz, bin_code, comp_code);
        gidobj_t** arrmco = llmco2arrmco(llmco);
        char arrmcofname[PATHLEN];
        sprintf(arrmcofname,"%s/%d.mco.%d",dist_mco_dir,bin_code, comp_code);
        int write_lines_count =  write_unit_arrmco_file(arrmcofname,arrmco);
				fprintf(logfp,"arrmco file %s write out %d lines\n",arrmcofname,write_lines_count );
				printf("tid=%d\tarrmco file %s write out %d lines\n",tidm,arrmcofname,write_lines_count );
				free_unit_arrmco(arrmco);
      }
			
			free((char*)dist_co_dir );
			free((char*)dist_rslt_dir);
			free((char*)dist_mco_dir);
			free((char*)mco_dstat_fpath);
}

static ctx_obj_ct_t initial_dist[BIN_SZ];
static int ref_seq_num,qry_seq_num,kmerlen,dim_reduct_len;
void mco_co_dist( char *refmco_dname, char *qryco_dname, const char *distout_dir, int p_fit_mem)
{
	fprintf(logfp,"run mco_co_dist(), %d threads used\n",p_fit_mem);
	printf("run mco_co_dist(), %d threads used\n",p_fit_mem);
	FILE *refmco_dstat_fp, *qryco_dstat_fp;
	char *refmco_dstat_fpath = malloc(PATHLEN*sizeof(char));	
	char *qryco_dstat_fpath = malloc(PATHLEN*sizeof(char));
	sprintf(refmco_dstat_fpath,"%s/%s",refmco_dname,mco_dstat);
	sprintf(qryco_dstat_fpath,"%s/%s",qryco_dname,co_dstat); 

	if( (refmco_dstat_fp = fopen(refmco_dstat_fpath,"rb")) == NULL ) 
		err(errno,"need provied mco dir path for mco_co_dist() arg 1. refmco_dstat_fpath");
	if( (qryco_dstat_fp = fopen(qryco_dstat_fpath,"rb")) == NULL ) 
		err(errno,"need provied co dir path for mco_co_dist() arg 2.  qryco_dstat_fpath");
		
	mco_dstat_t mco_dstat_readin ;
  co_dstat_t co_dstat_readin ;
 	fread(&mco_dstat_readin,sizeof(mco_dstat_readin),1,refmco_dstat_fp);
	fread(&co_dstat_readin,sizeof(co_dstat_readin),1,qryco_dstat_fp);

	//read ctx_ct_list
	unsigned int * qry_ctx_ct_list = malloc(co_dstat_readin.infile_num * sizeof(unsigned int));
	unsigned int * ref_ctx_ct_list = malloc(mco_dstat_readin.infile_num * sizeof(unsigned int)); 
	fread(qry_ctx_ct_list,sizeof(unsigned int),co_dstat_readin.infile_num,qryco_dstat_fp);
	fread(ref_ctx_ct_list,sizeof(unsigned int),mco_dstat_readin.infile_num,refmco_dstat_fp);
	// read filenames	
	char (*cofname)[PATHLEN] = malloc(co_dstat_readin.infile_num * PATHLEN);
	char (*mcofname)[PATHLEN] = malloc(mco_dstat_readin.infile_num * PATHLEN);	
	fread(cofname,PATHLEN,co_dstat_readin.infile_num,qryco_dstat_fp);
	fread(mcofname,PATHLEN,mco_dstat_readin.infile_num,refmco_dstat_fp);	
  
	fclose(refmco_dstat_fp);
  fclose(qryco_dstat_fp);

	if( !(mco_dstat_readin.comp_num == co_dstat_readin.comp_num) )
    err(errno,"query args not match ref args: ref.comp_num = %d vs. %d = qry.comp_num",
    mco_dstat_readin.comp_num, co_dstat_readin.comp_num);			
	if(!(mco_dstat_readin.shuf_id == co_dstat_readin.shuf_id))
		err(errno,"query args not match ref args: ref.shuf_id = %d vs. %d = qry.shuf_id",
    mco_dstat_readin.shuf_id, co_dstat_readin.shuf_id);	
	
	int ref_bin_num = mco_dstat_readin.infile_num / BIN_SZ;
	int binsz; //comp_sz = (1 << 4*COMPONENT_SZ) ;	
  // set when use mco_co_dist_core() 
	//mco_co_dist_t *shared_ctx_num = malloc( sizeof(mco_co_dist_t) * p_fit_mem );
	//mco_co_dist_t *diff_obj_num = malloc( sizeof(mco_co_dist_t) * p_fit_mem );
	for(int i=0; i<=ref_bin_num;i++ ){
		if( i == ref_bin_num ){
    	binsz =  mco_dstat_readin.infile_num % BIN_SZ;
          if(binsz == 0) continue;
    }else binsz = BIN_SZ;
		//create dist file
#pragma omp parallel for  num_threads(p_fit_mem) schedule(guided)
		for( int k=0; k< co_dstat_readin.infile_num; k++){

			if(qry_ctx_ct_list[k]==0){
				warnx("%dth co file is empty",k);
				continue;
			}
			char dist_fcode[PATHLEN];
			sprintf(dist_fcode,"%s/%d.%d.dist", distout_dir, i, k );

			FILE *distfp;
			if ( (distfp = fopen(dist_fcode,"wb")) == NULL) err(errno,"mco_co_dist()::%s",dist_fcode);
			fwrite(initial_dist, sizeof(ctx_obj_ct_t), binsz, distfp);				
			fclose(distfp);
		}
		fprintf(logfp,"all %d co files' distance file initialized\n",co_dstat_readin.infile_num);		

		for ( int j = 0; j < mco_dstat_readin.comp_num; j++ ){

			char mco_fcode[PATHLEN];
			sprintf(mco_fcode,"%s/%d.mco.%d",refmco_dname,i,j);					
			gidobj_t** unit_arrmco_readin = read_unit_arrmco_file(mco_fcode);

#pragma omp parallel for  num_threads(p_fit_mem) schedule(guided)
			for(int k=0; k< co_dstat_readin.infile_num; k++){

				if(qry_ctx_ct_list[k]==0) continue;

        int tid = 0;
#ifdef _OPENMP
        tid = omp_get_thread_num();
#endif
				char co_fcode[PATHLEN]; char dist_fcode[PATHLEN];
				sprintf(co_fcode,"%s/%d.%d.co.%d",qryco_dname,k/BIN_SZ, k % BIN_SZ, j);
				sprintf(dist_fcode,"%s/%d.%d.dist",distout_dir, i , k );

				int fd;
				if( ( (fd = open(dist_fcode,O_RDWR, 0600) ) == -1) )
         err(errno,"mco_co_dist()::distfile = %s[tid = %d]",dist_fcode,tid);
				
				ctx_obj_ct_t * ctx_obj_ct = mmap(NULL, binsz*(sizeof(ctx_obj_ct_t)), PROT_READ | PROT_WRITE, MAP_SHARED, fd, 0); 		
				if(ctx_obj_ct == MAP_FAILED) err(errno,"ctx_obj_ct mmap error");				
				
				ctx_obj_ct = mco_co_mmpdist_core(unit_arrmco_readin,co_fcode,ctx_obj_ct);
				//mco_co_dist_core(unit_arrmco_readin,co_fcode,binsz,shared_ctx_num[0],diff_obj_num[0]);
				if ( msync( ctx_obj_ct, binsz*(sizeof(ctx_obj_ct_t)), MS_SYNC ) < 0 ) 
					err(errno,"mco_co_dist()::ctx_obj_ct msync failed"); 

				munmap(ctx_obj_ct,binsz*(sizeof(ctx_obj_ct_t)));
				close(fd);
			}
			free_unit_arrmco(unit_arrmco_readin);
		}
	}
	//=== distance output======// 
	/* global var iniital */
	ref_seq_num = mco_dstat_readin.infile_num ;
  qry_seq_num = co_dstat_readin.infile_num ;
	kmerlen = co_dstat_readin.kmerlen;
  dim_reduct_len = co_dstat_readin.dim_rd_len;
	char distf[PATHLEN];
	sprintf(distf, "%s/distance.out", distout_dir);
	fprintf(logfp,"distance output to : %s\n",distf);
	printf("distance output to : %s\n",distf);
	FILE *distfp;
	if( (distfp = fopen(distf,"a")) == NULL ) err(errno,"mco_co_dist():%s",distf);
	for(int i=0; i<=ref_bin_num;i++ ){
		for ( int k = 0; k < co_dstat_readin.infile_num; k++ ){
			if(qry_ctx_ct_list[k]>0)
				fname_dist_print(i,k,distout_dir,ref_ctx_ct_list,qry_ctx_ct_list,mcofname,cofname,distfp);
		}
	}
	fclose(distfp);

	free(ref_ctx_ct_list);
	free(qry_ctx_ct_list);
	free(mcofname);
	free(cofname);
}//func end


static inline ctx_obj_ct_t * mco_co_mmpdist_core(gidobj_t** unit_arrmco, char *co_fcode_in, ctx_obj_ct_t * ctx_obj_ct_in )
{
  mmp_uint_t mmpcofile;
  mmpcofile = mmp_uint_arr(co_fcode_in);

  unsigned int ind,mcogid;
  int ctx_num = mmpcofile.fsize/sizeof(unsigned int);
  
	for(int n = 0; n < ctx_num; n++){
    ind = mmpcofile.mmpco[n];
    if(unit_arrmco[ind] != NULL){
      for(unsigned int k = 1; k< unit_arrmco[ind][0] + 1; k++ ){
        mcogid = unit_arrmco[ind][k] ;
        ctx_obj_ct_in[mcogid]++ ; //ctx_obj_ct_t [0] is ctx_ct, [1] is obj_c
     }
    }
  }
  munmap(mmpcofile.mmpco, mmpcofile.fsize);
	return ctx_obj_ct_in;
}

static inline void mco_co_dist_core(gidobj_t** unit_arrmco, char *co_fcode_in, int bin_sz, 
						mco_co_dist_t shared_ctx_num_in )
{
	mmp_uint_t mmpcofile;
	mmpcofile = mmp_uint_arr(co_fcode_in);

	memset(shared_ctx_num_in,0, bin_sz);

	unsigned int ind,mcogid; 
	int ctx_num = mmpcofile.fsize/sizeof(unsigned int);	
	for(int n = 0; n < ctx_num; n++){	
		ind = mmpcofile.mmpco[n] ;		
		if(unit_arrmco[ind] != NULL){
			for(int k = 1; k< unit_arrmco[ind][0] + 1; k++ ){
				mcogid = unit_arrmco[ind][k];			
			  shared_ctx_num_in[mcogid]++ ;
			} 			
		}	
	}
	munmap(mmpcofile.mmpco, mmpcofile.fsize);
	//printf dist
	printf("%s:",co_fcode_in);
	for(int e=0; e < bin_sz; e++ ){
	//	printf ("%f\t", (double)diff_obj_num_in[e]/shared_ctx_num_in[e]);
	}
	printf("\n");	
}

void dist_print( const char *distf, FILE *dist_fp )
{
	ctx_obj_ct_t *ctx_obj_ct; 
	int fd;
	struct stat s;
	fd = open (distf, O_RDONLY);
	check (fd < 0, "open %s failed: %s", distf, strerror (errno));
	fstat (fd, & s);

	ctx_obj_ct = mmap(0, s.st_size , PROT_READ, MAP_PRIVATE, fd, 0);
	check ( ctx_obj_ct == MAP_FAILED, "mmap %s failed: %s", distf, strerror (errno));

  fprintf(dist_fp,"output %s\n",distf);
	for(int i = 0;i < s.st_size/sizeof(ctx_obj_ct_t); i++)
	//	fprintf(dist_fp,"%s\t%d\t%d\t%d\t%f\n",distf,i,ctx_obj_ct[i],(double)ctx_obj_ct[i][1]/ctx_obj_ct[i][0]);	
	close(fd);
	munmap(ctx_obj_ct, s.st_size);
}

// move in function if do multiple threads distance output 
char full_distfcode[PATHLEN];
//int ctx_len = 16; // 2*k, depend on k ,but might not affect if mashD is propotional to real mashD 
void fname_dist_print(int ref_bin_code, int qry_fcode, const char *distout_dir, unsigned int*ref_ctx_ct_list, 
			unsigned int*qry_ctx_ct_list, char (*refname)[PATHLEN], char (*qryfname)[PATHLEN], FILE *dout_fp)
{
	sprintf(full_distfcode,"%s/%d.%d.dist",distout_dir,ref_bin_code,qry_fcode);
  ctx_obj_ct_t *ctx_obj_ct;	 
	int fd;
	struct stat s;
	fd = open (full_distfcode, O_RDONLY);
	check (fd < 0, "open %s failed: %s", full_distfcode, strerror (errno));
	fstat (fd, & s);

	ctx_obj_ct = mmap(0, s.st_size , PROT_READ, MAP_PRIVATE, fd, 0);
	check ( ctx_obj_ct == MAP_FAILED, "mmap %s failed: %s", full_distfcode, strerror (errno));
//	fprintf(dout_fp,"output %s\n",full_distfcode);	
	double jac_ind,contain_ind,Dm,Da,P_K_in_X_XnY, P_K_in_Y_XnY,
	j_prim, c_prim, Dm_prim, Da_prim, sd_j_prim, sd_c_prim,
	CI95_j_prim1,	CI95_j_prim2, CI95_c_prim1, CI95_c_prim2,
	CI95_Dm_prim1,CI95_Dm_prim2,CI95_Da_prim1,CI95_Da_prim2;
	int Min_XY_size, X_size, Y_size, XnY_size, XuY_size, X_XnY_size, Y_XnY_size ;
	int alp_size = 4; // dna 	

	for(int i = 0;i < s.st_size/sizeof(ctx_obj_ct_t); i++) {

		X_size = ref_ctx_ct_list[ref_bin_code*BIN_SZ + i];
		Y_size = qry_ctx_ct_list[qry_fcode];
		Min_XY_size = X_size < Y_size ? X_size : Y_size ;
		XnY_size = ctx_obj_ct[i];
		XuY_size =  X_size + Y_size - XnY_size ;
		X_XnY_size = X_size - XnY_size;
		Y_XnY_size = Y_size- XnY_size;		

		jac_ind = (double)XnY_size / XuY_size;		
		contain_ind = (double)XnY_size / Min_XY_size ;
		Dm = jac_ind ==1? 0: -log(2*jac_ind/(1+jac_ind)) / kmerlen ;
		Da = contain_ind==1? 0: -log(contain_ind) / kmerlen ;
				
		P_K_in_X_XnY = 1 - pow( (1- 1/pow(alp_size,(kmerlen - dim_reduct_len) )), X_XnY_size );
		P_K_in_Y_XnY = 1 - pow( (1- 1/pow(alp_size,(kmerlen - dim_reduct_len) )), Y_XnY_size );
		double rs = P_K_in_X_XnY * P_K_in_Y_XnY * ( X_XnY_size + Y_XnY_size )
							/(P_K_in_X_XnY + P_K_in_Y_XnY - 2*P_K_in_X_XnY * P_K_in_Y_XnY);
		j_prim = ((double)XnY_size - rs) / XuY_size ;
		c_prim = ((double)XnY_size - rs) / Min_XY_size ;
		Dm_prim = j_prim == 1? 0:-log(2*j_prim/(1+j_prim)) / kmerlen ;
		Da_prim = c_prim==1? 0:-log(c_prim) / kmerlen ;
		
		sd_j_prim = pow(j_prim*(1 - j_prim) / XuY_size, 0.5) ;
		sd_c_prim = pow(c_prim*(1 - c_prim) / Min_XY_size,0.5) ;
		CI95_j_prim1 = j_prim - 1.96*sd_j_prim;
		CI95_j_prim2 = j_prim + 1.96*sd_j_prim;
		CI95_c_prim1 = c_prim - 1.96*sd_c_prim;
		CI95_c_prim2 = c_prim + 1.96*sd_c_prim;
		CI95_Dm_prim1 = CI95_j_prim2 == 1? 0:-log(2*CI95_j_prim2/(1+CI95_j_prim2)) / kmerlen ;
		CI95_Dm_prim2 = CI95_j_prim1 == 1? 0:-log(2*CI95_j_prim1/(1+CI95_j_prim1)) / kmerlen ;
		CI95_Da_prim1 = CI95_c_prim2 == 1? 0:-log(CI95_c_prim2) / kmerlen ; 
		CI95_Da_prim2 = CI95_c_prim1 == 1? 0:-log(CI95_c_prim1) / kmerlen ;
			
		double pv_j_prim = 0.5 * erfc( j_prim / sd_j_prim * pow(0.5,0.5)  );
		double pv_c_prim = 0.5 * erfc( c_prim/sd_c_prim * pow(0.5,0.5)  );
		double qv_j_prim = pv_j_prim * ref_seq_num*qry_seq_num ;
		double qv_c_prim = pv_c_prim * ref_seq_num*qry_seq_num ;

	//	if(Dm_prim < mutdistance_thre)
			fprintf(dout_fp,"%s\t%s\t%u-%u|%u|%u\t%lf\t%lf\t%lf\t%lf\t%lf[%lf,%lf]\t%lf[%lf,%lf]\t%lf[%lf,%lf]\t%lf[%lf,%lf]\t%E\t%E\t%E\t%E\n", 
				qryfname[qry_fcode],refname[ref_bin_code*BIN_SZ + i],XnY_size,(unsigned int)rs,X_size,Y_size,jac_ind,Dm,
        contain_ind,Da,j_prim,CI95_j_prim1,CI95_j_prim2,Dm_prim,CI95_Dm_prim1,CI95_Dm_prim2,c_prim,CI95_c_prim1,
				CI95_c_prim2,Da_prim,CI95_Da_prim1,CI95_Da_prim2,pv_j_prim,pv_c_prim,qv_j_prim,qv_c_prim);			
	}
	close(fd);
  munmap(ctx_obj_ct, s.st_size);
}
// collect and orgnize seq files names from --list,
// remaining args(including folder or files) into a string array
infile_tab_t* dist_organize_infiles (dist_opt_val_t *opt_val) 
{
	//do it if fpath is not ""
	if( strcmp( opt_val->fpath, "" ) !=0 ) 
	{
		return organize_infile_list(opt_val->fpath);
	}
	else if( opt_val->num_remaining_args > 0 )
	{
		return organize_infile_frm_arg(opt_val->num_remaining_args, opt_val->remaining_args);
	}
	else
	{
		perror("please specify the input/query files");
    return NULL;
	}
};

extern const char *acpt_infile_fmt[ACPT_FMT_SZ] ;

infile_tab_t* dist_organize_refpath( dist_opt_val_t *opt_val){
	struct stat path_stat;
	if( stat(opt_val->refpath, &path_stat) < 0) 
		err(errno,"dist_organize_refpath():%s",opt_val->refpath); 
	if( S_ISDIR(path_stat.st_mode) || isOK_fmt_infile(opt_val->refpath, acpt_infile_fmt,ACPT_FMT_SZ) ){
		// char **tmp_arg = malloc(sizeof(char **)); 
		char * tmp_arg[] = {opt_val->refpath};		
		return organize_infile_frm_arg(1, tmp_arg);
	}
	else if(S_ISREG(path_stat.st_mode))
		return organize_infile_list(opt_val->refpath);
	else
		return NULL;		
}




