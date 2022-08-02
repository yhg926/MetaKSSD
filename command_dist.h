#ifndef COMMAND_DIST
#define COMMAND_DIST

#include "global_basic.h"
#include "command_dist_wrapper.h"
#include "command_shuffle.h"
#include <string.h>

// memory management
#define DISM_MEM_PCT 0.25 //largest memory percent dismtrix can occupy
#define CO_MEM_PCT 0.125 //largest memory percent .co can occupy
extern const char co_dstat[];
typedef struct mem_dispatch
{
  bool has_ref;
  bool has_arg;
  bool dm_inmem; //determine wether distance matrix in memory.
  bool keep_coindisk;  //determine wether keep intermedia .co file in file.
  bool keep_mcoinmem; //determine wether keep .mco in mem when make .co file.
  bool alphabet_group; //wether group objs by alphabet.
  bool multiple_composition;

  enum{
    OBJ_COMPATIBLE = 4096, //obj compatible mode, 4 bit for obj,12 bits for gid
    OBJ_TRADITION = 16384, //obj incompatible,no alphabet group, 2bit obj, 14 bits gid
    OBJ_GROUP = 65536, // objs grouped by alphabet, 0bits obj,16 bits gid
  } enum_bin_sz ;

  enum{
     SLIM, //use small self define binsize due to memory constraints
     STD,  //use BIN_SZ
     ALL,  //no bin and all genomes in, for reference seqs only, may use int size genome id.
  } gbin  ;

  enum {
    MAKE_REF, //make and output ref .mco only, no distance compute and query search ( -r no -l/args)
    REF_QRY, //ouput top n, or distance < x ref seqs for each input query (need both -r and -l/args)
    DIST, //output distance matrix/tab distance (no -r, use -l/args, use -t if tab distance )
  } usage_mode;

  bool mco_info; //add info column for mco file

} mem_dispatch_t ;


typedef struct mem_usage_stat
{
  llong shuffled_subctx_arr_sz;
  llong input_file_name_sz;
  llong others; // stack usage, threads ...
} mem_usage_stat_t;

typedef long long int llint;
/* move to global.h bcs command_set use it
typedef struct co_dirstat
{	
	unsigned int shuf_id;
	bool koc; //kmer occurence or not
	int kmerlen; //2*k 
	int dim_rd_len; //2*drlevel
  int comp_num; // components number
  int infile_num; // .co file num, .co file with many components only count as 1
  //int comp_sz[comp_num * infile_num]; // do this if storage all .co per mco bin in one file
	llong all_ctx_ct;	
} co_dstat_t;
*/
typedef struct mco_dirstat
{
	unsigned int shuf_id;
	//bool koc; //kmer occurence or not
	int kmerlen; //2*k
	int dim_rd_len; //2*drlevel
  int comp_num; // components number
  int infile_num; // use to compute last bin size
} mco_dstat_t;

typedef struct output 
{
	int X_size;
	int Y_size;
	int XnY_size;
	double metric;
	double dist;
 	double CI95_mtrc_1;  // use prime
	double CI95_mtrc_2;
	double CI95_dist_1;
	double CI95_dist_2;
	double pv;          // use prime
} output_t;

typedef struct print_ctrl
{
	//real print control
	MTRIC metric;  //metric // make sure are first two element in enum { Jcd, Ctm, Bth } MTRIC (in command_dist_wrap.h);
	PFIELD pfield; // print out filed
	bool correction; // use rs correction or not, default not
	double dthreshold ; // distance threshold, make sure default initialized to 1

	// some unchanged variable during each qury to all reference
  unsigned int Y_size; //query size
  llong cmprsn_num; //total comparison, for qvalue snprintf
  char *qname ; //queryname
	int qry_len; //qryname string length
}	print_ctrl_t;


//for use in command_decomopse and iseqroco.c only, not for shuffle.c, do not include in shuffle.c
extern dim_shuffle_t* dim_shuffle;  
extern unsigned int hashsize;
extern int component_num;

int dist_dispatch(struct dist_opt_val* opt_val);
infile_tab_t* dist_organize_infiles (dist_opt_val_t *opt_val);
infile_tab_t* dist_organize_refpath (dist_opt_val_t *opt_val);
//const char* get_co_dstat_fpath(const char *refpath); 
//const char* test_get_fullpath(const char *parent_path, const char *dstat_f);
dim_shuffle_t *get_dim_shuffle( dist_opt_val_t *opt_val_in );
int get_hashsz(dim_shuffle_t *dim_shuffle_in );
const char * run_stageI (dist_opt_val_t *opt_val,infile_tab_t *seqfile_stat, 
				int* shuffled_seqfname_ind, const char *co_dir, int p_fit_mem);

void run_stageII(const char * co_dstat_fpath, const char* dist_mco_dir, int p_fit_mem);
//void run_stageII(const char * co_dstat_fpath, int p_fit_mem);
void mco_cbdco_nobin_dist(dist_opt_val_t *opt_val_in); //currently used version(20220802) 
void mco_co_dist( char *refmco_dname, char *qryco_dname, const char *distout_dir, int p_fit_mem);
void mco_cbd_co_dist(dist_opt_val_t *opt_val_in);
void mco_cbd_koc_compatible_dist (dist_opt_val_t *opt_val_in);
//void mco_cbd_co_dist ( char *refmco_dname, char *qryco_dname,const char *distout_dir, int p_fit_mem, llong mem_limit );
void dist_print( const char *distf, FILE *dist_fp );
void fname_dist_print(int ref_bin_code, int qry_fcode, const char *distout_dir, unsigned int*ref_ctx_ct_list,
      unsigned int*qry_ctx_ct_list, char (*refname)[PATHLEN], char (*qryfname)[PATHLEN], FILE *dout_fp);

void dist_print_nobin ( const char *distout_dir,unsigned int ref_num, unsigned int qry_num, unsigned int*ref_ctx_ct_list,
      unsigned int*qry_ctx_ct_list, int num_cof_batch, char (*refname)[PATHLEN], char (*qryfname)[PATHLEN],dist_opt_val_t *opt_val);

void koc_dist_print_nobin ( const char *distout_dir,unsigned int ref_num, unsigned int qry_num, unsigned int*ref_ctx_ct_list,
      unsigned int*qry_ctx_ct_list, int num_cof_batch, char (*refname)[PATHLEN], char (*qryfname)[PATHLEN]);

const char * combine_queries(dist_opt_val_t *opt_val);
#endif 






