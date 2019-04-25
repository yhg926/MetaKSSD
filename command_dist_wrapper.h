#ifndef  DIST_WRAPPER
#define DIST_WRAPPER 

#include "global_basic.h"
#include <stdbool.h>
#include <argp.h>

int cmd_dist(struct argp_state* state);
const char *mk_dist_rslt_dir (const char *parentdirpath, const char * outdirpath );

typedef struct dist_opt_val
{
  int k ; //halfctx len
  int p ; //threads counts
	int dr_level;//dimension reduction level 
	char dr_file[PATHLEN];//dimension reduction file
  double mmry; //maxMemory;
  char fmt[10]; 
  char refpath[PATHLEN]; //reference sequences path
  char fpath[PATHLEN]; //query files path 
  char outdir[PATHLEN]; // results dir
	//fastq
  int kmerocrs;
  int kmerqlty;
	bool keepco; 
	bool stage2; //input is intermedia .co or not
	int num_neigb; // distance output filter: num of nearest neight 
	double mut_dist_max;// distance output filter: maximun allow mutation distance
	bool abundance;
  int num_remaining_args;
  char **remaining_args;
  //char outp[128]= "./";
} dist_opt_val_t;

#endif 
