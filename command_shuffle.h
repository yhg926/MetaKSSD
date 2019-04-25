#ifndef COMMAND_SHUFFLE
#define COMMAND_SHUFFLE

typedef struct dim_shuffle_stat
{
	int id; //random shuffle file id 
	int k; //half context len
	int subk; //half suncontex len
	int drlevel; // dimension reduction level 	
} dim_shuffle_stat_t; //dimension shuffle type

typedef struct dim_shuffle
{
	dim_shuffle_stat_t dim_shuffle_stat;
	int *shuffled_dim;
} dim_shuffle_t;


//minimal allowed subcontext
#define MIN_SUBCTX_DIM_SMP_SZ 4096 //256 

int * shuffle ( int arr[], int len_arr);
int * shuffleN ( int n, int base);//return shuffled $base based natural number array
int write_dim_shuffle_file(dim_shuffle_stat_t* dim_shuffle_stat, char *outfile_prefix);
dim_shuffle_t* read_dim_shuffle_file(char *dim_shuffle_file);

int add_len_drlevel2subk(void);
	#include <argp.h>
	int cmd_shuffle(struct argp_state* state);
#endif 
