#ifndef COMMAND_COMPONENT
#define COMMAND_COMPONENT
#include <argp.h>
#include "global_basic.h"

extern const char binVec_suffix[];
extern const char abunMtx_suffix[];
extern const char binVec_dirname[];
extern const char abunMtx_idx_suffix[];
extern const char abunMtx_name_suffix[];
extern const char y_l2n_suffix[];
typedef struct binVec
{
	int ref_idx; //reference species order idx
	float pct;
} binVec_t;

typedef struct composite_opt
{
	int b; //write out abundance binary vector (1) or not (0)
	int i; // index abundance binary vectors (1) or not (0)
	int s;
	int d; // read abv or not
	int p; //threads
	char refdir[PATHLEN]; 
	char qrydir[PATHLEN];
	char outdir[PATHLEN];
	int num_remaining_args;
  char **remaining_args;
} composite_opt_t; //dimension shuffle type


int cmd_composite(struct argp_state* state);
int comparator_idx (const void*, const void*);
int comparator_measure (const void *a, const void *b);
int comparator (const void*, const void*);

int get_species_abundance (composite_opt_t *);
int index_abv (composite_opt_t *);
int abv_search  (composite_opt_t *composite_opt);
int read_abv (composite_opt_t *composite_opt);

#endif 
