#ifndef COMMAND_COMPONENT
#define COMMAND_COMPONENT
#include <argp.h>
#include "global_basic.h"

typedef struct component_opt
{
	int p; //threads
	char refdir[PATHLEN]; 
	char qrydir[PATHLEN];
	char outdir[PATHLEN];
} component_opt_t; //dimension shuffle type


int cmd_component(struct argp_state* state);
int comparator_idx (const void*, const void*);

int comparator (const void*, const void*);

int get_species_abundance (component_opt_t *);

#endif 
