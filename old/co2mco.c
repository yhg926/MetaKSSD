#include "co2mco.h"
#include "command_dist.h"
#include <err.h>
#include <errno.h>
//linked list mco per bin per componenet

mco_entry_stat_t** co2unitllmco(const char *codirname, int bin_sz, int bin_id, int component_id)
{	 
	unsigned int comp_sz = (1 << 4*COMPONENT_SZ) ;
	mco_entry_stat_t** mco = calloc(comp_sz, sizeof(mco_entry_stat_t*) );	
	gid_arr_llist_t* tmp; // for swap
	char cofname[PATHLEN];
	mmp_uint_t mmpcofile;
	unsigned int ind;
	int mod;
	for(unsigned int i = 0 ; i < bin_sz ; i++)
	{
		sprintf(cofname,"%s/%d.%d.co.%d",codirname,bin_id,i,component_id);
		mmpcofile = mmp_uint_arr(cofname); 
		int ctx_num = mmpcofile.fsize/sizeof(unsigned int); 

		for(int j = 0; j < ctx_num; j++ )
		{ 
			ind = mmpcofile.mmpco[j]; 
			if(mco[ind] == NULL) mco[ind] = calloc( 1, sizeof(mco_entry_stat_t) );
			mod = mco[ind]->g_num % GID_ARR_SZ ;		
			if(mod == 0)
			{				
				if ( (tmp = malloc(sizeof(gid_arr_llist_t)) ) == NULL) err(errno,"co2unitllmco()") ;
				tmp->next = mco[ind]->next;
				mco[ind]->next = tmp ;
			}
			mco[ind]->next->gidobj[mod] = i ;
			mco[ind]->g_num++;
		};	
		munmap(mmpcofile.mmpco, mmpcofile.fsize);
	}
	return mco;
}

//regularize unit llmco(per bin per componenet)
gidobj_t** llmco2arrmco(mco_entry_stat_t** llmco)
{
	unsigned int comp_sz = (1 << 4*COMPONENT_SZ) ;
	gidobj_t** arrmco = calloc(comp_sz, sizeof(gidobj_t*));
	gidobj_t* current_blkpos_mapin_arrmco;
 
	gid_arr_llist_t *tmpblk, *tmpptr;
	
  for(unsigned int i = 0; i< comp_sz ; i++ ){

		if(llmco[i] == NULL) continue; 
		int arr_len = (int)llmco[i]->g_num + 1; //num of gid + each objalphbet count + total count
  	arrmco[i] = malloc( arr_len * sizeof(gidobj_t) );	
		current_blkpos_mapin_arrmco = arrmco[i] + arr_len;
		arrmco[i][0] = llmco[i]->g_num;
//		memcpy(arrmco[i], llmco[i]->obj_stat, ( (1 << OBJ_BITS) + 1) * sizeof(gidobj_t) );	
		int blk_num = arrmco[i][0] % GID_ARR_SZ == 0 ? arrmco[i][0]/GID_ARR_SZ :arrmco[i][0]/GID_ARR_SZ + 1;

		for(int blk = 0; blk < blk_num; blk++){
			if( blk == 0 ){
				tmpblk = llmco[i]->next;			
				if(arrmco[i][0] % GID_ARR_SZ == 0){
					current_blkpos_mapin_arrmco -= GID_ARR_SZ; //make sure not current_blkpos_mapin_arrmco -= GID_ARR_SZ*sizeof(gidobj_t) ?
					memcpy(current_blkpos_mapin_arrmco, tmpblk->gidobj, GID_ARR_SZ * sizeof(gidobj_t));		
				}
				else {
					current_blkpos_mapin_arrmco -= arrmco[i][0] % GID_ARR_SZ ;	
					memcpy(current_blkpos_mapin_arrmco, tmpblk->gidobj, (arrmco[i][0] % GID_ARR_SZ) * sizeof(gidobj_t));			
				}
				free(llmco[i]);
			}			
			else {
				tmpptr = tmpblk;
				tmpblk = tmpblk->next;
				current_blkpos_mapin_arrmco -= GID_ARR_SZ;
        memcpy(current_blkpos_mapin_arrmco, tmpblk->gidobj,GID_ARR_SZ * sizeof(gidobj_t));
				free(tmpptr);
			}								
		}	//end blk			
  } // end all entry 
	free(llmco);	
	return arrmco;
} // end func


unsigned int write_unit_arrmco_file(const char* unitmcofname, gidobj_t** arrmco)
{
  FILE *outf;
  if( (outf = fopen(unitmcofname,"wb") ) == NULL ) err(errno,"write_unit_arrmco_file()");
  unsigned int comp_sz = (1 << 4*COMPONENT_SZ);
	unsigned int validrow=0; 
  for(unsigned int i = 0; i< comp_sz ; i++ ){
    if(arrmco[i] != NULL){
			validrow++;
			fwrite(&i,sizeof(i),1,outf);
			fwrite(arrmco[i], sizeof(gidobj_t), (unsigned int)arrmco[i][0] + 1,  outf);
		} 
  }
	fclose(outf);
	return validrow;
}

gidobj_t** read_unit_arrmco_file(const char *mco_fncode)
{
	FILE *inf;
	if( (inf = fopen(mco_fncode ,"rb") ) == NULL ) err(errno,"read_unit_arrmco_file()");
	unsigned int comp_sz = (1 << 4*COMPONENT_SZ);	
	unsigned int ind ;
		
	gidobj_t** arrmco = calloc(comp_sz, sizeof(gidobj_t*));
  gidobj_t gid_arr_len;
	
	while( fread(&ind,sizeof(ind),1,inf) == 1){ //while( fread(&ind,sizeof(ind),1,inf) == 1 ){
		fread(&gid_arr_len, sizeof(gidobj_t), 1 , inf);		
		arrmco[ind] = malloc(sizeof(gidobj_t)* ( (unsigned int)gid_arr_len + 1));
		arrmco[ind][0] = gid_arr_len;		
		fread( arrmco[ind] + 1, sizeof(gidobj_t), (unsigned int)gid_arr_len , inf);
	}
	fclose(inf);
	return  arrmco ;
} 

void free_unit_arrmco(gidobj_t** unit_arrmco)
{
	unsigned int comp_sz = (1 << 4*COMPONENT_SZ);
	for(unsigned int i = 0; i < comp_sz ; i++){
		if(unit_arrmco[i] != NULL)
			free(unit_arrmco[i]);
	}
	free(unit_arrmco);
}; 

size_t est_unitllmco_mem(void)
{
  size_t mem_sz = 0;
  unsigned int comp_sz = (1U << 4*COMPONENT_SZ) ;
  mem_sz = comp_sz
          *(  sizeof(mco_entry_stat_t*)
            + sizeof(mco_entry_stat_t)
              //asymptotic num. of gid_arr_llist_t block
            + ( (unsigned int)( ( (double) BIN_SZ / ( 1U << CTX_SPC_USE_L ) ) / GID_ARR_SZ ) + 1 )
              * ( sizeof(gidobj_t) * GID_ARR_SZ + sizeof(gid_arr_llist_t *) )
            );
  return mem_sz;
};

size_t precise_est_unitllmco_mem(const char *co_dstat_fpath)
{
	FILE *co_stat_fp;
	if( ( co_stat_fp = fopen(co_dstat_fpath,"rb")) == NULL )
		 err(errno,"precise_est_unitllmco_mem():%s",co_dstat_fpath);

	co_dstat_t co_dstat_readin;
	fread( &co_dstat_readin, sizeof(co_dstat_t),1,co_stat_fp );

	unsigned int comp_sz = (1U << 4*COMPONENT_SZ) ;
	double ctx_spc_use_rate = (double)co_dstat_readin.all_ctx_ct
		/co_dstat_readin.infile_num/co_dstat_readin.comp_num/comp_sz ;
	printf("ctx_spc_use_rate=%lf\n",ctx_spc_use_rate);
	size_t mem_sz = comp_sz
								*(  sizeof(mco_entry_stat_t*)
								+ sizeof(mco_entry_stat_t)
								+ ( (unsigned int)( ( (double) BIN_SZ * ctx_spc_use_rate ) / GID_ARR_SZ ) + 1 )
									* ( sizeof(gidobj_t) * GID_ARR_SZ + sizeof(gid_arr_llist_t *) )
								);
	
	return mem_sz;
}



