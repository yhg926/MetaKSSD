#include "co2mco.h"
#include "command_dist.h"
#include <err.h>
#include <errno.h>
#include <math.h>
#ifdef _OPENMP
  #include <omp.h>
#endif
typedef struct kmerdb_index
{
	size_t *row_offset;
	unsigned int *row_gnum;
	unsigned int *row_bin_gnum;
} kmerdb_index_t;
void cdb_kmerf2kmerdb(const char *mcodirname, const char *codirname, int cofnum, int comp_num, int p_fit_mem)
{	
	kmerdb_index_t mco_map;		
	unsigned int comp_sz = (1 << 4*COMPONENT_SZ) ;
	int binnum = ceil((double)cofnum/BIN_SZ);	
	// gnum in each bin each row
	mco_map.row_bin_gnum = malloc(comp_sz*binnum*sizeof(unsigned int));
	mco_map.row_gnum = malloc(comp_sz * sizeof(unsigned int));
	mco_map.row_offset = malloc(comp_sz*sizeof(size_t) );	
	mco_map.row_offset[0] = 0;

	//make sure fsize of combco.index.i = sizeof(size_t)*(cofnum + 1);
	size_t *cbdcoindex = malloc( sizeof(size_t)*(cofnum + 1) );
	char cbdcofname[PATHLEN]; char cbdcoindexf[PATHLEN];
	char mcofname[PATHLEN]; char mcoindexf[PATHLEN];
	mmp_uint_t mmpcbd_cofile;

	gid_arr_llist_t** mco = malloc(comp_sz* sizeof(gid_arr_llist_t*));
	for(unsigned int i = 0; i< comp_num; i++){
		//reset count to 0 for each component
		memset(mco_map.row_bin_gnum,0,comp_sz*binnum*sizeof(unsigned int) );
		memset(mco_map.row_gnum,0,comp_sz * sizeof(unsigned int));

		FILE *cbdfp, *cbdindexfp;
		sprintf(cbdcoindexf,"%s/combco.index.%d",codirname,i);		
		sprintf(cbdcofname,"%s/combco.%d",codirname,i);
    if( (cbdfp = fopen(cbdcofname,"rb")) == NULL) err(errno,"%s",cbdcofname);
    if( (cbdindexfp = fopen(cbdcoindexf,"rb")) == NULL) err(errno,"%s",cbdcoindexf);	
		//get index
		fread(cbdcoindex,sizeof(size_t),cofnum + 1,cbdindexfp);			
		mmpcbd_cofile = mmp_uint_arr(cbdcofname);
				
/********cbdco to llmco block**********/
		for(int j = 0;j< cofnum; j++ ){
#pragma omp parallel for num_threads(p_fit_mem) schedule(guided)
			for(size_t k = cbdcoindex[j]; k< cbdcoindex[j+1]; k++){
				unsigned int ind = mmpcbd_cofile.mmpco[k];	
				unsigned int mod = mco_map.row_gnum[ind] % GID_ARR_SZ ;
				gid_arr_llist_t* tmp;
				if(mod == 0){
					tmp = mco[ind];
					mco[ind] = malloc(sizeof(gid_arr_llist_t));
					if (mco[ind] == NULL) err(errno,"cdb_kmerf2kmerdb()::mco[ind]") ;
					mco[ind]->next = tmp; 
				}								
				mco[ind]->gidobj[mod] = j % BIN_SZ ;
				mco_map.row_gnum[ind]++;
				mco_map.row_bin_gnum[ind*binnum + j/BIN_SZ]++; //#define arr[x][y] arr[x*binnum+y/BIN_SZ], x:row_num, y:bin_num
			} 
		}
		munmap(mmpcbd_cofile.mmpco, mmpcbd_cofile.fsize);

		for(int n=1; n<comp_sz; n++)
			mco_map.row_offset[n] = mco_map.row_offset[n-1] + mco_map.row_gnum[n-1];
	
/**************write index file**************/
		sprintf(mcoindexf,"%s/mco.index.%d",mcodirname,i);
		FILE *arrmco_index_fp = fopen(mcoindexf,"wb");		
		if( arrmco_index_fp  == NULL) err(errno,"%s",mcoindexf);
		fwrite(mco_map.row_offset,sizeof(size_t),comp_sz,arrmco_index_fp);
		fwrite(mco_map.row_bin_gnum,sizeof(unsigned int),comp_sz*binnum, arrmco_index_fp);
		fclose(arrmco_index_fp);	
/**************write llmco to arrmco file**************/
		sprintf(mcofname,"%s/mco.%d",mcodirname,i);
		int arrmco_fp = open(mcofname,O_RDWR|O_CREAT,0600) ;
		if (arrmco_fp == -1) err(errno,"cdb_kmerf2kmerdb()::%s",mcofname);
		size_t mco_comp_fsize = sizeof(gidobj_t)*(mco_map.row_offset[comp_sz-1] + mco_map.row_gnum[comp_sz-1]) ;		
		if(ftruncate(arrmco_fp,mco_comp_fsize) == -1) err(errno,"cdb_kmerf2kmerdb()::ftruncate");
		// mcofile mmaped
		gidobj_t*	mco_mmpf = mmap(NULL,mco_comp_fsize,PROT_WRITE,MAP_SHARED,arrmco_fp,0);	
		close(arrmco_fp);
#pragma omp parallel for  num_threads(p_fit_mem) schedule(guided)
		for(unsigned int s = 0; s< comp_sz ; s++ ){
			if( mco_map.row_gnum[s] == 0 ) continue; //arr_len
				gid_arr_llist_t *tmpblk; // for swap
				//get end address of s row
				gidobj_t* current_blkpos_mapin_arrmco = mco_mmpf + mco_map.row_offset[s] + mco_map.row_gnum[s] ;
				
				int blk_num = mco_map.row_gnum[s]/GID_ARR_SZ;
				int remainder = mco_map.row_gnum[s] % GID_ARR_SZ;				
				if( remainder > 0 ) blk_num+=1;
				int blk_len;

				for(int blk = 0; blk < blk_num; blk++){

					if((blk==0) && (remainder > 0) ) blk_len = remainder;
					else blk_len = GID_ARR_SZ;

					current_blkpos_mapin_arrmco -= blk_len;
					memcpy(current_blkpos_mapin_arrmco, mco[s]->gidobj, blk_len * sizeof(gidobj_t));

					tmpblk = mco[s];
					mco[s] = mco[s]->next;
					free(tmpblk);

				}//end blk						
		} // end all rows
		if ( msync( mco_mmpf, mco_comp_fsize, MS_ASYNC ) < 0 )
          err(errno,"cdb_kmerf2kmerdb()::msync failed");		
	 	munmap(mco_mmpf,mco_comp_fsize);			
		fclose(cbdfp);
		fclose(cbdindexfp);
	}  //componets loop end

  free(mco_map.row_bin_gnum); 
  free(mco_map.row_gnum); 
  free(mco_map.row_offset);
  free(cbdcoindex);	
}

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
	fclose(co_stat_fp);	
	return mem_sz;
}



