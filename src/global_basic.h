#ifndef GLOBAL_BASIC
#define GLOBAL_BASIC
#include <stdbool.h>
/*CTX/OBJ/Kmer maskers*/
#define _64MASK 0xffffffffffffffffLLU
//#define TUPMASK (_64MASK >> (64-BITTL))
#define BIT1MASK 0x0000000000000001LLU
/*compile mode for different alphabet: 1:nt reduction,2:AA*/
#if ALPHABET == 1 // objs compatible mode
	#define DEFAULT 15
	#define OBJ_ALPH 16
	#define OBJ_BITS 4
	#define BIN_SZ 4096 // maximum sequences allowed in one bin 
//	#include <stdbool.h>
  extern const bool Objdist[16][16];
  #define IS_DIFF(X,Y) ( Objdist[(X)][(Y)] )

#elif ALPHABET == 2 //AA seq
	#define DEFAULT (-1)
	#define OBJ_BITS 5
	#define BIN_SZ 2048
	#define IS_DIFF(X,Y) ( (X) != (Y) )
#else
	#define DEFAULT (-1)
	#define OBJ_ALPH 4
	#define OBJ_BITS 2
	#define BIN_SZ 65536 //50 //for test //16384
	#define IS_DIFF(X,Y) ( (X) != (Y) )
#endif
/**/
#define LMAX 4096 // characters limit per line for any documents might be analysis 
#define PATHLEN 256 // limit for input file path length
#define MCO_BUF_S 4096 // buf size of reading/analysis mco 
	//COMPONENT_SZ control .co file divided how many components 
#ifndef COMPONENT_SZ
#define COMPONENT_SZ 8 //6 or 7,default unit component dimension size = 16^(COMPONENT_SZ) or 1<<4*(COMPONENT_SZ)
#endif

#ifndef CTX_SPC_USE_L
#define CTX_SPC_USE_L 8 // 8 for small mem.(< 1g), 4 for larger mem.  //ctx space occupy rate limit = 1/(1<<CTX_SPC_USE_L)  
#endif

#define CTX_DR_LMT 100 //limit for least CTX after dimensionality reduction
#define LD_FCTR 0.6 //hash function load factor

//*****core data struct *****//
//gid array linked list: .mco genome id+obj array size
#define GID_ARR_SZ 16 //32 //short arr[GID_ARR_SZ] 
#define BBILLION 1073741824  //binary billion
/*type define */
typedef unsigned long long int llong;
/*argp wrapper*/
	#include <argp.h>
	/** Local Prototypes **/
	 struct arg_global {	int verbosity; };
	
	void log_printf(struct arg_global* g, int level, const char* fmt, ...);
	#define ARGP_KEY_INVALID 16777219	


//const char co_dstat[] = "cofiles.stat"; //stat file name in co dir
//const char mco_dstat[] = "mcofiles.stat";

//*****basic function prototypes ****/

FILE * fpathopen (const char *dpath, const char *fname, const char *mode );
double get_sys_mmry(void);
/*completeReverse*/
#define SWAP2  0x3333333333333333ULL
#define SWAP4  0x0F0F0F0F0F0F0F0FULL
#define SWAP8  0x00FF00FF00FF00FFULL
#define SWAP16 0x0000FFFF0000FFFFULL
#define SWAP32 0x00000000FFFFFFFFULL
static inline llong crvs64bits(llong n) {

  n = ((n >> 2 ) & SWAP2 ) | ((n & SWAP2 ) << 2 );
  n = ((n >> 4 ) & SWAP4 ) | ((n & SWAP4 ) << 4 );
  n = ((n >> 8 ) & SWAP8 ) | ((n & SWAP8 ) << 8 );
  n = ((n >> 16) & SWAP16) | ((n & SWAP16) << 16);
  n = ((n >> 32) & SWAP32) | ((n & SWAP32) << 32);
  return ~n;
}

/*basemap for different alphabet */
extern const int Basemap[128];
extern const char Mapbase[];
extern const unsigned int primer[25];
llong find_lgst_primer_2pow(int w);
int nextPrime(int);
/* orgnized infile table*/

typedef struct infile_entry { 
	size_t fsize; 
	char* fpath;
} infile_entry_t ;

typedef struct infile_tab {
  int infile_num;
  infile_entry_t* organized_infile_tab; // infile_entry_t arr
} infile_tab_t ;

#define BASENAME_LEN 128 //largest allowed bytes for input file basename
typedef struct bin_stat	{
	//estimated kmer count sum accross files in the bin (dimension reduction rate not comsidered here)
	llong est_kmc_bf_dr;
	char (*seqfilebasename)[BASENAME_LEN];
	llong AllcoMem; //exact co file size in mem. estimate after stage I
} bin_stat_t;


infile_tab_t * organize_infile_list(char* list_path,int fmt_ck);
infile_tab_t * organize_infile_frm_arg (int num_remaining_args, char ** remaining_args,int fmt_ck);
//per bin get file basename and estimate .co files (before dimension reduction) memory usage.
bin_stat_t * get_bin_basename_stat(infile_entry_t* organized_infile_tab, int *shuffle_arr,int binsz);

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

//********** input file formats test ********************//

/*input file formats array size */
#define ACPT_FMT_SZ 7
#define FAS_FMT_SZ 4
#define FQ_FMT_SZ	2
#define CO_FMT_SZ 1
#define MCO_FMT_SZ 1
#define CMPRESS_FMT_SZ 2


/* input file formats */
extern const char
*acpt_infile_fmt[ACPT_FMT_SZ],
*fasta_fmt[FAS_FMT_SZ],
*fastq_fmt[FQ_FMT_SZ],
*co_fmt[CO_FMT_SZ],
*mco_fmt[MCO_FMT_SZ],
*compress_fmt[CMPRESS_FMT_SZ];

/* format test fun.*/
#include <string.h>
static inline int isCompressfile(char *fname)
{
	int ret = 0;
	for(int i=0; i < CMPRESS_FMT_SZ;i++)
	{
		int basename_len = strlen(fname) - strlen(compress_fmt[i]);
		if( strcmp( ( fname + basename_len ),compress_fmt[i]) == 0 )
			return 1;
	}
	return ret;
}

static inline int isOK_fmt_infile (char* fname, const char *test_fmt[], int test_fmt_arr_sz)
{
  int ret = 0 ;
  char suftmp[10];

	for(int i=0; i < CMPRESS_FMT_SZ;i++ ){
  	int basename_len = strlen(fname) - strlen(compress_fmt[i]);
  	// if infile suffix with .gz or other compress format
  	if( strcmp( ( fname + basename_len ),compress_fmt[i]) == 0 ){
    	char cp_fname[PATHLEN];
    	strcpy(cp_fname, fname);
    	*(cp_fname + basename_len) = '\0';
    	fname = cp_fname;
			break;
  	};
	};

  for(int i=0; i< test_fmt_arr_sz; i++){
    sprintf(suftmp,".%s",test_fmt[i]);
    if ( strcmp((char *)(fname+strlen(fname) - strlen(suftmp)), suftmp) == 0  ){
      return 1;
    }
  };
  return ret;
};


#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <sys/mman.h>
#include <fcntl.h> // open function
#include <unistd.h> //close function

static inline void check (int test, const char * message, ...)
{
    if (test) {
        va_list args;
        va_start (args, message);
        vfprintf (stderr, message, args);
        va_end (args);
        fprintf (stderr, "\n");
        exit (EXIT_FAILURE);
    }
}; 

typedef struct mmpco // member type match writeco2file()
{
	size_t fsize;
	unsigned int *mmpco;
} mmp_uint_t;

static inline mmp_uint_t mmp_uint_arr (char *cofname) 
{
	mmp_uint_t cofilemmp;//char* fpath; llong fsize
	int fd;
	struct stat s;
	fd = open (cofname, O_RDONLY);
	check (fd < 0, "open %s failed: %s", cofname, strerror (errno));
	fstat (fd, & s);
	cofilemmp.fsize = s.st_size;
	cofilemmp.mmpco = mmap(NULL, s.st_size , PROT_READ, MAP_PRIVATE, fd, 0);	
	check ( cofilemmp.mmpco == MAP_FAILED, "mmap %s failed: %s", cofname, strerror (errno));
	close(fd);
	return cofilemmp;
};

typedef struct mmpany
{
  size_t fsize;
  void *mmp;
} mmp_any_t;

static inline mmp_any_t mmp_any (char *fname)
{
  mmp_any_t filemmp;
  int fd;
  struct stat s;
  fd = open (fname, O_RDONLY);
  check (fd < 0, "open %s failed: %s", fname, strerror (errno));
  fstat (fd, & s);
  filemmp.fsize = s.st_size;
  filemmp.mmp = mmap(NULL, s.st_size , PROT_READ, MAP_PRIVATE, fd, 0);
  check ( filemmp.mmp == MAP_FAILED, "mmap %s failed: %s", fname, strerror (errno));
  close(fd);
  return filemmp;
};


void replaceChar(char *str, char oldChar, char newChar);

int str_suffix_match(char *str, const char *suf); 
const char * get_pathname(const char *fullpath, const char *suf);
const char* test_get_fullpath(const char *parent_path, const char *dstat_f);
// infile fmt count struct
typedef struct
{
  int fasta;
  int fastq;
  int co;
	int mco;
} infile_fmt_count_t ;
// infile fmt count function 
infile_fmt_count_t *infile_fmt_count ( infile_tab_t * infile_tab );

extern const char co_dstat[];
extern const char skch_prefix[];
extern const char idx_prefix[];
extern const char pan_prefix[];
extern const char uniq_pan_prefix[];

extern const char mco_dstat[];
extern const char mco_gids_prefix[];
extern const char mco_idx_prefix[];

 
typedef unsigned int ctx_obj_ct_t;

#define H1(K,HASH_SZ) ((K)%(HASH_SZ))
#define H2(K,HASH_SZ) ( 1 + (K) % ( (HASH_SZ) - 1 ) )
#define HASH(K,I,HASH_SZ) ( ( H1(K,HASH_SZ) + I * H2(K,HASH_SZ) ) % HASH_SZ )
#define LOG2(X) ((unsigned) (8*sizeof (unsigned long long) - __builtin_clzll((X)) - 1))

#endif


