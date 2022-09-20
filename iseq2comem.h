#ifndef SEQ2CO_H
#define SEQ2CO_H
#include "global_basic.h"
#define READSEQ_BUFFSZ 65536 
//for koc file //kmer occrence
#define OCCRC_BIT 16 //24 // make sure it is 64 bits machine
#define OCCRC_MAX 0xffffLLU // 0xffffffLLU //must be 1<<OCCRC_BIT - 1

extern void seq2co_global_var_initial(void);

//llong * bgzfasta2co(char* seqfname,char seqin_buff[READSEQ_BUFFSZ], llong *co);
llong * mmpfasta2co(char* seqfname, llong *co);
llong * fasta2co(char* seqfname,llong *co,char * pipecmd);
llong * uniq_fasta2co(char* seqfname,llong *co,char * pipecmd);
llong * fastq2co(char* seqfname, llong *co, char * pipecmd,int Q, int M );
llong * fastq2koc (char* seqfname, llong *co, char * pipecmd, int Q);
llong * mt_shortreads2koc (char* seqfname, llong *co, char *pipecmd,int p);

llong write_fqco2file(char* cofilename, llong *co);
llong wrt_co2cmpn_use_inn_subctx(char* cofilename, llong *co);
llong writeco2file(char* cofilename, llong *co);
llong write_fqkoc2file(char* cofilename, llong *co);
unsigned int write_fqkoc2files(char* cofilename, llong *co);

int reads2mco(char* seqfname,const char *co_dir, char * pipecmd);

#endif
