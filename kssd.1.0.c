#include "global_wrapper.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "global_basic.h"
  /** Global Options **/
const char *argp_program_version =  "kssd version 1.0";
const char *argp_program_bug_address = "yhg926@gmail.com";
//argp_err_exit_status = 1;

int main(int argc, char** argv)
{
	
  setvbuf(stdout, (char *)NULL, _IOLBF, 0);
  setvbuf(stderr, (char *)NULL, _IOLBF, 0);
  // initial domain and long domain(with subcommand)
	const char subcommand_n[] = "<subcommand>";
  domain = (char*) malloc(strlen(argv[0])+1);
  strcpy(domain,argv[0]);
  long_domain = (char*) malloc( strlen(domain) + strlen(subcommand_n) + 1);
  snprintf(long_domain,strlen(domain) + strlen(subcommand_n) + 2,"%s %s",domain,subcommand_n);

  cmd_global(argc, argv);
  return 0;
}
