#!/bin/bash
for i in `ls *.c`; do gcc -fpreprocessed -dD -E -P $i > ../public_kssd/$i ; done
for i in `ls *.h`; do gcc -fpreprocessed -dD -E -P $i > ../public_kssd/$i ; done
cp Makefile ../public_kssd/ 
