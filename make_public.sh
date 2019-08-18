#!/bin/bash
for i in `ls *.c`; do gcc -fpreprocessed -dD -E -P $i > /tmp/$i ; cat COPYRIGHT.txt /tmp/$i > ../public_kssd/$i ; done
for i in `ls *.h`; do gcc -fpreprocessed -dD -E -P $i > /tmp/$i ; cat COPYRIGHT.txt /tmp/$i > ../public_kssd/$i ; done
cp Makefile LICENSE.txt ../public_kssd/ 
