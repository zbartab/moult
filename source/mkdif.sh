#!/bin/bash

# process the dump_repval files produced by moult program

awk '{if($1!=52) print $0}' repval_dump1.txt > rv.tmp
awk '{print $7}' rv.tmp > rv1.txt
awk '{print $1, $2, $3, $4, $5, $6}' rv.tmp > states.txt
awk '{print $8, $9, $10, $11}' rv.tmp > pol1.txt

awk '{if($1!=52) print $0}' repval_dump2.txt > rv.tmp
awk '{print $7}' rv.tmp > rv2.txt
awk '{print $8, $9, $10, $11}' rv.tmp > pol2.txt

paste rv1.txt rv2.txt | awk '{print $2-$1}' > rv_diff.txt
