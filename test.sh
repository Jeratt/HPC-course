#!/bin/bash

prog="a.out"
ans_file="out.txt"
tmp_file=".tmp_ans"

truncate -s 0 $tmp_file
while IFS=' ' read -r c1 c2 c3 c4 ; do
  `./$prog $c1 $c2 $c3 $c4 >> $tmp_file`
  echo "" >> $tmp_file 
done <in.txt

cmp -b $tmp_file $ans_file
