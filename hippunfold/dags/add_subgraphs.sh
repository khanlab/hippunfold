#!/bin/bash


script=./dag_subgraphs.py

#pipe stdin into a temp file
tmp_dag=`mktemp`
cat > $tmp_dag

#print header lines
grep -v '^	' $tmp_dag | head -n -1

#print subgraphs
cat $tmp_dag | python $script 

#print edges
grep '^	' $tmp_dag

#print footer (ie closing brace)
grep -v '^	' $tmp_dag | tail -n 1

#echo $part1 $part2 $part3  $part4

rm $tmp_dag
