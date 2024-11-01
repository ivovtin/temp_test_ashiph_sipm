#!/bin/bash

indirs=()

count=0
while read str
do
        count=$[ count + 1 ]
        indirs+=($str)
        echo "Line contents are : $str "
done < runs_list.dat
echo "count: $count "

list=( 0 2 4 )
#list=( 2 4 )

for (( i = 0; i < $count; i++ ))
do
	for n in ${list[@]}
	do
		root -l -b -q 'histo.C("'${indirs[i]}'",'${n}')'
	done
done
