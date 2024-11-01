#!/bin/bash

while read -a line; 
do 
	#echo -e "${line[1]}";
        if [[ ${line[1]} > 55.44 ]]
	then
		echo -e "${line[1]}";
  		echo 'Caution - danger! Temperature exceeded.' | mutt -s "Temperature ASHIPH-SiPM tests" "i.v.ovtin@inp.nsk.su"
	fi
done < res.dat
