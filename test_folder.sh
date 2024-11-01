#!/bin/bash

#du -h 2024-05-31_13-20-17

if [[ $(du 2024-05-31_13-20-17 | cut -f1) -lt 150000 ]];
then
	echo 'WARNING! File size of $file is small' | mutt -s "Temperature ASHIPH-SiPM tests" "i.v.ovtin@inp.nsk.su"
fi
