#!/bin/bash

while true
do
WDIR='/home/daq/runs/temp_test_ashiph_sipm/'
CMD='/home/daq/dist/wavedump-3.10.0/src/wavedump'
CFGfile="$WDIR/WaveDumpConfig.txt"
tmpfile='delme.txt'
dir="$WDIR/`date +"%Y-%m-%d_%H-%M-%S"`"
mkdir -v $dir
cd $dir
cp -vp $CFGfile $dir
#
ln -s /dev/null wave_1.dat
ln -s /dev/null wave_3.dat
ln -s /dev/null wave_5.dat
ln -s /dev/null wave_6.dat
ln -s /dev/null wave_7.dat
#
../vkbd.sh | $CMD $CFGfile
chmod 444 TR_0_0.dat wave_0.dat wave_2.dat wave_4.dat
rm -v wave_[13567].dat PlotData.txt
cd ../

if [[ $(du $dir | cut -f1) -lt 50000 ]];
then
        echo 'WARNING! File size of '${dir}' is small' | mutt -s "Temperature ASHIPH-SiPM tests" "i.v.ovtin@inp.nsk.su"
fi

sleep 600
done
