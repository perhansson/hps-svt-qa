#!/bin/bash
if [ -z "$1" ]
then
	echo "cal.sh <timing data dir> <baseline data file> <output dir> <output filename>"
	exit
fi
pushd "${0/%`basename $0`/}"
binpath="${PWD/%root/}/bin"
popd

mkdir -p $3
$binpath/meeg_baseline $2 -n -t -o $3/$4
$binpath/meeg_baseline $2 -mn -o $3/$4_mux

$binpath/meeg_tp $1/*cal_?.bin -frn -o $3/$4
mkdir -p $3/fits
mv $3/$4_tp_fit* $3/fits

$binpath/meeg_t0res $3/$4.base $3/$4 $1/*cal_?.bin -n -o $3/$4_anfit

mkdir -p $3/neg
mv $3/*_neg.png $3/neg
