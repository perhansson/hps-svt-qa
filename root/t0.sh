#!/bin/bash
if [ -z "$1" ]
then
	echo "tp.sh <timing data dir> <baseline data file> <output dir> <output filename>"
	exit
fi
pushd "${0/%`basename $0`/}"
binpath="${PWD/%root/}/bin"
popd

mkdir -p $3
$binpath/meeg_baseline $2 -n -t -o $3/$4
$binpath/meeg_baseline $2 -mn -o $3/$4_mux

$binpath/meeg_tp $1/*cal?_[1-8]x3_125ns.bin -frn -s20 -o $3/$4
mkdir -p $3/fits
mv $3/$4_tp_fit* $3/fits

$binpath/meeg_t0res $3/$4.base $3/$4 $1/*cal?_1x3_125ns.bin -n -o $3/$4_anfit_1 &
$binpath/meeg_t0res $3/$4.base $3/$4 $1/*cal?_2x3_125ns.bin -n -o $3/$4_anfit_2 &
$binpath/meeg_t0res $3/$4.base $3/$4 $1/*cal?_3x3_125ns.bin -n -o $3/$4_anfit_3 &
$binpath/meeg_t0res $3/$4.base $3/$4 $1/*cal?_4x3_125ns.bin -n -o $3/$4_anfit_4
$binpath/meeg_t0res $3/$4.base $3/$4 $1/*cal?_5x3_125ns.bin -n -o $3/$4_anfit_5 &
$binpath/meeg_t0res $3/$4.base $3/$4 $1/*cal?_6x3_125ns.bin -n -o $3/$4_anfit_6 &
$binpath/meeg_t0res $3/$4.base $3/$4 $1/*cal?_7x3_125ns.bin -n -o $3/$4_anfit_7 &
$binpath/meeg_t0res $3/$4.base $3/$4 $1/*cal?_8x3_125ns.bin -n -o $3/$4_anfit_8
$binpath/meeg_t0res $3/$4.base $3/$4 $1/*cal?_1x3_125ns.bin -n -o $3/$4_linfit_1 -u &
$binpath/meeg_t0res $3/$4.base $3/$4 $1/*cal?_2x3_125ns.bin -n -o $3/$4_linfit_2 -u &
$binpath/meeg_t0res $3/$4.base $3/$4 $1/*cal?_3x3_125ns.bin -n -o $3/$4_linfit_3 -u &
$binpath/meeg_t0res $3/$4.base $3/$4 $1/*cal?_4x3_125ns.bin -n -o $3/$4_linfit_4 -u
$binpath/meeg_t0res $3/$4.base $3/$4 $1/*cal?_5x3_125ns.bin -n -o $3/$4_linfit_5 -u &
$binpath/meeg_t0res $3/$4.base $3/$4 $1/*cal?_6x3_125ns.bin -n -o $3/$4_linfit_6 -u &
$binpath/meeg_t0res $3/$4.base $3/$4 $1/*cal?_7x3_125ns.bin -n -o $3/$4_linfit_7 -u &
$binpath/meeg_t0res $3/$4.base $3/$4 $1/*cal?_8x3_125ns.bin -n -o $3/$4_linfit_8 -u

$binpath/meeg_t0res $3/$4.base $3/$4 $1/*cal?_[1-8]x3_125ns.bin -n -o $3/$4_anfit &
$binpath/meeg_t0res $3/$4.base $3/$4 $1/*cal?_[1-8]x3_125ns.bin -n -o $3/$4_linfit -u 

mkdir -p $3/neg
mv $3/*_neg.png $3/neg
