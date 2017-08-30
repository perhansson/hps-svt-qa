#!/bin/bash
if [ -z "$1" ]
then
	echo "tp.sh <QA dir, e.g. /u1/data/l123_qa/hm11/7>"
	exit
fi
pushd "${0/%`basename $0`/}"
scriptpath="$PWD"
binpath="${PWD/%scripts/}/bin"
popd

runnum=$(basename $1)
device=$(basename $(dirname $1))
datadir="$1/data/"
plotsdir="$1/plots/"

mkdir -p $datadir $plotsdir

tpfiles=""
for calgroup in `seq 0 7`
	do for delay in `seq 1 8`
		do tpfiles="$tpfiles ${datadir}/${device}_${runnum}_cal_g${calgroup}_d${delay}.bin"
	done
done

$binpath/meeg_valid "${datadir}/${device}_${runnum}_baseline_dtrig.bin" $tpfiles -t1 > "${plotsdir}/${device}_${runnum}_summary.txt"

$binpath/meeg_baseline "${datadir}/${device}_${runnum}_baseline_dtrig.bin" -t1 -o "${plotsdir}/${device}_${runnum}"
$binpath/meeg_baseline "${datadir}/${device}_${runnum}_baseline_dtrig.bin" -t1 -m -o "${plotsdir}/${device}_${runnum}_mux"

$binpath/meeg_tp $tpfiles -t1 -fr -s20 -o "${plotsdir}/${device}_${runnum}"
