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

sed "s/DEVICE/$device/g;s/RUN/$runnum/g" $scriptpath/plots_template.html > $1/plots.html

python $scriptpath/run_cal.py -t1 -c3 "${datadir}/${device}_${runnum}"
#$binpath/meeg_baseline "${datadir}/${device}_${runnum}_baseline_dtrig.bin" -t1 -o "${plotsdir}/${device}_${runnum}"

#$binpath/meeg_baseline "${datadir}/${device}_${runnum}_baseline_dtrig.bin" -t1 -m -o "${plotsdir}/${device}_${runnum}_mux"

#$binpath/meeg_tp "${datadir}/${device}_${runnum}_cal_g?_d?.bin" -t1 -fr -s20 -o "${plotsdir}/${device}_${runnum}"
