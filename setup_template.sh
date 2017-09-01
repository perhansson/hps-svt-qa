
export QA_DIR=$PWD

# Base DAQ directory must exist
if [ ! -d "$PWD/daq" ]; then
    echo "Need to create symlink to a daq/ directory; e.g. \$ln -s svt-daq daq."
	return 1
fi
export BASE=$QA_DIR/daq

if [ "$MACHINE" == "" ]; then
    export MACHINE=`uname -m`
fi

if [ "$ROOTSYS" == "" ]; then
	source /afs/slac.stanford.edu/g/atlas/packages/root/root_v5.34.05/bin/thisroot.sh
fi



# Setup library path
#if ($?LD_LIBRARY_PATH) then
#   setenv LD_LIBRARY_PATH ${ROOTSYS}/lib:${QTDIR}/lib:${LD_LIBRARY_PATH}
#else
#   setenv LD_LIBRARY_PATH ${ROOTSYS}/lib:${QTDIR}/lib
#endif

