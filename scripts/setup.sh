#!/bin/sh

GCCVERSION=4.9.3
ARCH=x86_64-slc6-gcc49-opt
ROOTVERSION=6.05.02

source /afs/cern.ch/sw/lcg/external/gcc/${GCCVERSION}/${ARCH}/setup.sh
source /afs/cern.ch/sw/lcg/app/releases/ROOT/${ROOTVERSION}/${ARCH}/root/bin/thisroot.sh

export PATH=$PATH:/afs/cern.ch/project/eos/installation/cms/bin/
export LD_LIBRARY_PATH=$ROOTSYS/lib:$LD_LIBRARY_PATH:`pwd`/lib
