#!/bin/sh

GCCVERSION=4.9.3
ROOTVERSION=6.05.02
ARCH=x86_64-slc6-gcc49-opt

source /afs/cern.ch/sw/lcg/external/gcc/${GCCVERSION}/${ARCH}/setup.sh
source /afs/cern.ch/sw/lcg/app/releases/ROOT/${ROOTVERSION}/${ARCH}/root/bin/thisroot.sh

export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:`pwd`/lib

