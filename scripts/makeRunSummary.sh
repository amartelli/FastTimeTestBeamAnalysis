#!/bin/bash

root -l -b -q makeRunSummary.C+\(0,2\)
root -l -b -q makeRunSummary.C\(0,3\) &
root -l -b -q makeRunSummary.C\(0,4\) &
root -l -b -q makeRunSummary.C\(1,2\)&
root -l -b -q makeRunSummary.C\(1,3\)&
root -l -b -q makeRunSummary.C\(1,4\) &
root -l -b -q makeRunSummary.C\(2,2\) &
root -l -b -q makeRunSummary.C\(2,3\) &
root -l -b -q makeRunSummary.C\(2,4\) &
