#!/bin/bash
find analy -print -delete
find histo -print -delete
find output -print -delete
find rdata_* -print -delete
mkdir analy
mkdir histo
mkdir output
rm -fv core.* fort.* *_out.out *_lsf.out *_out.err
rm -fv budget_* thermo_state_* photo_state_*
