#!/bin/bash
ex_R=$1
Ref=$2
python read_analysis.py genome -i $ex_R -rg $Ref

python simulator.py genome -rg $Ref -n 1000
