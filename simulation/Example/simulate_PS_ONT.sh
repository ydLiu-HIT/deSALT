#!/bin/bash
Ex_R=$1
Ref=$2
pbsim --sample-fastq $ExR --depth 4 $Ref
