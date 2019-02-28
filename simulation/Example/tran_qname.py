#!/usr/bin/env python
# coding=utf-8

import sys

def trans(fqpath, prefix, wpath):
    fq = open(fqpath, 'r')
    fq2 = open(wpath, 'w')
    lines = fq.readlines()
    cnt = 4
    for line in lines:
        if (cnt%4 == 0) or (cnt%4 == 2):
            fq2.write(line.replace('S', prefix))
        else:
            fq2.write(line)
        cnt += 1

if __name__ == '__main__':
    fqpath = sys.argv[1] # previous fastq file
    prefix = sys.argv[2] # the prefix of read name by group
    wpath = sys.argv[3]  # fastq file after change name by group
    trans(fqpath, prefix, wpath)

