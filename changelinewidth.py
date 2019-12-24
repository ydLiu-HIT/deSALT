#!/usr/bin/env python
# coding=utf-8
import sys

def change_line(path1, path2, linewidth):
    f = open(path1, 'r')
    w = open(path2, 'w')

    lines = f.readlines()
    count = 0

    for line in lines:
        if '>' in line:
            w.write(line)
        else:
            length = len(line)
            while (count + linewidth) < length:
                w.write(line[count:(count + linewidth)])
                w.write('\n')
                count = count + linewidth
            w.write(line[count:length])
            count = 0

if __name__ == '__main__':
    input_fa = sys.argv[1]
    output_fa = sys.argv[2]
    linewidth = int(sys.argv[3])
    change_line(input_fa, output_fa, linewidth)
