#!/usr/bin/python
import sys
import operator
import random
from itertools import imap

def hamming(str1, str2):
    assert len(str1) == len(str2)
    return sum(imap(operator.ne, str1, str2))

def main():
  WT           = 'VDGV'
  #infile       = 'tmp/pathofinterest'
  infile       = 'analysis/LocalMaxEvolvePotWT_pair'
  infile       = open(infile,'r')
  for line in infile.xreadlines():
    if 'mut' in line: continue
    line = line.rstrip().rsplit("\t")
    mut  = line[0]
    HD   = int(line[1])
    fit  = line[2]
    step = int(line[3])
    path = line[4].rsplit('->')
    if step-HD != 1: continue
    if HD==4 and hamming(path[1],WT) == 1 and hamming(path[1],mut) == HD: print line
  infile.close()

if __name__ == '__main__':
  main()
