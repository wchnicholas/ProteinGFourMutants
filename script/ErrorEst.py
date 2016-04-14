#!/usr/bin/python
import os
import sys
import gzip

def hashingpos(pos,nuc,poscount):
  if poscount.has_key(pos):
    if poscount[pos].has_key(nuc):
      poscount[pos][nuc] += 1
    else:
      poscount[pos][nuc] = 1
  else:
    poscount[pos] = {}
    poscount[pos][nuc] = 1
  return poscount

def ErrorEst(infile):
  infile    = open(infile,'r')
  poscount  = {}
  countline = 0
  for line in infile.xreadlines():
    countline += 1
    if countline%1000000==0: print 'Processed %d lines' % countline
    seq = line.rstrip()
    for n in range(len(seq)):
      nuc = seq[n]
      pos = n+1
      poscount = hashingpos(pos,nuc,poscount)
  infile.close()
  return poscount

def main():
  infile = 'paired/Input'
  poscount = ErrorEst(infile)
  print poscount

if __name__ == '__main__':
  main()
