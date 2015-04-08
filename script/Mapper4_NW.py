#!/usr/bin/python
import os
import sys
import operator
from string import atof
from itertools import imap

def hamming(str1, str2):
    assert len(str1) == len(str2)
    return sum(imap(operator.ne, str1, str2))


def hashin(sample, countfile, counthash):
  print 'Reading sample: %s' % sample
  infile = open(countfile,'r')
  counthash[sample] = {}
  for line in infile.xreadlines():
    line = line.rstrip().rsplit(' ')
    while '' in line: line.remove('')
    count = line[0]
    mutID = line[1]
    counthash[sample][mutID] = count
  infile.close()
  return counthash

def output(WToutfile, Mutoutfile, counthash, samples, WT):
  WToutfile  = open(WToutfile,'w')
  Mutoutfile = open(Mutoutfile,'w')
  muts       = []
  [muts.extend(counthash[sample].keys()) for sample in samples]
  muts       = list(set(muts))
  header     = 'mut'+"\t"+'HD'+"\t"+"\t".join(samples)
  WToutfile.write(header+"\n")
  Mutoutfile.write(header+"\n")
  for mut in muts:
    out = [mut,str(hamming(WT,mut))]
    for sample in samples: 
      if counthash[sample].has_key(mut): out.append(counthash[sample][mut])
      else: out.append('0')
    out = "\t".join(out)
    Mutoutfile.write(out+"\n")
    if mut == WT: WToutfile.write(out+"\n")
  WToutfile.close()
  Mutoutfile.close()

def main():
  samples    = ['Input','IGG10','IGG20','IGG90']
  WT         = 'VDGV'
  WToutfile  = 'result/WTcount'
  Mutoutfile = 'result/Mutcount'
  counthash  = {}
  for sample in samples: 
    countfile = 'count/'+sample+'.count'
    counthash = hashin(sample, countfile, counthash)
  output(WToutfile, Mutoutfile, counthash, samples, WT)

if __name__ == '__main__':
  main()
