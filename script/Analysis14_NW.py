#!/usr/bin/python
import os
import sys
import glob
import operator
import numpy as np
from itertools import imap

def fithashing(infile):
  fithash = {}
  infile = open(infile,'r')
  for line in infile.xreadlines(): 
    if 'mut' in line: continue
    line = line.rstrip().rsplit("\t")
    mut  = line[0]
    fit  = line[1]
    fithash[mut] = fit
  infile.close()
  return fithash

def lenghashing(lengfiles,lenghash):
  countline = 0
  for infile in lengfiles:
    infile = open(infile,'r')
    for line in infile.xreadlines(): 
      if 'mut' in line: continue
      countline += 1
      #if countline == 1600000: break #FOR TESTING PURPOSE
      if countline%4000000 == 0: print 'Processed %d paths' % countline
      line = line.rstrip().rsplit("\t")
      mut  = line[0]
      leng = line[1]
      if lenghash.has_key(mut): lenghash[mut].append(float(leng))
      else: lenghash[mut] = [float(leng)]
    infile.close()
  return lenghash

def outcompile(lenghash,muts,fithash,outfile):
  outfile = open(outfile,'w')
  header = "\t".join(['mut','fit','avgdist2max','maxdist2max','mindist2max','sddist2max'])
  outfile.write(header+"\n")
  for mut in muts:
    paths = lenghash[mut]
    outfile.write("\t".join(map(str,[mut,fithash[mut],np.mean(paths),max(paths),min(paths),np.std(paths)]))+"\n")
  outfile.close()


def main():
  simtype      = 'random' #weight or random
  fitfile      = 'analysis/LocalMaxDes_'+simtype
  outfile      = 'analysis/LocalMaxDist_'+simtype
  lengfiles    = glob.glob('simulations/'+simtype+'/LocalMaxClimb_'+simtype+'*')
  lenghash     = {}
  lenghash     = lenghashing(lengfiles,lenghash)
  fithash      = fithashing(fitfile)
  muts         = lenghash.keys()
  assert(len(muts)==160000)
  outcompile(lenghash,muts,fithash,outfile)

if __name__ == '__main__':
  main()
