#!/usr/bim/python
import os
import sys
import random
import operator
import numpy as np
import networkx as nx
import matplotlib.pyplot as plt
from math import exp
from itertools import imap
from collections import Counter

def hamming(str1, str2):
    assert len(str1) == len(str2)
    return sum(imap(operator.ne, str1, str2))

def TsvWithHeader2Hash(fitfile):
  H = {}
  infile = open(fitfile,'r')
  countline = 0
  header = []
  for line in infile.xreadlines():
    countline += 1
    line = line.rstrip().rsplit("\t")
    if countline == 1: header = line; continue
    mut = line[0]
    #if mut[1] != 'D': continue ####################
    H[mut] = {}
    for i in range(1,len(line)): H[mut][header[i]] = line[i]
  infile.close()
  return H

def filterfithash(fithash,condition,fcutoff):
  for mut in fithash.keys():
    if '_' in mut or fithash[mut][condition] == 'NA' or float(fithash[mut][condition]) < float(fcutoff):
      del fithash[mut]
  return fithash

def fillinmissing(fithash,missfitfile,condition):
  infile = open(missfitfile,'r')
  for line in infile.xreadlines():
    if 'genotype_' in line: continue
    line = line.rstrip().rsplit("\t")
    mut  = line[0]
    #if mut[1] != 'D': continue ####################
    fit  = exp(float(line[1]))
    fithash[mut] = {}
    fithash[mut][condition] = fit
  infile.close()
  return fithash

def countingedges(samplenodes):
  countedge = 0 
  avgHD     = []
  for i in range(len(samplenodes)):
    for j in range(len(samplenodes)):
      if i < j:
        HD = hamming(samplenodes[i],samplenodes[j])
        avgHD.append(HD)
        if HD == 2: countedge += 1
  return countedge, np.mean(avgHD)


def main():
  WT          = 'VDGV'
  fitfile     = 'result/Mutfit'
  missfitfile = 'result/regression_missing'
  #missfitfile = 'result/regression_all_WT'
  condition   = 'I20fit'
  fcutoff     = -1
  fithash     = TsvWithHeader2Hash(fitfile)
  print "Total # of variants in the raw data: %d" % len(fithash.keys())
  fithash     = filterfithash(fithash,condition,fcutoff)
  print "Total # of variants pass filter of raw data: %d" % len(fithash.keys())
  fithash     = fillinmissing(fithash,missfitfile,condition) 
  muts        = fithash.keys()
  print "Total # of variants after fill in with regression: %d" % len(muts)
  print "# of mutant pass cutoff: %d" % len(muts)
  outfile     = 'analysis/RandomSampled15_All'
  #outfile     = 'analysis/RandomSampled15_Ben'
  #muts = [mut for mut in muts if float(fithash[mut][condition]) > 1]
  #print "# of mutant pass constraint: %d" % len(muts)
  outfile = open(outfile,'w')
  for i in range(0,1000000):
    if i%100000 == 0: print 'Completed %d bootstrap' % i
    samplenodes = random.sample(muts, 15)
    #samplenodes = ['IGQV','WNWY','WYGW','IYGC','LYGV','THCA','TYGM','FWGS','FYGN','ANLG','IWGF','VAAA','FWLG','ANCA','FWAA']
    countedge, avgHD = countingedges(samplenodes)
    outfile.write(str(avgHD)+"\n")
  outfile.close()
  

if __name__ == '__main__':
  main()
