#!/usr/bin/python
import sys
import operator
import random
import networkx as nx
from math import exp
from itertools import imap
from scipy.stats.stats import pearsonr

def hamming(str1, str2):
    assert len(str1) == len(str2)
    return sum(imap(operator.ne, str1, str2))

def TsvWithHeader2Hash(fitfile,condition):
  H = {}
  infile = open(fitfile,'r')
  countline = 0
  header = []
  for line in infile.xreadlines():
    countline += 1
    line = line.rstrip().rsplit("\t")
    if countline == 1: header = line; continue
    mut = line[0]
    if '_' in str(mut): continue
    for i in range(1,len(line)):
      if header[i] == condition and str(line[i]) != 'NA':
        H[mut] = {}
        H[mut][header[i]] = line[i] 
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
    fit  = exp(float(line[1]))
    #if mut in fithash.keys(): print 'Variant %d is not a missing data' % mut; sys.exit()
    fithash[mut] = {}
    fithash[mut][condition] = fit
  infile.close()
  return fithash

def genvarsneighbor(var):
  variants = []
  aas = ['E','D','R','K','H','Q','N','S','T','P','G','C','A','V','I','L','M','F','Y','W']
  [variants.append(aa+var[1]+var[2]+var[3]) for aa in aas]
  [variants.append(var[0]+aa+var[2]+var[3]) for aa in aas]
  [variants.append(var[0]+var[1]+aa+var[3]) for aa in aas]
  [variants.append(var[0]+var[1]+var[2]+aa) for aa in aas]
  while var in variants: variants.remove(var)
  return variants

def outcompile(muts,fithash,condition,outfile):
  l1 = []
  l2 = []
  countmut = 0
  outfile = open(outfile,'w')
  header  = "\t".join(['mut1','mut2','fit1','fit2'])
  outfile.write(header+"\n")
  for i in muts:
    countmut += 1
    if countmut%10000 == 0: print 'Finished processing %d variants' % countmut
    mutfit    = float(fithash[i][condition])
    variants  = genvarsneighbor(i)
    for var in variants:
      varfit = float(fithash[var][condition])
      l1.append(mutfit)
      l2.append(varfit)
      outfile.write("\t".join(map(str,[i,var,mutfit,varfit]))+"\n")
  outfile.close()
  print 'Pearson\'s Correlation = %f' % float(pearsonr(l1, l2)[0])
  
  

def main():
  fitfile      = 'result/Mutfit'
  #missfitfile  = 'result/regression_missing'
  missfitfile  = 'result/regression_all_WT'
  outfile      = 'analysis/LocalMaxNeighCor'+'_pair'
  fcutoff      = -1
  condition    = 'I20fit'
  Index2pos    = {0:39,1:40,2:41,3:54}
  fithash      = TsvWithHeader2Hash(fitfile,condition)
  print "Total # of variants in the raw data: %d" % len(fithash.keys())
  fithash      = fillinmissing(fithash,missfitfile,condition) 
  muts         = fithash.keys()
  print "Total # of variants after fill in with regression: %d" % len(muts)
  print "# of mutant pass cutoff: %d" % len(muts)
  print 'Start calculating pairwise fitness correlation'
  outcompile(muts,fithash,condition,outfile)


if __name__ == '__main__':
  main()
