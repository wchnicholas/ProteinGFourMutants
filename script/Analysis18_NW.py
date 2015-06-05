#!/usr/bin/python
import sys
import operator
from math import exp
from itertools import imap

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

def hashinbulkfile(infile):
  bhash = {}
  infile = open(infile,'r')
  for line in infile.xreadlines():
    if 'aa' in line: continue
    line = line.rstrip().rsplit("\t")
    aa   = line[0]
    mass = line[2]
    bhash[aa] = float(mass)
  infile.close()
  return bhash

def outcompile(muts,fithash,condition,masshash,outfile):
  outfile = open(outfile,'w')
  header  = "\t".join(['mut','HD','fit','coremass','mass39','mass40','mass41','mass54'])
  outfile.write(header+"\n")
  coreaa  = ['G','A','I','L','F','V']
  coreaa  = ['A','F']
  for mut in muts:
    resi39 = mut[0]
    resi40 = mut[1]
    resi41 = mut[2]
    resi54 = mut[3]
    mass39 = masshash[resi39]
    mass40 = masshash[resi40]
    mass41 = masshash[resi41]
    mass54 = masshash[resi54]
    #if resi39 in coreaa and resi41 in coreaa and resi54 in coreaa and resi40 == 'D':
    if True:
      coremass = mass39+mass41+mass54
      mutfit   = fithash[mut][condition]
      HD       = hamming(mut,'VDGV')
      outfile.write("\t".join(map(str,[mut,HD,mutfit,coremass,mass39,mass40,mass41,mass54]))+"\n")
  outfile.close()
      
def main():
  fitfile      = 'result/Mutfit'
  outfile      = 'analysis/BulkvsFit'
  massfile     = 'doc/AAmass'
  fcutoff      = -1
  condition    = 'I20fit'
  Index2pos    = {0:39,1:40,2:41,3:54}
  masshash     = hashinbulkfile(massfile)
  fithash      = TsvWithHeader2Hash(fitfile,condition)
  print "Total # of variants in the raw data: %d" % len(fithash.keys())
  muts         = fithash.keys()
  print "# of mutant pass cutoff: %d" % len(muts)
  outcompile(muts,fithash,condition,masshash,outfile)
  

if __name__ == '__main__':
  main()
