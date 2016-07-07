#!/usr/bin/python
import os
import sys
import glob
import numpy as np
from math import log, sqrt

def adjpos(muts):
  muts = muts.rsplit('-')
  newmuts = []
  for mut in muts:
    newmuts.append(mut[0]+str(int(mut[1:-1])+1)+mut[-1])
  return '-'.join(newmuts)

def convertmut(muts):
  muts    = muts.rsplit('-')
  mutid   = ['V','D','G','V']
  aahash  = {'39':0,'40':1,'41':2,'54':3}
  for mut in muts: 
    pos = mut[1:-1]
    aa  = mut[-1]
    mutid[aahash[pos]] = aa
  return ''.join(mutid)

def hashin(filename, H):
  infile = open(filename,'r')
  for line in infile.xreadlines():
    if 'Mut' in line: continue
    line = line.rstrip().rsplit("\t")
    mut = adjpos(line[0])
    DNAcount = float(line[1])
    DNAfreq  = DNAcount/float(line[7])
    Selfreq  = float(line[6])/float(line[12])
    fit = Selfreq/DNAfreq
    if DNAcount >= 10: #TUNABLE
      goodmut = 'yes'
      muts = mut.rsplit('-')
      for m in muts:
        if '39' not in m and '40' not in m and '41' not in m and '54' not in m: goodmut ='no'; break
      if goodmut == 'yes':
        adjmut = convertmut(mut)
        H[adjmut] = fit
  return H

def floorfit2Kd(fit,C):
  if fit < 0.01: return 0.01
  elif fit > C: return C*0.95
  else: return fit

def fit2relKd(fit,C):
  return 1.9858775*296/1000*log((C-1)/((C/floorfit2Kd(fit,C))-1))


def ComputeKd(fits,Bmaxs):
  Kd = []
  for i in range(len(fits)):
    Kd.append(fit2relKd(fits[i],Bmaxs[i]))
  return np.mean(Kd), np.std(Kd)/sqrt(len(Kd))

def fitting(Mutfile, WTIpt, WTI10, WTI20, WTI90, outfile, AOhash, Bmax10, Bmax20, Bmax90, BmaxAO):
  infile = open(Mutfile,'r')
  outfile = open(outfile,'w')
  countline = 0
  for line in infile.xreadlines():
    countline += 1
    line  = line.rstrip().rsplit("\t")
    if countline == 1: 
      outfile.write("\t".join(line)+"\t"+"\t".join(['I10fit', 'I20fitRaw', 'I90fit','AOfit','I20fit'])+"\n")
      continue
    Mut   = line[0]
    HD    = line[1]
    mIpt  = float(line[2])
    mI10  = float(line[3])
    mI20  = float(line[4])
    mI90  = float(line[5])
    fits  = []
    if mIpt < 10:
      I10_fit = 'NA'
      I20_fit = 'NA'
      I90_fit = 'NA'
      Kd_avg  = 'NA'
      Kd_std  = 'NA'
      fits    = [I10_fit, I20_fit, I90_fit]
    else:
      Iptfreq = mIpt/WTIpt
      I10_fit = (mI10/WTI10)/Iptfreq
      I20_fit = (mI20/WTI20)/Iptfreq
      I90_fit = (mI90/WTI90)/Iptfreq
      fits    = [I10_fit, I20_fit, I90_fit]
      #Kd_avg, Kd_std = ComputeKd(fits,[Bmax10,Bmax20,Bmax90])
    if AOhash.has_key(Mut):
      AOfit = AOhash[Mut]
      fits.append(AOfit)
      #Kd_AO_avg, Kd_AO_std = ComputeKd([AOhash[Mut]],[BmaxAO])
    else: 
      fits.append('NA')
      #Kd_AO_avg  = 'NA'
      #Kd_AO_std  = 'NA'
    if int(HD) == 1: I20_fitadj = float(AOfit)*1.15868
    else: I20_fitadj = I20_fit
    fits+=[I20_fitadj]
    out = "\t".join(map(str,line))+"\t"+"\t".join(map(str,fits))
    outfile.write(out+"\n")
  infile.close()
  outfile.close()

def main():
  WTinfo  = open('result/WTcount','r').readlines()[1].rstrip().rsplit("\t")
  WTIpt   = float(WTinfo[2])
  WTI10   = float(WTinfo[3])
  WTI20   = float(WTinfo[4])
  WTI90   = float(WTinfo[5])
  Bmax10  = 18.4
  Bmax20  = 7
  Bmax90  = 1.86
  BmaxAO  = 8
  Mutfile = 'result/Mutcount'
  outfile = 'result/Mutfit'
  AOhash  = {}
  print 'Starting Hashing Anders Data'
  AOhash  = hashin('doc/SMutList', AOhash)
  AOhash  = hashin('doc/DMutList', AOhash)
  print 'Starting Compiling data'
  fitting(Mutfile, WTIpt, WTI10, WTI20, WTI90, outfile, AOhash, Bmax10, Bmax20, Bmax90, BmaxAO)

if __name__ == '__main__':
  main()
