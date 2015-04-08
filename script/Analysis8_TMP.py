#!/usr/bin/python
import os
import sys
import glob
import networkx as nx
import operator
from itertools import imap

def floor(fit):
  if fit == 'NA': return 'NA'
  elif float(fit) < 0.01: return 0.01
  else: return fit

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
    H[mut] = {}
    for i in range(1,len(line)): H[mut][header[i]] = line[i]
  infile.close()
  return H

def filterfithash(fithash):
  for mut in fithash.keys():
    if fithash[mut]['I10fit'] == 'NA' or '_' in mut:
      del fithash[mut]
  return fithash

def expectingfit(fithash, condition, muts):
  fits = []
  for mut in muts:
    if fithash.has_key(mut): fits.append(float(fithash[mut][condition]))
    else: return 'NA'
  if min(fits) < 0.01: return float(0.01)
  else:
    expfit = float(1)
    for fit in fits: expfit = fit*expfit
    return expfit

def callepistaticterm(mut,WT,fithash,condition,outE2,outE3,outE4,outE3_AO,outE4_AO):
  assert(hamming(mut,WT)==4)
  S1   = mut[0]+  WT[1]+  WT[2]+  WT[3]
  S2   =  WT[0]+ mut[1]+  WT[2]+  WT[3]
  S3   =  WT[0]+  WT[1]+ mut[2]+  WT[3]
  S4   =  WT[0]+  WT[1]+  WT[2]+ mut[3]
  D12  = mut[0]+ mut[1]+  WT[2]+  WT[3]
  D13  = mut[0]+  WT[1]+ mut[2]+  WT[3]
  D14  = mut[0]+  WT[1]+  WT[2]+ mut[3]
  D23  =  WT[0]+ mut[1]+ mut[2]+  WT[3]
  D24  =  WT[0]+ mut[1]+  WT[2]+ mut[3]
  D34  =  WT[0]+  WT[1]+ mut[2]+ mut[3]
  T123 = mut[0]+ mut[1]+ mut[2]+  WT[3]
  T124 = mut[0]+ mut[1]+  WT[2]+ mut[3]
  T134 = mut[0]+  WT[1]+ mut[2]+ mut[3]
  T234 =  WT[0]+ mut[1]+ mut[2]+ mut[3]
  Q1234= mut
  for i in [S1,S2,S3,S4,D12,D13,D14,D23,D34,T123,T124,T134,T234,mut]: 
    if not fithash.has_key(i): return 'NA'
  S1fit    = float(fithash[S1][condition])
  S2fit    = float(fithash[S2][condition])
  S3fit    = float(fithash[S3][condition])
  S4fit    = float(fithash[S4][condition])
  D12fit   = float(fithash[D12][condition])
  D13fit   = float(fithash[D13][condition])
  D14fit   = float(fithash[D14][condition])
  D23fit   = float(fithash[D23][condition])
  D24fit   = float(fithash[D24][condition])
  D34fit   = float(fithash[D34][condition])
  T123fit  = float(fithash[T123][condition])
  T124fit  = float(fithash[T124][condition])
  T134fit  = float(fithash[T134][condition])
  T234fit  = float(fithash[T234][condition])
  Q1234fit = float(fithash[Q1234][condition])

  S1expfit    = expectingfit(fithash, condition, [S1])
  S2expfit    = expectingfit(fithash, condition, [S2])
  S3expfit    = expectingfit(fithash, condition, [S3])
  S4expfit    = expectingfit(fithash, condition, [S4])
  D12expfit   = expectingfit(fithash, condition, [S1,S2])
  D13expfit   = expectingfit(fithash, condition, [S1,S3])
  D14expfit   = expectingfit(fithash, condition, [S1,S4])
  D23expfit   = expectingfit(fithash, condition, [S2,S3])
  D24expfit   = expectingfit(fithash, condition, [S2,S4])
  D34expfit   = expectingfit(fithash, condition, [S3,S4])
  T123expfit  = expectingfit(fithash, condition, [S1,S2,S3])
  T124expfit  = expectingfit(fithash, condition, [S1,S2,S4])
  T134expfit  = expectingfit(fithash, condition, [S1,S3,S4])
  T234expfit  = expectingfit(fithash, condition, [S2,S3,S4])
  Q1234expfit = expectingfit(fithash, condition, [S1,S2,S3,S4])

  D12epi   = D12fit/D12expfit
  D13epi   = D13fit/D13expfit
  D14epi   = D14fit/D14expfit
  D23epi   = D23fit/D23expfit
  D24epi   = D24fit/D24expfit
  D34epi   = D34fit/D34expfit
  T123epi  = T123fit/T123expfit
  T124epi  = T124fit/T124expfit
  T134epi  = T134fit/T134expfit
  T234epi  = T234fit/T234expfit
  Q1234epi = Q1234fit/Q1234expfit

  T123epi_AO  = T123fit/T123expfit*(D12epi*D13epi*D23epi)
  T124epi_AO  = T124fit/T124expfit*(D12epi*D14epi*D24epi)
  T134epi_AO  = T134fit/T134expfit*(D12epi*D13epi*D24epi)
  T234epi_AO  = T234fit/T234expfit*(D23epi*D24epi*D34epi)
  Q1234epi_AO = Q1234fit/Q1234expfit*(D12epi*D13epi*D14epi*D23epi*D24epi*D34epi*T123epi_AO*T124epi_AO*T134epi_AO*T234epi_AO)

  outE2.write(mut+"\t"+"\t".join(map(str,[D23epi,D24epi,D34epi]))+"\n")
  outE3.write(mut+"\t"+"\t".join(map(str,[T123epi,T124epi,T134epi,T234epi]))+"\n")
  outE4.write(mut+"\t"+"\t"+str(Q1234epi)+"\n")
  outE3_AO.write(mut+"\t"+"\t".join(map(str,[T123epi_AO,T124epi_AO,T134epi_AO,T234epi_AO]))+"\n")
  outE4_AO.write(mut+"\t"+"\t".join(map(str,[Q1234epi_AO]))+"\n")

def main():
  fitfile    = 'result/Mutfit'
  conditions = ['I10fit','I20fit','I90fit']
  condition  = conditions[1]
  WT         = 'VDGV'
  outE2      = 'result/Epi2'
  outE3      = 'result/Epi3'
  outE4      = 'result/Epi4'
  outE3_AO   = 'result/Epi3AO'
  outE4_AO   = 'result/Epi4AO'
  poshash    = {0:39,1:40,2:41,3:54}
  fithash    = TsvWithHeader2Hash(fitfile)
  fithash    = filterfithash(fithash)
  muts       = fithash.keys()
  outE2 = open(outE2+'_'+condition,'w')
  outE3 = open(outE3+'_'+condition,'w')
  outE4 = open(outE4+'_'+condition,'w')
  outE3_AO = open(outE3_AO+'_'+condition,'w')
  outE4_AO = open(outE4_AO+'_'+condition,'w')
  for mut in muts:
    HD = int(fithash[mut]['HD'])
    if HD == 4:
      callepistaticterm(mut,WT,fithash,condition,outE2,outE3,outE4,outE3_AO,outE4_AO)
  outE2.close()
  outE3.close()
  outE4.close()
  outE3_AO.close()
  outE4_AO.close()

if __name__ == '__main__':
  main()
