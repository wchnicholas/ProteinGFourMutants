#!/usr/bin/python
import os
import sys
import operator
import numpy as np
from math import log
from itertools import imap

def floor(fit):
  if fit == 'NA': return 'NA'
  elif float(fit) < 0.01: return 0.01
  else: return float(fit)

def mut2index(mut, pos2index, var): #mutation for converting var1 to var2
  muts  = mut.rsplit('-')
  var   = list(var)
  for m in muts:
    pos = pos2index[int(m[1:-1])]
    assert(var[pos]==m[0])
    var[pos] = m[-1]
  return ''.join(var)

def hamming(str1, str2):
    assert len(str1) == len(str2)
    return sum(imap(operator.ne, str1, str2))

def genvars(var,WT):
  variants = []
  aas = ['E','D','R','K','H','Q','N','S','T','P','G','C','A','V','I','L','M','F','Y','W']
  if var[0] == WT[0] and var[1] == WT[1]:
    for aa1 in aas: [variants.append(aa1+aa2+var[2]+var[3]) for aa2 in aas]
  if var[0] == WT[0] and var[2] == WT[2]:
    for aa1 in aas: [variants.append(aa1+var[1]+aa2+var[3]) for aa2 in aas]
  if var[0] == WT[0] and var[3] == WT[3]:
    for aa1 in aas: [variants.append(aa1+var[1]+var[2]+aa2) for aa2 in aas]
  if var[1] == WT[1] and var[2] == WT[2]:
    for aa1 in aas: [variants.append(var[0]+aa1+aa2+var[3]) for aa2 in aas]
  if var[1] == WT[1] and var[3] == WT[3]:
    for aa1 in aas: [variants.append(var[0]+aa1+var[2]+aa2) for aa2 in aas]
  if var[2] == WT[2] and var[3] == WT[3]:
    for aa1 in aas: [variants.append(var[0]+var[1]+aa1+aa2) for aa2 in aas]
  return variants

def callmut(var1, var2, Index2pos): #mutation for converting var1 to var2
  mut = []
  assert(len(var1)==len(var2))
  for i in range(len(var1)):
    if var1[i] != var2[i]:
      mut.append(var1[i]+str(Index2pos[i])+var2[i])
  return '-'.join(mut)

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

def genepihash(WT, Index2pos):
  aas = ['E','D','R','K','H','Q','N','S','T','P','G','C','A','V','I','L','M','F','Y','W']
  epihash = {}
  for i in range(len(WT)):
    for aa1 in aas:
      if aa1 != WT[i]:
        for j in range(len(WT)):
          for aa2 in aas:
            if aa2 != WT[j] and i < j:
              m1 = WT[i]+str(Index2pos[i])+aa1
              m2 = WT[j]+str(Index2pos[j])+aa2
              epihash[m1+'-'+m2] = {'RE':0,
                                    'ME':0,
                                    'SE':0}
  return epihash

def epistasiscal(var1fit,var2fit,dfit,varm1rawfit,varm2rawfit,vardrawfit):
  if var1fit < 0.01 or var2fit < 0.01 or varm1rawfit < 0.01 or varm2rawfit < 0.01: cap = 'NoNeg'
  if dfit < 0.01 or vardrawfit < 0.01: cap = 'NoPos'
  else: cap ='All'
  absepi = float(dfit) - floor(var1fit*var2fit) #Absolute Epistasis model
  relepi = float(floor(dfit))/floor(var1fit*var2fit) #Relative Epistasis model (default)
  if cap == 'NoNeg' and relepi < 1: relepi = 1
  elif cap == 'NoPos' and relepi > 1: relepi = 1
  elif varm1rawfit < 0.01 and varm2rawfit < 0.01 and vardrawfit < 0.01: relepi = 1
  else: relepi = relepi
  return absepi, relepi, cap

def analysis(fithash,epihash,WT,condition,Index2pos,pos2index,outfile):
  print 'start analysis2. Condition = %s WT = %s TotalVariant = %d' % (condition, WT, len(fithash.keys()))
  countvar = 0
  outfile  = open(outfile,'w')
  for var1 in fithash.keys():
    countvar += 1
    if countvar%10000 == 0: print 'Processed %d variants' % countvar
    variants = genvars(var1,WT)
    for var2 in variants:
      if hamming(var1, var2) != 2: continue
      mut = callmut(var1, var2, Index2pos)
      m1  = mut.rsplit('-')[0]
      m2  = mut.rsplit('-')[1]
      varm1 = mut2index(m1, pos2index, var1)
      varm2 = mut2index(m2, pos2index, var1)
      vard  = mut2index(mut, pos2index, var1)
      if fithash.has_key(varm1) and fithash.has_key(varm2) and fithash.has_key(vard):
        varfit   = float(fithash[var1][condition])
        varm1rawfit = float(fithash[varm1][condition])
        varm2rawfit = float(fithash[varm2][condition])
        vardrawfit  = float(fithash[vard][condition])
        allfits     = [varfit,varm1rawfit,varm2rawfit,vardrawfit]
        varm1fit = floor(varm1rawfit)/floor(varfit)
        varm2fit = floor(varm2rawfit)/floor(varfit)
        vardfit  = floor(vardrawfit)/floor(varfit)
        varabsepi, varrelepi, cap = epistasiscal(varm1fit,varm2fit,vardfit,varm1rawfit,varm2rawfit,vardrawfit) 
        if cap == 'All' and max(allfits) > 1 and min(allfits) < 0.2:
          outfile.write("\t".join(map(str, [mut, var1, var2, varfit, varm1rawfit, varm2rawfit, vardrawfit, varrelepi]))+"\n")
  outfile.close()

def main():
  fitfile    = 'result/Mutfit'
  outfile    = 'analysis/AllPairwiseEpi'
  condition  = 'I20fit'
  WT         = 'VDGV'
  Index2pos  = {0:39,1:40,2:41,3:54}
  pos2index  = {39:0,40:1,41:2,54:3}
  fithash    = TsvWithHeader2Hash(fitfile)
  fithash    = filterfithash(fithash)  
  epihash    = genepihash(WT,Index2pos)
  analysis(fithash,epihash,WT,condition,Index2pos,pos2index,outfile)

if __name__ == "__main__":
  main()
