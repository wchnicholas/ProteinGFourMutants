#!/usr/bin/python
#Find fitness effect of a particular mutant in different backgrounds
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

def hamming(str1, str2):
    assert len(str1) == len(str2)
    return sum(imap(operator.ne, str1, str2))

def callmut(var1, var2, Index2pos): #mutation for converting var1 to var2
  mut = []
  assert(len(var1)==len(var2))
  for i in range(len(var1)):
    if var1[i] != var2[i]:
      mut.append(var1[i]+str(Index2pos[i])+var2[i])
  return '-'.join(mut)

def callcommon(var1, var2, Index2pos): #mutation for converting var1 to var2
  bg = []
  assert(len(var1)==len(var2))
  for i in range(len(var1)):
    if var1[i] == var2[i]:
      bg.append(str(Index2pos[i])+var1[i])
  return '-'.join(bg)

def mut2index(mut, pos2index, var): #mutation for converting var1 to var2
  muts  = mut.rsplit('-')
  var   = list(var)
  for m in muts:
    pos = pos2index[int(m[1:-1])]
    assert(var[pos]==m[0])
    var[pos] = m[-1]
  return ''.join(var)

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

def genepihash(daalist,epiofint):
  epihash = {}
  epiofint_aa1 = epiofint.rsplit('-')[0]
  epiofint_aa2 = epiofint.rsplit('-')[1]
  epiofint_pos = set([epiofint_aa1[1:-1],epiofint_aa2[1:-1]])
  for daa_bg in daalist: 
    bg_aa1  = daa_bg.rsplit('-')[0]
    bg_aa2  = daa_bg.rsplit('-')[1]
    bg_pos  = set([bg_aa1[1:-1], bg_aa2[1:-1]])
    if len(list(epiofint_pos.intersection(bg_pos))) != 0: continue
    daa_bg  = bg_aa1[1::]+'-'+bg_aa2[1::]
    epihash[daa_bg] = {}
    for daa_epi in daalist:
      if daa_epi != epiofint: continue
      epi_aa1 = daa_epi.rsplit('-')[0]
      epi_aa1 = daa_epi.rsplit('-')[0]
      epi_aa2 = daa_epi.rsplit('-')[1]
      epi_pos = set([epi_aa1[1:-1], epi_aa2[1:-1]])
      if epi_aa1[0]==epi_aa1[-1] or epi_aa2[0]==epi_aa2[-1]: continue
      if len(list(epi_pos.intersection(bg_pos))) == 0: epihash[daa_bg][daa_epi] = 99
  return epihash

def genaalist(WT, Index2pos):
  daalist = []
  saalist = []
  aas = ['E','D','R','K','H','Q','N','S','T','P','G','C','A','V','I','L','M','F','Y','W']
  for i in range(len(WT)):
    for aa1 in aas:
      m1 = WT[i]+str(Index2pos[i])+aa1
      saalist.append(m1)
      for j in range(len(WT)):
        if i < j: 
          for aa2 in aas:
            m2 = WT[j]+str(Index2pos[j])+aa2
            daalist.append(m1+'-'+m2)
  return saalist, daalist

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

def analysis(fithash,epihash,WT,condition,Index2pos,pos2index,epiofint):
  print 'start analysis2. Condition = %s WT = %s TotalVariant = %d' % (condition, WT, len(fithash.keys()))
  countvar = 0
  for var1 in fithash.keys():
    countvar += 1
    if countvar%10000 == 0: print 'Processed %d variants' % countvar
    variants = genvars(var1,WT)
    varfit   = float(fithash[var1][condition])
    for var2 in variants:
      if hamming(var1, var2) != 2: continue
      mut = callmut(var1, var2, Index2pos)
      bg  = callcommon(var1, var2, Index2pos)
      if mut != epiofint: continue
      m1  = mut.rsplit('-')[0]
      m2  = mut.rsplit('-')[1]
      varm1 = mut2index(m1, pos2index, var1)
      varm2 = mut2index(m2, pos2index, var1)
      vard  = mut2index(mut, pos2index, var1)
      assert(vard == var2)
      if fithash.has_key(varm1) and fithash.has_key(varm2) and fithash.has_key(vard):
        varm1rawfit = float(fithash[varm1][condition])
        varm2rawfit = float(fithash[varm2][condition])
        vardrawfit  = float(fithash[vard][condition])
        varm1fit = floor(varm1rawfit)/floor(varfit)
        varm2fit = floor(varm2rawfit)/floor(varfit)
        vardfit  = floor(vardrawfit)/floor(varfit)
        varabsepi, varrelepi, cap = epistasiscal(varm1fit,varm2fit,vardfit,varm1rawfit,varm2rawfit,vardrawfit)
        epihash[bg][mut] = log(varrelepi)
  return epihash

def heatmapping2(epihash,outfile,epiofint,saalist):
  epiofint_aa1 = epiofint.rsplit('-')[0]
  epiofint_aa2 = epiofint.rsplit('-')[1]
  epiofint_pos1 = epiofint_aa1[1:-1]
  epiofint_pos2 = epiofint_aa2[1:-1]
  outfile = open(outfile,'w')
  header  = [saa[1::] for saa in saalist if epiofint_pos1 not in saa and epiofint_pos2 not in saa]
  header1 = header[0:len(header)/2]
  header2 = header[len(header)/2::]
  outfile.write('mut'+"\t"+"\t".join(header2)+"\n")
  for i in range(len(header1)):
    out = [header[i]]
    for j in range(len(header2)):
      bg1  = header1[i]
      bg2  = header2[j]
      bg   = bg1+'-'+bg2
      assert(len(epihash[bg].keys())==1)
      out.append(epihash[bg][epiofint])
    outfile.write("\t".join(map(str,out))+"\n")
  outfile.close()

def main():
  epiofint   = 'G41F-V54A'
  fitfile    = 'result/Mutfit'
  outfile    = 'analysis/Heatmap_'+epiofint.replace('-','')
  condition  = 'I20fit'
  WT         = 'VDGV'
  Index2pos  = {0:39,1:40,2:41,3:54}
  pos2index  = {39:0,40:1,41:2,54:3}
  fithash    = TsvWithHeader2Hash(fitfile)
  fithash    = filterfithash(fithash)
  saalist, daalist = genaalist(WT, Index2pos)
  epihash    = genepihash(daalist,epiofint)
  epihash    = analysis(fithash,epihash,WT,condition,Index2pos,pos2index,epiofint)
  heatmapping2(epihash,outfile,epiofint,saalist)

if __name__ == '__main__':
  main()

