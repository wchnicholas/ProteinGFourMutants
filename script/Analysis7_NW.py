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
              epihash[m1+'-'+m2] = {'WTepi':'NA',
                                    'WTm1fit':'NA',
                                    'WTm2fit':'NA',
                                    'WTDfit':'NA',
                                    'count':0,
                                    'PosCount':0,
                                    'NegCount':0,
                                    'maxepimut':'NA',
                                    'maxepi':-1,
                                    'maxm1fit':'NA',
                                    'maxm2fit':'NA',
                                    'maxDfit':'NA',
                                    'minepimut':'NA',
                                    'minepi':9999,
                                    'minm1fit':'NA',
                                    'minm2fit':'NA',
                                    'minDfit':'NA',
                                    'ALL':[]}
  return epihash

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

def analysis(fithash,epihash,WT,condition,Index2pos,pos2index):
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
      m1  = mut.rsplit('-')[0]
      m2  = mut.rsplit('-')[1]
      WTm1  = mut2index(m1, pos2index, WT)
      WTm2  = mut2index(m2, pos2index, WT)
      WTd   = mut2index(mut, pos2index, WT)
      varm1 = mut2index(m1, pos2index, var1)
      varm2 = mut2index(m2, pos2index, var1)
      vard  = mut2index(mut, pos2index, var1)
      assert(vard == var2)
      if fithash.has_key(WTm1) and fithash.has_key(WTm2) and fithash.has_key(WTd):
        WTm1fit = floor(float(fithash[WTm1][condition]))
        WTm2fit = floor(float(fithash[WTm2][condition]))
        WTdfit  = floor(float(fithash[WTd][condition]))
        WTabsepi, WTrelepi, cap = epistasiscal(WTm1fit,WTm2fit,WTdfit,WTm1fit,WTm2fit,WTdfit)
        epihash[mut]['WTDfit']  = WTdfit
        epihash[mut]['WTm1fit'] = WTm1fit
        epihash[mut]['WTm2fit'] = WTm2fit
        epihash[mut]['WTepi']   = WTrelepi
      if fithash.has_key(varm1) and fithash.has_key(varm2) and fithash.has_key(vard):
        varm1rawfit = float(fithash[varm1][condition])
        varm2rawfit = float(fithash[varm2][condition])
        vardrawfit  = float(fithash[vard][condition])
        varm1fit = floor(varm1rawfit)/floor(varfit)
        varm2fit = floor(varm2rawfit)/floor(varfit)
        vardfit  = floor(vardrawfit)/floor(varfit)
        varabsepi, varrelepi, cap = epistasiscal(varm1fit,varm2fit,vardfit,varm1rawfit,varm2rawfit,vardrawfit)
        epihash[mut]['count'] += 1
        signepi = ''
        #if varm1fit < 1 and varm2fit < 1 and vardfit > 1: signepi = 'pos'
        #if varm1fit > 1 and varm2fit > 1 and vardfit < 1: signepi = 'neg'
        if varrelepi > 1: signepi = 'pos'; epihash[mut]['PosCount'] += 1
        if varrelepi < 1: signepi = 'neg'; epihash[mut]['NegCount'] += 1
        if signepi == signepi:
          epihash[mut]['ALL'].append(varrelepi)
          if varrelepi > epihash[mut]['maxepi']:# and signepi == 'pos': 
            epihash[mut]['maxepimut'] = var1+'->'+var2
            epihash[mut]['maxepi']    = varrelepi
            epihash[mut]['maxm1fit']  = varm1fit
            epihash[mut]['maxm2fit']  = varm2fit
            epihash[mut]['maxDfit']   = vardfit
          if varrelepi < epihash[mut]['minepi']:#and signepi == 'neg': 
            epihash[mut]['minepimut'] = var1+'->'+var2
            epihash[mut]['minepi'] = varrelepi
            epihash[mut]['minm1fit']  = varm1fit
            epihash[mut]['minm2fit']  = varm2fit
            epihash[mut]['minDfit']   = vardfit
  return epihash

def output(epihash, outfile, condition):
  outfile = open(outfile+condition,'w')
  header = ['Mut','WTEpi','WTm1Fit','WTm2Fit','WTDFit','Count','PosCount','NegCount',
            'MaxEpiMut','MaxEpi','MaxM1Fit','MaxM2Fit','MaxDFit',
            'MinEpiMut','MinEpi','MinM1Fit','MinM2Fit','MinDFit','EpiSD','EpiRange']
  outfile.write("\t".join(header)+"\n")
  for mut in epihash.keys():
    WTepi       = epihash[mut]['WTepi']
    WTm1fit     = epihash[mut]['WTm1fit']
    WTm2fit     = epihash[mut]['WTm2fit']
    WTDfit      = epihash[mut]['WTDfit']
    count       = epihash[mut]['count']
    PosCount    = epihash[mut]['PosCount']
    NegCount    = epihash[mut]['NegCount']
    maxepimut   = epihash[mut]['maxepimut']
    maxepi      = epihash[mut]['maxepi']
    maxm1fit    = epihash[mut]['maxm1fit']
    maxm2fit    = epihash[mut]['maxm2fit']
    maxDfit     = epihash[mut]['maxDfit']
    minepimut   = epihash[mut]['minepimut']
    minepi      = epihash[mut]['minepi']
    minm1fit    = epihash[mut]['minm1fit']
    minm2fit    = epihash[mut]['minm2fit']
    minDfit     = epihash[mut]['minDfit']
    if len(epihash[mut]['ALL']) <= 1: episd = 0
    else: episd = np.std(map(log,epihash[mut]['ALL']))
    if maxepi > -1 and minepi < 9999: epirange = log(float(maxepi))-log(float(minepi))
    else: epirange = 'NA'
    out = "\t".join(map(str,[mut, WTepi,WTm1fit,WTm2fit,WTDfit,count,PosCount,NegCount,
                             maxepimut,maxepi,maxm1fit,maxm2fit,maxDfit,
                             minepimut,minepi,minm1fit,minm2fit,minDfit,episd,epirange]))
    outfile.write(out+"\n")
  outfile.close()

def main():
  fitfile    = 'result/Mutfit'
  outfile    = 'analysis/EpiDiffBG'
  condition  = 'I20fit'
  WT         = 'VDGV'
  Index2pos  = {0:39,1:40,2:41,3:54}
  pos2index  = {39:0,40:1,41:2,54:3}
  fithash    = TsvWithHeader2Hash(fitfile)
  fithash    = filterfithash(fithash)
  epihash    = genepihash(WT,Index2pos)
  epihash    = analysis(fithash,epihash,WT,condition,Index2pos,pos2index)
  output(epihash,outfile,condition)

if __name__ == '__main__':
  main()

