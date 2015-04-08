#!/usr/bin/python
#Find fitness effect of a particular mutant in different backgrounds
import os
import sys
import operator
from itertools import imap

def floor(fit):
  if fit == 'NA': return 'NA'
  elif float(fit) < 0.01: return 0.01
  else: return fit

def hamming(str1, str2):
    assert len(str1) == len(str2)
    return sum(imap(operator.ne, str1, str2))

def callmut(var1, var2, poshash): #mutation for converting var1 to var2
  mut = []
  assert(len(var1)==len(var2))
  for i in range(len(var1)):
    if var1[i] != var2[i]:
      mut.append(var1[i]+str(poshash[i])+var2[i])
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

def genmuthash(WT, poshash):
  aas = ['E','D','R','K','H','Q','N','S','T','P','G','C','A','V','I','L','M','F','Y','W']
  muthash = {}
  for i in range(len(WT)):
    for aa in aas:
      if aa != WT[i]:
        m = WT[i]+str(poshash[i])+aa
        muthash[m] = {'count':0,
                      'usefulcount':0,
                      'maxfit':-1,
                      'maxfitmut':'NA',
                      'minfit':9999,
                      'minfitmut':'NA',
                      'CrypticBenFit':-1,
                      'CrypticVar2Fit':-1,
                      'CrypticBenMut':'NA'}
  return muthash

def genvars(var,WT):
  variants = []
  aas = ['E','D','R','K','H','Q','N','S','T','P','G','C','A','V','I','L','M','F','Y','W']
  if var[0] == WT[0]: [variants.append(aa+var[1]+var[2]+var[3]) for aa in aas]
  if var[1] == WT[1]: [variants.append(var[0]+aa+var[2]+var[3]) for aa in aas]
  if var[2] == WT[2]: [variants.append(var[0]+var[1]+aa+var[3]) for aa in aas]
  if var[3] == WT[3]: [variants.append(var[0]+var[1]+var[2]+aa) for aa in aas]
  while var in variants: variants.remove(var)
  return variants

def analysis(fithash,muthash,WT,condition,poshash,HD):
  print 'start analysis1. Condition = %s WT = %s TotalVariant = %d. HD = %d' % (condition, WT, len(fithash.keys()), HD)
  countvar = 0
  fits = {}
  for var1 in fithash.keys():
    countvar += 1
    if countvar%10000 == 0: print 'Processed %d variants' % countvar
    variants = genvars(var1,WT)
    if hamming(var1,WT) != HD-1: continue
    for var2 in variants:
      if hamming(var1, var2) == 1 and fithash.has_key(var2) and hamming(var2,WT) == HD:
        mut = callmut(var1, var2, poshash)
        muthash[mut]['count'] += 1
        var1fit = float(fithash[var1][condition])
        var2fit = float(fithash[var2][condition])
        if var1fit < 0.01 and var2fit < 0.01: continue
        if var2fit < 0.01: mutfit = 0.01
        else: mutfit  = floor(var2fit)/floor(var1fit)
        if mutfit < muthash[mut]['minfit']: muthash[mut]['minfit'] = mutfit; muthash[mut]['minfitmut'] = var1+'->'+var2
        if mutfit > muthash[mut]['maxfit']: muthash[mut]['maxfit'] = mutfit; muthash[mut]['maxfitmut'] = var1+'->'+var2
        if mutfit > muthash[mut]['CrypticBenFit'] and var2fit > 1: 
          muthash[mut]['CrypticBenFit']  = mutfit
          muthash[mut]['CrypticVar2Fit'] = var2fit
          muthash[mut]['CrypticBenMut']  = var1+'->'+var2
        fits[var1+'->'+var2] = mutfit
  return muthash, fits

def output(muthash, outfile,condition, HD):
  for mut in muthash.keys():
    count  = muthash[mut]['count']
    maxfit = floor(muthash[mut]['maxfit'])
    maxmut = muthash[mut]['maxfitmut']
    minfit = floor(muthash[mut]['minfit'])
    minmut = muthash[mut]['minfitmut']
    CrypticBenFit  = muthash[mut]['CrypticBenFit']
    CrypticBenMut  = muthash[mut]['CrypticBenMut']
    CrypticVar2Fit = muthash[mut]['CrypticVar2Fit']
    if count == 0:
      maxfit = 'NA'
      minfit = 'NA'
    out = "\t".join(map(str,[mut,HD,count,maxfit,maxmut,minfit,minmut,CrypticBenFit,CrypticBenMut,CrypticVar2Fit]))
    outfile.write(out+"\n")

#MAIN
def main():
  fitfile    = 'result/Mutfit' 
  outfile    = 'analysis/MutDiffBG'
  conditions = ['I10fit','I20fit','I90fit']
  conditions = ['I20fit']
  HDs        = [1,2,3,4]
  WT         = 'VDGV'
  poshash    = {0:39,1:40,2:41,3:54}
  fithash    = TsvWithHeader2Hash(fitfile)
  fithash    = filterfithash(fithash)
  fitout     = open('analysis/FitDiffBG','w')
  print >>fitout, "condition\tmut\tHD\tfit"
  for condition in conditions:
    outfile = open(outfile+condition,'w')
    header = ['Mut','HD','CountBG','MaxFit','MaxFitVar','MinFit','MinFitVar','CrypticBenFit','CrypticBenMut','CrypticVar2Fit']
    outfile.write("\t".join(header)+"\n")
    for HD in HDs:
      muthash       = genmuthash(WT, poshash)
      muthash, fits = analysis(fithash,muthash,WT,condition,poshash, HD)
      output(muthash,outfile,condition,HD)
      for mut in fits.keys(): 
        fit = fits[mut]
        print >>fitout, condition+"\t"+str(mut)+"\t"+str(HD)+"\t"+str(floor(fit))
    outfile.close()
  fitout.close()
  
if __name__ == '__main__':
  main()
