#!/usr/bin/python
import os
import sys
import glob

def floor(fit):
  if fit == 'NA': return 'NA'
  elif float(fit) < 0.01: return 0.01
  else: return fit

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
  cap  = ''
  for mut in muts:
    if fithash.has_key(mut): fits.append(float(fithash[mut][condition]))
    else: return 'NA', 'NA', ['NA','NA','NA','NA']
  if min(fits) < 0.01: cap = 'NoNeg'
  else: cap = 'All'
  expfit = float(1)
  for fit in fits: expfit = fit*expfit
  return expfit, cap, fits

def epical(fithash,muts,WT,condition,outfile):
  for mut in muts:
    HD = int(fithash[mut]['HD'])
    mutfit = float(fithash[mut][condition])
    S1     = mut[0]+  WT[1]+  WT[2]+  WT[3]
    S2     =  WT[0]+ mut[1]+  WT[2]+  WT[3]
    S3     =  WT[0]+  WT[1]+ mut[2]+  WT[3]
    S4     =  WT[0]+  WT[1]+  WT[2]+ mut[3]
    expfit, cap, fits = expectingfit(fithash, condition, [S1, S2, S3, S4])
    if cap in ['NoNeg','All']: 
      epis   = floor(mutfit)/floor(expfit)
      if cap == 'NoNeg': epis = sorted([epis,1])[1]
    else: assert(cap=='NA')
    if 'NA' not in cap:
      outfile.write("\t".join(map(str,[condition, mut,HD,mutfit,expfit,fits[0],fits[1],fits[2],fits[3],cap,epis]))+"\n")
   
def main():
  WT         = 'VDGV'
  fitfile    = 'result/Mutfit'
  outfile    = 'result/AllEpi'
  fithash    = TsvWithHeader2Hash(fitfile)
  fithash    = filterfithash(fithash)
  muts       = fithash.keys()
  print 'Generating file: %s' % outfile
  outfile = open(outfile,'w')
  header = ['condition','mut','HD','mutfit','expfit','S1fit','S2fit','S3fit','S4fit','cap','epis']
  outfile.write("\t".join(header)+"\n")
  epical(fithash,muts,WT,'I10fit',outfile)
  epical(fithash,muts,WT,'I20fit',outfile)
  epical(fithash,muts,WT,'I90fit',outfile)
  outfile.close()

if __name__ == '__main__':
  main()
