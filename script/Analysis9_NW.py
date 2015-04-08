#!/usr/bin/python
import os
import sys
import glob
import networkx as nx
import operator
from itertools import imap, combinations

def floor(fit):
  if fit == 'NA': return 'NA'
  elif float(fit) < 0.01: return 0.01
  else: return fit

def hamming(str1, str2):
    assert len(str1) == len(str2)
    return sum(imap(operator.ne, str1, str2))

def mut2index(mut, Pos2Index, var): #mutation for converting var1 to var2
  muts  = mut.rsplit('-')
  var   = list(var)
  for m in muts:
    pos = Pos2Index[int(m[1:-1])]
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

def GenMutfromPos(MutPos,Index2Pos,WT):
  aas   = ['E','D','R','K','H','Q','N','S','T','P','G','C','A','V','I','L','M','F','Y','W']
  MutID = []
  MutPos1 = MutPos[0]
  MutPos2 = MutPos[1]
  MutPos3 = MutPos[2]
  for aa1 in aas:
    for aa2 in aas:
      for aa3 in aas:
        if aa1 != WT[MutPos1] and aa2 != WT[MutPos2] and aa3 != WT[MutPos3]: 
          M1 = WT[MutPos1]+str(Index2Pos[MutPos1])+aa1
          M2 = WT[MutPos2]+str(Index2Pos[MutPos2])+aa2
          M3 = WT[MutPos3]+str(Index2Pos[MutPos3])+aa3
          MutID.append('-'.join([M1,M2,M3])) 
  return MutID

def epistasiscal(var1fit,var2fit,dfit):
  if var1fit < 0.01 or var2fit < 0.01: cap = 'NoNeg'
  else: cap ='All'
  absepi = float(dfit) - floor(var1fit*var2fit) #Absolute Epistasis model
  relepi = float(floor(dfit))/floor(var1fit*var2fit) #Relative Epistasis model (default)
  if cap == 'NoNeg' and relepi < 1: relepi = 1
  else: relepi = relepi
  return absepi, relepi

def ThreeWayAnalysis(muts,WT,Index2Pos,Pos2Index,fithash,condition,outfile):
  countmut = 0
  for mut in muts:
    countmut += 1
    if countmut%10000 == 0: print 'Processed %d variants' % countmut
    if hamming(mut,WT) <= 1:
      MutPosList = []
      [MutPosList.append(n) for n in range(len(WT)) if WT[n] == mut[n]]
      MutPosList = combinations(MutPosList,3)
      mutfit = float(fithash[mut][condition])
      for MutPos in MutPosList: 
        assert(len(MutPos)==3)
        MutIDList = GenMutfromPos(sorted(list(MutPos)),Index2Pos,WT)
        for MutID in MutIDList:
          MutID = MutID.rsplit('-')
          for M1 in MutID:
            MutIDtmp = MutID[:]
            MutIDtmp.remove(M1)
            M2 = MutIDtmp[0]
            M3 = MutIDtmp[1]
            M1ID = mut2index(M1,Pos2Index,mut)
            M2ID = mut2index(M2,Pos2Index,mut)
            M3ID = mut2index(M3,Pos2Index,mut)
            M12ID = mut2index('-'.join([M1,M2]),Pos2Index,mut)
            M13ID = mut2index('-'.join([M1,M3]),Pos2Index,mut)
            M23ID = mut2index('-'.join([M2,M3]),Pos2Index,mut)
            M123ID = mut2index('-'.join([M1,M2,M3]),Pos2Index,mut)
            checkpoint = 'PASS'
            for m in [M1ID, M2ID, M3ID, M12ID, M13ID, M23ID, M123ID]: 
              if not fithash.has_key(m): checkpoint = 'NA'
            if checkpoint == 'NA': continue
            M1fit   = floor(float(fithash[M1ID][condition]))/floor(mutfit)
            M2fit   = floor(float(fithash[M2ID][condition]))/floor(mutfit)
            M3fit   = floor(float(fithash[M3ID][condition]))/floor(mutfit)
            M12fit  = floor(float(fithash[M12ID][condition]))/floor(mutfit)
            M13fit  = floor(float(fithash[M13ID][condition]))/floor(mutfit)
            M23fit  = floor(float(fithash[M23ID][condition]))/floor(mutfit)
            M123fit = floor(float(fithash[M123ID][condition]))/floor(mutfit)
            epi12abs, epi12rel     = epistasiscal(M1fit,M2fit,M12fit)
            epi13abs, epi13rel     = epistasiscal(M1fit,M3fit,M13fit)
            epi1_23abs, epi1_23rel = epistasiscal(M1fit,M23fit,M123fit)
            outfile.write("\t".join(map(str,[mut, mutfit, M1, M2, M3, 
                                             fithash[M1ID][condition], fithash[M2ID][condition], fithash[M3ID][condition],
                                             fithash[M12ID][condition], fithash[M13ID][condition], fithash[M23ID][condition],
                                             fithash[M123ID][condition], epi12rel, epi13rel, epi1_23rel]))+"\n")

def main():
  fitfile     = 'result/Mutfit'
  outfile     = 'analysis/ThreeWayEpi'
  WT          = 'VDGV'
  Index2Pos   = {0:39,1:40,2:41,3:54}
  Pos2Index   = {39:0,40:1,41:2,54:3}
  condition   = 'I20fit'
  fithash     = TsvWithHeader2Hash(fitfile)
  fithash     = filterfithash(fithash)
  muts        = fithash.keys()
  outfile     = open(outfile,'w')
  header      = "\t".join(['Mut','Mutfit','M1','M2','M3','M1fit','M2fit','M3fit',
                           'M12fit','M13fit','M23fit','M123fit','Epi1vs2','Epi1vs3','Epi1vs23'])
  outfile.write(header+"\n")
  ThreeWayAnalysis(muts,WT,Index2Pos,Pos2Index,fithash,condition,outfile)
  outfile.close()

if __name__ == '__main__':
  main()
