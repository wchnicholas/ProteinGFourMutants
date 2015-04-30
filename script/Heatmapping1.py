#!/usr/bin/python
#Formatting the heatmap for the EpiSD and EpiRange for individual double-mutations
import os
import sys

def TsvWithHeader2Hash(infile):
  H = {}
  infile = open(infile,'r')
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

def gendmutlist(WT, Index2pos):
  dmutlist = []
  smutlist = []
  aas = ['E','D','R','K','H','Q','N','S','T','P','G','C','A','V','I','L','M','F','Y','W']
  for i in range(len(WT)):
    for aa1 in aas:
      if aa1 != WT[i]:
        m1 = WT[i]+str(Index2pos[i])+aa1
        smutlist.append(m1)
        for j in range(len(WT)):
          for aa2 in aas:
            if aa2 != WT[j] and i < j:
              m2 = WT[j]+str(Index2pos[j])+aa2
              dmutlist.append(m1+'-'+m2)
  return smutlist, dmutlist

def heatmapping1(epihash,smutlist,outfile1,outfile2):
  blankvalue = 99
  outfile1 = open(outfile1,'w')
  outfile2 = open(outfile2,'w')
  header  = 'mut'+"\t"+"\t".join(smutlist)
  outfile1.write(header+"\n")
  outfile2.write(header+"\n")
  for i in range(len(smutlist)):
    out1 = [smutlist[i]]
    out2 = [smutlist[i]]
    for j in range(len(smutlist)):
      m1   = smutlist[i]
      m2   = smutlist[j]
      pos1 = int(m1[1:-1])
      pos2 = int(m2[1:-1])
      if pos1 == pos2: out1.append(blankvalue); out2.append(blankvalue);continue
      elif pos1 < pos2: dmut = m1+'-'+m2
      elif pos1 > pos2: dmut = m2+'-'+m1
      if i == j: out1.append(blankvalue); out2.append(blankvalue); print 'Something is wrong'; continue
      else: out1.append(epihash[dmut]['EpiRange']); out2.append(epihash[dmut]['EpiSD']);
    outfile1.write("\t".join(map(str,out1))+"\n")
    outfile2.write("\t".join(map(str,out2))+"\n")
  outfile1.close()
  outfile2.close()
            
def main():
  epifile    = 'analysis/EpiDiffBGI20fit'
  outfile1   = 'analysis/HeatMapEpiRange'
  outfile2   = 'analysis/HeatMapEpiSD'
  WT         = 'VDGV'
  Index2pos  = {0:39,1:40,2:41,3:54}
  pos2index  = {39:0,40:1,41:2,54:3}
  epihash    = TsvWithHeader2Hash(epifile)
  smutlist, dmutlist   = gendmutlist(WT, Index2pos)
  heatmapping1(epihash,smutlist,outfile1,outfile2)

if __name__ == '__main__':
  main()
