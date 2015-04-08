#!/usr/bin/python
import operator
from itertools import imap

def adjpos(muts):
  muts = muts.rsplit('-')
  newmuts = []
  for mut in muts:
    newmuts.append(mut[0]+str(int(mut[1:-1])+1)+mut[-1])
  return '-'.join(newmuts)

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

def callmut(var1, var2, Index2pos): #mutation for converting var1 to var2
  mut = []
  assert(len(var1)==len(var2))
  for i in range(len(var1)):
    if var1[i] != var2[i]:
      mut.append(var1[i]+str(Index2pos[i])+var2[i])
  if len(mut) == 0: return 'WT'
  else: return '-'.join(mut)

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

def hamming(str1, str2):
    assert len(str1) == len(str2)
    return sum(imap(operator.ne, str1, str2))

def callfit(mut,condition,fithash,WT,HD):
  if hamming(mut,WT) != int(HD): return 'WT'
  if fithash.has_key(mut): return fithash[mut][condition]
  else: return 'NA'

def formatting(muts,condition,fithash,AOhash,outfile,Index2Pos,WT):
  outfile1 = open(outfile+'HD1','w')
  outfile2 = open(outfile+'HD2','w')
  outfile3 = open(outfile+'HD3','w')
  outfile4 = open(outfile+'HD4','w')
  header1  = "\t".join(['mut','ID1','HD','Input','fit',
                        'AOfit'])
  header2  = "\t".join(['mut','ID1','ID2','HD','Input','fit',
                        'S1','S2','AO1fit','AO2fit','AODfit'])
  header3  = "\t".join(['mut','ID1','ID2','ID3','HD','Input','fit',
                        'S1','S2','S3','D12','D13','D23'])
  header4  = "\t".join(['mut','ID1','ID2','ID3','ID4','HD','Input','fit',
                        'S1','S2','S3','S4','D12','D13','D14','D23','D24','D34','T123','T124','T134','T234'])
  outfile1.write(header1+"\n")
  outfile2.write(header2+"\n")
  outfile3.write(header3+"\n")
  outfile4.write(header4+"\n")
  for mut in muts:
    ID    = callmut(WT, mut, Index2Pos)
    HD    = fithash[mut]['HD']
    Input = fithash[mut]['Input']
    fit   = fithash[mut][condition]
    S1    = callfit(mut[0]+  WT[1]+  WT[2]+  WT[3], condition,fithash,WT,1)
    S2    = callfit( WT[0]+ mut[1]+  WT[2]+  WT[3], condition,fithash,WT,1)
    S3    = callfit( WT[0]+  WT[1]+ mut[2]+  WT[3], condition,fithash,WT,1)
    S4    = callfit( WT[0]+  WT[1]+  WT[2]+ mut[3], condition,fithash,WT,1)
    D12   = callfit(mut[0]+ mut[1]+  WT[2]+  WT[3], condition,fithash,WT,2)
    D13   = callfit(mut[0]+  WT[1]+ mut[2]+  WT[3], condition,fithash,WT,2)
    D14   = callfit(mut[0]+  WT[1]+  WT[2]+ mut[3], condition,fithash,WT,2)
    D23   = callfit( WT[0]+ mut[1]+ mut[2]+  WT[3], condition,fithash,WT,2)
    D24   = callfit( WT[0]+ mut[1]+  WT[2]+ mut[3], condition,fithash,WT,2)
    D34   = callfit( WT[0]+  WT[1]+ mut[2]+ mut[3], condition,fithash,WT,2)
    T123  = callfit(mut[0]+ mut[1]+ mut[2]+  WT[3], condition,fithash,WT,3)
    T124  = callfit(mut[0]+ mut[1]+  WT[2]+ mut[3], condition,fithash,WT,3)
    T134  = callfit(mut[0]+  WT[1]+ mut[2]+ mut[3], condition,fithash,WT,3)
    T234  = callfit( WT[0]+ mut[1]+ mut[2]+ mut[3], condition,fithash,WT,3)
    if int(HD) == 4:
      outfit = [fit, S1, S2, S3, S4, D12, D13, D14, D23, D24, D34, T123, T124, T134, T234]
      while 'WT' in outfit: outfit.remove('WT')
      assert(len(outfit)==15)
      output = "\t".join(map(str,[mut,ID.replace('-',"\t"),HD,Input]))+"\t"+"\t".join(map(str,outfit))
      outfile4.write(output+"\n")
    if int(HD) == 3:
      outfit = [fit, S1, S2, S3, S4, D12, D13, D14, D23, D24, D34]
      while 'WT' in outfit: outfit.remove('WT')
      assert(len(outfit)==7)
      output = "\t".join(map(str,[mut,ID.replace('-',"\t"),HD,Input]))+"\t"+"\t".join(map(str,outfit))
      outfile3.write(output+"\n")
    if int(HD) == 2: 
      outfit = [fit, S1, S2, S3, S4]
      while 'WT' in outfit: outfit.remove('WT')
      assert(len(outfit)==3)
      M1 = convertmut(ID.rsplit('-')[0])
      M2 = convertmut(ID.rsplit('-')[1])
      if AOhash.has_key(mut): AODfit = AOhash[mut]
      else: AODfit = 'NA'
      output = "\t".join(map(str,[mut,ID.replace('-',"\t"),HD,Input]))+"\t"+"\t".join(map(str,outfit))+"\t"+"\t".join(map(str,[AOhash[M1],AOhash[M2],AODfit]))
      outfile2.write(output+"\n")
    if int(HD) == 1: 
      outfit = [fit]
      while 'WT' in outfit: outfit.remove('WT')
      assert(len(outfit)==1)
      output = "\t".join(map(str,[mut,ID.replace('-',"\t"),HD,Input]))+"\t"+"\t".join(map(str,outfit))+"\t"+"\t".join(map(str,[AOhash[mut]]))
      outfile1.write(output+"\n")
  outfile1.close()
  outfile2.close()
  outfile3.close()
  outfile4.close()
       
      

def main():
  WT         = 'VDGV'
  fitfile    = 'result/Mutfit'
  outfile    = 'transfer/NW2AOset'
  condition  = 'I20fit'
  Index2Pos  = {0:39,1:40,2:41,3:54}
  fithash    = TsvWithHeader2Hash(fitfile)
  fithash    = filterfithash(fithash)
  muts       = fithash.keys()
  AOhash  = {}
  AOhash  = hashin('/home/wchnicholas/RSun/Robust/Doc/SMutList', AOhash)
  AOhash  = hashin('/home/wchnicholas/RSun/Robust/Doc/DMutList', AOhash)
  formatting(muts,condition,fithash,AOhash,outfile,Index2Pos,WT)

if __name__ == '__main__':
  main()

