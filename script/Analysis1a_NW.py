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
  for mut in muts:
    if fithash.has_key(mut): fits.append(float(fithash[mut][condition]))
    else: return 'NA'
  if min(fits) < 0.01: return float(0.01)
  else: 
    expfit = float(1)
    for fit in fits: expfit = fit*expfit
    return expfit

def epical(fithash,muts,WT,condition,outfile):
  print 'Generating file: %s' % outfile
  outfile = open(outfile,'w')
  header = ['mut','mutfit','Sexpfit','Dexpfit_12vs34','Dexpfit_13vs24','Dexpfit_23vs14',
            'Texpfit_1vs234','Texpfit_2vs134','Texpfit_3vs124','Texpfit_4vs123',
            'Mexpfit_12vs3vs4','Mexpfit_13vs2vs4','Mexpfit_14vs2vs3',
            'Mexpfit_23vs1vs4','Mexpfit_24vs1vs3','Mexpfit_34vs1vs2','Minexpfit','Maxexpfit']
  outfile.write("\t".join(header)+"\n")
  for mut in muts:
    HD = int(fithash[mut]['HD'])
    mutfit = float(fithash[mut][condition])
    if HD == 4:
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
      Sexpfit = expectingfit(fithash, condition, [S1,S2,S3,S4])
      Dexpfit_12vs34   = expectingfit(fithash, condition, [D12,D34])
      Dexpfit_13vs24   = expectingfit(fithash, condition, [D13,D24])
      Dexpfit_23vs14   = expectingfit(fithash, condition, [D23,D14])
      Texpfit_1vs234   = expectingfit(fithash, condition, [S1,T234])
      Texpfit_2vs134   = expectingfit(fithash, condition, [S2,T134])
      Texpfit_3vs124   = expectingfit(fithash, condition, [S3,T124])
      Texpfit_4vs123   = expectingfit(fithash, condition, [S4,T123])
      Mexpfit_12vs3vs4 = expectingfit(fithash, condition, [D12,S3,S4])
      Mexpfit_13vs2vs4 = expectingfit(fithash, condition, [D13,S2,S4])
      Mexpfit_14vs2vs3 = expectingfit(fithash, condition, [D14,S2,S3])
      Mexpfit_23vs1vs4 = expectingfit(fithash, condition, [D23,S1,S4])
      Mexpfit_24vs1vs3 = expectingfit(fithash, condition, [D24,S1,S3])
      Mexpfit_34vs1vs2 = expectingfit(fithash, condition, [D34,S1,S2])

      outfit = [Sexpfit,Dexpfit_12vs34,Dexpfit_13vs24,Dexpfit_23vs14,
                Texpfit_1vs234,Texpfit_2vs134,Texpfit_3vs124,Texpfit_4vs123,
                Mexpfit_12vs3vs4,Mexpfit_13vs2vs4,Mexpfit_14vs2vs3,
                Mexpfit_23vs1vs4,Mexpfit_24vs1vs3,Mexpfit_34vs1vs2]
      maxexp = str(floor(max(outfit)))
      minexp = str(floor(min(outfit)))
      if 'NA' not in outfit: outfile.write(mut+"\t"+str(floor(mutfit))+"\t"+"\t".join(map(str,map(floor,outfit)))+"\t"+minexp+"\t"+maxexp+"\n")
  outfile.close()
   
def main():
  WT         = 'VDGV'
  fitfile    = 'result/Mutfit'
  out10      = 'result/HD4EpiIGG10'
  out20      = 'result/HD4EpiIGG20'
  out90      = 'result/HD4EpiIGG90'
  fithash    = TsvWithHeader2Hash(fitfile)
  fithash    = filterfithash(fithash)
  muts       = fithash.keys()
  epical(fithash,muts,WT,'I10fit',out10)
  epical(fithash,muts,WT,'I20fit',out20)
  epical(fithash,muts,WT,'I90fit',out90)

if __name__ == '__main__':
  main()
