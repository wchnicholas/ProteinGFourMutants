#!/usr/bin/python
import os
import sys
import glob
import operator
import numpy as np
from math import exp
from itertools import imap

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

def hamming(str1, str2):
    assert len(str1) == len(str2)
    return sum(imap(operator.ne, str1, str2))

def filterfithash(fithash,condition,fcutoff):
  for mut in fithash.keys():
    if '_' in mut or fithash[mut][condition] == 'NA' or float(fithash[mut][condition]) < float(fcutoff):
      del fithash[mut]
  return fithash

def fillinmissing(fithash,missfitfile,condition):
  infile = open(missfitfile,'r')
  for line in infile.xreadlines():
    if 'genotype_' in line: continue
    line = line.rstrip().rsplit("\t")
    mut  = line[0]
    fit  = exp(float(line[1]))
    fithash[mut] = {}
    fithash[mut][condition] = fit
  infile.close()
  return fithash

def desthashing(probhash,pathhash,lenghash,filename):
  infile = open(filename,'r')
  for line in infile.xreadlines():
    if 'mut' in line: continue 
    line = line.rstrip().rsplit("\t")
    mut  = line[0] 
    leng = line[1]
    path = line[2].replace('->TheEnd','')
    des  = line[3]
    #lenghash[mut].append(float(leng))
    probhash[mut][des] += 1 
    if pathhash[mut].has_key(path): pathhash[mut][path] += 1
    else: pathhash[mut][path] = 1
  infile.close()
  return probhash, pathhash, lenghash

def genprobhash(localmaxs):
  probhash = {}
  pathhash = {}
  lenghash = {}
  aas = ['E','D','R','K','H','Q','N','S','T','P','G','C','A','V','I','L','M','F','Y','W']
  for aa1 in aas: 
    for aa2 in aas:
      for aa3 in aas:
        for aa4 in aas:
          mut = aa1+aa2+aa3+aa4
          probhash[mut] = {}
          pathhash[mut] = {}
          lenghash[mut] = []
          for localmax in localmaxs: 
            probhash[mut][localmax] = 0
  return probhash, pathhash, lenghash

def pathdistancecal(Steps_path1,Steps_path2, dists):
  for step_path1 in Steps_path1:
    mindist = 4
    for step_path2 in Steps_path2:
      stepdist = hamming(step_path1,step_path2)
      if stepdist < mindist: mindist = stepdist
    dists.append(mindist)
  return dists

def pathdiversitycal(pathhash, mut):
  paths = sorted(pathhash[mut].keys())
  diversity   = 0
  for i in range(len(paths)):
    for j in range(len(paths)):
      if i != j:
        path1 = paths[i]
        path2 = paths[j]
        Prob_path1  = float(pathhash[mut][path1])/float(sum(pathhash[mut].values()))
        Prob_path2  = float(pathhash[mut][path2])/float(sum(pathhash[mut].values()))
        Steps_path1 = path1.rsplit('->')
        Steps_path2 = path2.rsplit('->')
        dists   = []
        dists   = pathdistancecal(Steps_path1,Steps_path2, dists)
        dists   = pathdistancecal(Steps_path2,Steps_path1, dists)
        avgdist = float(np.mean(dists))
        diversity += Prob_path1*Prob_path2*avgdist
  return diversity
        
def compileout(probhash,pathhash,lenghash,outf,lmaxs,fithash,condition):
  outfile = open(outf,'w')
  outfile.write('mut'+"\t"+'fit'+"\t"+"\t".join(lmaxs)+"\t"+'PathReprod'+"\n")
  for mut in probhash.keys():
    out = [mut,fithash[mut][condition]]
    for lmax in lmaxs:
      out.append(probhash[mut][lmax])
    out.append(sum([(float(rep)/float(sum(pathhash[mut].values())))**2 for rep in pathhash[mut].values()]))
    #out.append(pathdiversitycal(pathhash, mut))
    outfile.write("\t".join(map(str,out))+"\n")
  print 'done'
  outfile.close()

def main():
  simtype = 'weight' #weight or random or greedy
  lmaxs = [l.rstrip().rsplit("\t")[0] for l in open('analysis/LocalMaxMuts','r').readlines()]
  files = sorted(glob.glob('simulations/'+simtype+'/LocalMaxClimb_'+simtype+'*'))
  #files = sorted(glob.glob('analysis/LocalMaxClimb_greedy'))
  outf  = 'analysis/LocalMaxDes_'+simtype
  fitfile     = 'result/Mutfit'
  missfitfile = 'result/regression_missing'
  #missfitfile = 'result/regression_all_WT'
  fcutoff     = -1
  condition   = 'I20fit'
  Index2pos   = {0:39,1:40,2:41,3:54}
  fithash     = TsvWithHeader2Hash(fitfile)
  print "Total # of variants in the raw data: %d" % len(fithash.keys())
  fithash     = filterfithash(fithash,condition,fcutoff)
  fithash     = fillinmissing(fithash,missfitfile,condition)
  muts        = fithash.keys()
  print "Total # of variants after fill in with regression: %d" % len(muts)

  probhash, pathhash, lenghash = genprobhash(lmaxs)
  print 'A hash of %d variants with %d local max is created' % (len(probhash.keys()),len(lmaxs))
  counttrial = 0
  for filename in files: 
    counttrial += 1
    if counttrial%100 == 0: print "Finished parsing %d files" % counttrial
    probhash, pathhash, lenghash = desthashing(probhash,pathhash,lenghash,filename)
  compileout(probhash,pathhash,lenghash,outf,lmaxs,fithash,condition)
  
if __name__ == '__main__':
  main()
