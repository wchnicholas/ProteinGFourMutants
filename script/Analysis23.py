#!/usr/bin/python
import os
import sys
import operator
from itertools import imap

def hamming(str1, str2):
    assert len(str1) == len(str2)
    return sum(imap(operator.ne, str1, str2))

def peakdictinit(peaks):
  H = {}
  for peak in peaks:
    H[peak] = [0,0] #[Count of all accessible nodes, Count of accessible nodes by direct paths]
  return H

def accessdictinit(peaknum):
  H = {}
  for i in range(peaknum+1):
    H[i] = [0,0] #[All paths, direct paths only]
  return H

def dirpathing(mut, path2peak):
  l = []
  for p, s in path2peak:
    if s == hamming(mut,p):
      l.append('dir')
    elif s > hamming(mut,p):
      l.append('indir')
    elif s == -1: 
      l.append('na')
    else: print "ERROR: Something weird about the path lengths"; sys.exit()
  return l

def pathanalysis(infile):
  infile = open(infile,'r')
  peaks  = []
  peakdict   = {}
  accessdict = {}
  for line in infile.xreadlines():
    line  = line.rstrip().rsplit("\t")
    if 'mut' == line[0]:
      peaks      = line[1:-1]
      peakdict   = peakdictinit(peaks)
      accessdict = accessdictinit(len(peakdict.keys()))
      continue
    else: 
      mut = line[0]
      pathlengths = map(int, line[1:-1])
      path2peak = zip(peaks, pathlengths)
      dirpaths  = dirpathing(mut, path2peak)
      for p, c in zip(peaks, dirpaths):
        if c == 'dir' or c == 'indir': peakdict[p][0] += 1
        if c == 'dir': peakdict[p][1] += 1
      numaccesspeaks_all = dirpaths.count('dir')+dirpaths.count('indir')
      numaccesspeaks_dir = dirpaths.count('dir')
      accessdict[numaccesspeaks_all][0]+=1
      accessdict[numaccesspeaks_dir][1]+=1
  infile.close()
  return accessdict, peakdict
    

def main():
  infile = 'analysis/LocalMaxPathLen'
  accessdict, peakdict = pathanalysis(infile)
  print "\t".join(['NumOfPeaks','Count_AccessibleFromNodes','Count_AccessibleFromNodesWithDirPath'])
  for numpeak in accessdict:
    print "\t".join(map(str,[numpeak, accessdict[numpeak][0], accessdict[numpeak][1]]))
  
  print "\t".join(['Peaks','Count_AccessibleFromNodes','Count_AccessibleFromNodesWithDirPath'])
  for peak in peakdict:
    print "\t".join(map(str,[peak, peakdict[peak][0], peakdict[peak][1]]))

if __name__ == "__main__":
  main()
