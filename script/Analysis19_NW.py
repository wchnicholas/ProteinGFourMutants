#!/usr/bin/python
import sys
import glob
import operator
from itertools import imap

def hamming(str1, str2):
    assert len(str1) == len(str2)
    return sum(imap(operator.ne, str1, str2))

def classifysteps(infile,classhash):
  infile = open(infile,'r')
  for line in infile.xreadlines():
    if 'mut' in line: continue
    line = line.rstrip().rsplit("\t")
    path = line[2] #COLUMN THAT CONTAINS THE PATH INFORMATION
    #path = line[4] #WHEN INPUT FILE = analysis/LocalMaxEvolvePotWT
    if path == 'NA': continue
    path = path.replace('->TheEnd','').rsplit('->')
    end  = path[-1]
    ptype       = [0,0,0]
    for n in range(1,len(path)):
      dist_before = hamming(path[n-1],end)
      dist_after  = hamming(path[n],end)
      if dist_after-dist_before == -1:  classhash['Forward'] += 1; ptype[0] = 0 
      elif dist_after-dist_before == 0: classhash['Neutral'] += 1; ptype[1] = 1
      elif dist_after-dist_before == 1: classhash['Reverse'] += 1; ptype[2] = 2
    if sum(ptype) == 3 and hamming(path[0],end)==4: print path
  infile.close()
  return classhash

def main():
  #infiles   = glob.glob('simulations/weight/LocalMaxClimb_*')
  infiles   = glob.glob('analysis/LocalMaxClimb_greedy')
  #infiles   = ['analysis/LocalMaxEvolvePotWT']
  print "Total # of files = %d" % len(infiles)
  classhash = {'Forward':0, 'Neutral':0, 'Reverse': 0}
  countfile = 0
  for infile in infiles:
    countfile += 1
    if countfile%10 == 0: print "Processed %d files" % countfile
    classhash = classifysteps(infile,classhash)
  print classhash
  #RESULTS COPIED IN analysis/StepsInfo
 
if __name__ == '__main__':
  main()
