#!/usr/bin/python
import os
import sys
import glob
from math import log

def HashingRepoducibility(infile):
  infile = open(infile,'r')
  rephash = {}
  for line in infile.xreadlines():
    if 'mut' in line: continue
    line = line.rstrip().rsplit("\t")
    mut    = line[0]
    fit    = float(line[1])
    if fit < 1: continue
    count  = map(float,line[2:-1])
    freq   = map(lambda x:x/sum(count), count)
    freq   = filter(lambda x:x!=0, freq)
    rep    = sum(map(lambda x:-x*log(x), freq))
    rephash[mut] = {'fit':fit,'rep':rep}
  infile.close()
  return rephash

def TrajRepTrend(infile, outfile, rephash):
  infile  = open(infile,'r')
  for line in infile.xreadlines():
    line  = line.rstrip().rsplit("\t")
    mut   = line[0]
    leng  = line[1]
    path  = line[2]
    if not rephash.has_key(mut): continue
    steps = path.rsplit('->')
    steps.remove('TheEnd')
    reps  = map(lambda x:rephash[x]['rep'], steps)
    reps  = reps+map(str,[-1]*(20-int(leng)))
    outfile.write(mut+"\t"+leng+"\t"+"\t".join(map(str,reps))+"\n")
  infile.close()

def main():
  destfile  = "analysis/LocalMaxDes_weight"
  outfile   = "analysis/RepTraj_weight"
  trajfiles = sorted(glob.glob("simulations/weight/LocalMaxClimb_weight*"))
  rephash   = HashingRepoducibility(destfile)
  outfile = open(outfile,'w')
  for trajfile in trajfiles:
    print trajfile
    TrajRepTrend(trajfile, outfile, rephash)
  outfile.close()
  
if __name__ == "__main__":
  main()
