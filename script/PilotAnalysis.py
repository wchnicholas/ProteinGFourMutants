#!/usr/bin/python
import os
import sys
import glob
from collections import Counter

def epihashing(fitfile):
  H = {}
  infile = open(fitfile,'r')
  countline = 0
  for line in infile.xreadlines():
    countline += 1
    if countline == 1: continue
    line = line.rstrip().rsplit("\t")
    mut = line[0]
    epi = line[1]
    H[mut] = float(epi)
  infile.close()
  return H

def candidating(epicand):
  H = {}
  epicand = [str(int(cand.rsplit('-')[0][1:-1])+1)+'-'+str(int(cand.rsplit('-')[1][1:-1])+1) for cand in epicand]
  epicand = Counter(epicand)
  return epicand

def drawgraph(outfile,candhash):
  outfile=open(outfile,'w')
  outfile.write('graph{'+"\n"+"\t"+'layout=dot'+
                         "\n"+"\t"+'rankdir=LR'+
                         "\n"+"\t"+'node [shape=box]'+"\n")
  nodes = []
  [nodes.extend(cand.rsplit('-')) for cand in candhash.keys()]
  nodes = list(set(nodes))
  for var in nodes:
    outfile.write("\t"+var+' [fillcolor="white", color=black, style="filled,rounded", width=0.5];'+"\n")
  for cand in candhash.keys():
    var1 = cand.rsplit('-')[0]
    var2 = cand.rsplit('-')[1]
    for n in range(candhash[cand]):
      outfile.write("\t"+var1+'--'+var2+' [color=black];'+"\n")
  outfile.write('}'+"\n")
  outfile.close()

def main():
  infile   = 'doc/Epistasis'
  outfile  = 'xdot/Pilot.dot'
  epihash  = epihashing(infile)
  epicand  = sorted(epihash.keys(),key=lambda x:epihash[x])[::-1]
  candhash = candidating(epicand[0:20])
  print candhash
  #drawgraph(outfile,candhash)
  #os.system('dot -Tpng %s -o %s' % (outfile,outfile.replace('.dot','.png')))

if __name__ == '__main__':
  main()
