#!/usr/bim/python
import os
import sys
import random
import operator
import numpy as np
import networkx as nx
import matplotlib.pyplot as plt
from math import exp
from itertools import imap
from collections import Counter

def all_shortest_paths(G, source, target, weight=None):
    if weight is not None:
        pred,dist = nx.dijkstra_predecessor_and_distance(G,source,weight=weight)
    else:
        pred = nx.predecessor(G,source)
    if target not in pred:
        raise nx.NetworkXNoPath()
    stack = [[target,0]]
    top = 0
    while top >= 0:
        node,i = stack[top]
        if node == source:
            yield [p for p,n in reversed(stack[:top+1])]
        if len(pred[node]) > i:
            top += 1
            if top == len(stack):
                stack.append([pred[node][i],0])
            else:
                stack[top] = [pred[node][i],0]
        else:
            stack[top-1][1] += 1
            top -= 1

def connected_components(G):
    if G.is_directed():
        raise nx.NetworkXError("""Not allowed for directed graph G.
              Use UG=G.to_undirected() to create an undirected graph.""")
    seen={}
    components=[]
    for v in G:      
        if v not in seen:
            c=nx.single_source_shortest_path_length(G,v)
            components.append(list(c.keys()))
            seen.update(c)
    components.sort(key=len,reverse=True)            
    return components 

def hamming(str1, str2):
    assert len(str1) == len(str2)
    return sum(imap(operator.ne, str1, str2))

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
    #if mut[1] != 'D': continue ####################
    H[mut] = {}
    for i in range(1,len(line)): H[mut][header[i]] = line[i]
  infile.close()
  return H

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
    #if mut[1] != 'D': continue ####################
    fit  = exp(float(line[1]))
    fithash[mut] = {}
    fithash[mut][condition] = fit
  infile.close()
  return fithash

def buildgraph(muts):
  G=nx.Graph()
  for m1 in muts:  G.add_node(m1)
  for m1 in muts:
    for m2 in muts:
      if hamming(m1,m2) == 1: G.add_edge(m1,m2)
  return G

def genvarsneighbor(var):
  variants = []
  aas = ['E','D','R','K','H','Q','N','S','T','P','G','C','A','V','I','L','M','F','Y','W']
  [variants.append(aa+var[1]+var[2]+var[3]) for aa in aas]
  [variants.append(var[0]+aa+var[2]+var[3]) for aa in aas]
  [variants.append(var[0]+var[1]+aa+var[3]) for aa in aas]
  [variants.append(var[0]+var[1]+var[2]+aa) for aa in aas]
  while var in variants: variants.remove(var)
  assert(len(variants)==76)
  return variants

def varfitdist(muts, fithash, condition,outfile):
  outfile = open(outfile,'w')
  process = 0
  for mut in muts:
    process += 1
    if process%10000 == 0: print 'Finished step fitness computing initiating from %s variants' % process
    mutfit   = float(fithash[mut][condition])
    variants = genvarsneighbor(mut)
    neighfit = []
    for var in variants:
      varfit = float(fithash[var][condition])
      outfile.write("\t".join(map(str,sorted([mut,var])))+"\t"+str(abs(varfit-mutfit))+"\n")
  outfile.close()

  

def main():
  climbtype   = 'greedy' #"random", "greedy" or "weight"
  WT          = 'VDGV'
  fitfile     = 'result/Mutfit'
  missfitfile = 'result/regression_missing'
  outfile     = 'analysis/FitStepDist'
  fcutoff     = -1
  condition   = 'I20fit'
  Index2pos   = {0:39,1:40,2:41,3:54}
  fithash     = TsvWithHeader2Hash(fitfile)
  print "Total # of variants in the raw data: %d" % len(fithash.keys())
  fithash     = filterfithash(fithash,condition,fcutoff)
  print "Total # of variants pass filter of raw data: %d" % len(fithash.keys())
  fithash     = fillinmissing(fithash,missfitfile,condition) 
  muts        = fithash.keys()
  print "Total # of variants after fill in with regression: %d" % len(muts)
  print "# of mutant pass cutoff: %d" % len(muts)
  varfitdist(muts, fithash, condition,outfile)

if __name__ == '__main__':
  main()
