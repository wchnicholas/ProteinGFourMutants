#!/usr/bim/python
import os
import sys
import random
import operator
import numpy as np
import networkx as nx
from math import exp
from itertools import imap

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

def findlocalmax(muts, fithash, condition):
  localmax_count = 0
  count    = 0
  process  = 0
  localmax_muts = {}
  for mut in muts:
    process += 1
    if process%10000 == 0: print 'Finished %s localmax classificiation' % process
    mutfit   = float(fithash[mut][condition])
    variants = genvarsneighbor(mut)
    neighfit = []
    for var in variants:
      if fithash.has_key(var): neighfit.append(float(fithash[var][condition])); 
      else: neighfit.append('NA'); print 'problematic mutant: %s' % mut; break
    if 'NA' not in neighfit:
      assert(len(neighfit)==76)
      count += 1
      if max(neighfit) <= mutfit: 
        localmax_count += 1
        localmax_muts[mut] = str(mutfit)
  return count, localmax_count, localmax_muts

def climbing_greedy(mut,fithash,condition):
  mutfit   = float(fithash[mut][condition])
  neighmaxfit = float(mutfit)
  neighmaxmut = 'TheEnd'
  variants = genvarsneighbor(mut)
  assert(len(variants)==76)
  for var in variants:
    neighfit = float(fithash[var][condition])
    if neighfit > neighmaxfit: 
      neighmaxfit = neighfit
      neighmaxmut = var
  return neighmaxmut

def climbing_random(mut,fithash,condition):
  mutfit   = float(fithash[mut][condition])
  variants = genvarsneighbor(mut)
  assert(len(variants)==76)
  neighmaxmuts = [var for var in variants if float(fithash[var][condition]) > mutfit]
  if len(neighmaxmuts) >= 1: return random.choice(neighmaxmuts) 
  else: return 'TheEnd'
  
def climbing_weight(mut,fithash,condition):
  mutfit   = float(fithash[mut][condition])
  variants = genvarsneighbor(mut)
  assert(len(variants)==76)
  neighhash = {}
  for var in variants: 
    varfit = float(fithash[var][condition])
    if varfit > mutfit: neighhash[var] = varfit-mutfit
  if len(neighhash.keys()) >= 1:
    normscale = sum(neighhash.values())
    neighs = neighhash.keys()
    probs  = [(float(neighhash[var])/float(normscale)) for var in neighs]
    return np.random.choice(neighs,p=probs)
  else: return 'TheEnd' 
  

def pathtomax(muts,fithash,condition,outfile,climbtype):
  outfile  = open(outfile,'w')
  header   = "\t".join(['mut','steps','path','localmax'])
  outfile.write(header+"\n")
  process  = 0
  assert(len(muts)==160000)
  for mut in muts: 
    process += 1
    if process%10000 == 0: print 'Finished %d mutants: %s search for local max' % (process, climbtype)
    step  = mut
    path  = [mut]
    while step != 'TheEnd': 
      if   climbtype == 'greedy': step = climbing_greedy(step,fithash,condition)
      elif climbtype == 'random': step = climbing_random(step,fithash,condition)
      elif climbtype == 'weight': step = climbing_weight(step,fithash,condition)
      else: print 'Bad Option for ClimbType'; sys.exit()
      path.append(step)
    outfile.write("\t".join(map(str,[mut,len(path)-2,'->'.join(path),path[-2]]))+"\n")
  outfile.close()

def main():
  seed        = sys.argv[1]
  climbtype   = 'weight' #"random", "greedy" or "weight"
  WT          = 'VDGV'
  fitfile     = 'result/Mutfit'
  missfitfile = 'result/regression_missing'
  climbout    = 'sim/LocalMaxClimb_'+climbtype+'_'+seed
  fcutoff     = -1
  condition   = 'I20fit'
  Index2pos   = {0:39,1:40,2:41,3:54}
  fithash     = TsvWithHeader2Hash(fitfile)
  fithash     = filterfithash(fithash,condition,fcutoff)
  fithash     = fillinmissing(fithash,missfitfile,condition) 
  muts        = fithash.keys()
  Varwithallneighbors, localmax_count, localmax_muts = findlocalmax(muts, fithash, condition)
  pathtomax(muts,fithash,condition,climbout,climbtype)

if __name__ == '__main__':
  main()
