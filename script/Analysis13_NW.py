#!/usr/bim/python
import sys
import matplotlib.pyplot as plt
import networkx as nx
import operator
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
  return variants

def buildDigraph(G,muts,fithash,condition):
  for mut in muts:
    mutfit = float(fithash[mut][condition])
    G.add_node(mut)
    variants = genvarsneighbor(mut)
    assert(len(variants)==76)
    [G.add_edge(mut,var) for var in variants if float(fithash[var][condition]) > mutfit]
  return G

def searchingpath(G,localmaxs,muts,outfile):
  outfile = open(outfile,'w')
  outfile.write('mut'+"\t"+"\t".join(localmaxs)+"\n")
  countmut = 0
  for mut in muts:
    countmut += 1
    if countmut%1000 == 0: print 'Finish searching path for %d variants' % countmut
    pathlengths = []
    for localmax in localmaxs:
      if nx.has_path(G,mut,localmax): pathlengths.append(int(nx.shortest_path_length(G, mut, localmax)))
      else: pathlengths.append(-1)
    outfile.write(mut+"\t"+"\t".join(map(str,pathlengths))+"\n")
  outfile.close()

def main():
  WT           = 'VDGV'
  mut          = 'WNWY'
  fitfile      = 'result/Mutfit'
  #missfitfile  = 'result/regression_missing'
  missfitfile  = 'result/regression_all_WT'
  localmaxfile = 'analysis/LocalMaxCompile_greedy'+'_pair'
  outfile      = 'analysis/LocalMaxPathLen'+'_pair'
  fcutoff      = -1
  condition    = 'I20fit'
  Index2pos    = {0:39,1:40,2:41,3:54}
  fithash      = TsvWithHeader2Hash(fitfile)
  print "Total # of variants in the raw data: %d" % len(fithash.keys())
  fithash      = filterfithash(fithash,condition,fcutoff)
  fithash      = fillinmissing(fithash,missfitfile,condition) 
  muts         = fithash.keys()
  print "Total # of variants after fill in with regression: %d" % len(muts)
  print "# of mutant pass cutoff: %d" % len(muts)
  localmaxs = []
  for line in open(localmaxfile,'r').xreadlines(): 
    if 'mut' in line: continue
    if float(line.rstrip().rsplit("\t")[2]) > 1: localmaxs.append(line.rstrip().rsplit("\t")[0])
  localmaxs.append('VDGV')
  print 'Total # of localmax = %s' % (len(localmaxs))
  G  = nx.DiGraph()
  G  = buildDigraph(G,muts,fithash,condition)
  print 'Finish building digraph with %d nodes and %d edges' % (len(G.nodes()), len(G.edges()))
  searchingpath(G,localmaxs,muts,outfile)

if __name__ == '__main__':
  main()
