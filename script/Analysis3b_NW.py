#!/usr/bin/python
import os
import sys
import operator
import networkx as nx
import numpy as np
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

def filterfithash(fithash, condition):
  for mut in fithash.keys():
    if fithash[mut][condition] == 'NA' or '_' in mut:
      del fithash[mut]
  return fithash

def generatenodes(mut,WT,fithash):
  pathhash = {}
  nodes    = []
  for i in range(4): pathhash[i] = [WT[i],mut[i]]
  for i in range(16): 
    index = str(bin(i)).rsplit('b')[1]
    index = ''.join(map(str,[0]*(4-len(index))))+index
    node  = ''
    for p in range(len(index)):
      node+=pathhash[p][int(index[p])]
    if fithash.has_key(node): nodes.append(node)
    else: nodes.append('NA')
  return nodes

def buildgraph(nodes):
  G = nx.Graph()
  [G.add_node(node) for node in nodes] 
  for n1 in nodes:
    for n2 in nodes: 
      if hamming(n1,n2) == 1: G.add_edge(n1,n2)
  return G

def removenodes(nodes, fithash, cutoff, condition):
  cleannodes = []
  for node in nodes:
    if float(fithash[node][condition]) >= cutoff:
      cleannodes.append(node)
  return cleannodes

def stucking(path,fithash,condition,mut,WT):
  mutfit = float(fithash[mut][condition])
  wtfit  = float(fithash[WT][condition])
  for step in path:
    stepfit = float(fithash[step][condition])
    if step != WT and step != mut and stepfit > mutfit and stepfit > wtfit: return 1
  return 0

def monoincr(path,fithash,condition,mut,WT):
  for n in range(1,len(path)):
    stepPrev = path[n-1]
    stepCurr = path[n]
    if float(fithash[stepPrev][condition]) > float(fithash[stepCurr][condition]): return 0
  return 1

def localmaxing(G,WT,fithash,condition):
  localmaxs = []
  for node in G.nodes():
    fit = float(fithash[node][condition])
    neighfits = []
    [neighfits.append(float(fithash[neigh][condition])) for neigh in G.neighbors(node)]
    if max(neighfits) < fit: localmaxs.append(node)
  return localmaxs

def pathwayanalysis(fithash,muts,WT,condition,outfile):
  print 'Pathway analysis started'
  for mut in muts:
    HD = int(fithash[mut]['HD'])
    mutfit = float(fithash[mut][condition])
    if HD == 4 and mutfit < 1 and mutfit > 0.01:
      nodes = generatenodes(WT,mut,fithash)
      if 'NA' in nodes: continue
      G     = buildgraph(nodes)
      localmaxs      = localmaxing(G,WT,fithash,condition)
      if len(localmaxs) != 1 or localmaxs[0] != WT: continue
      assert(max([float(fithash[node][condition]) for node in nodes])==1)
      paths = all_shortest_paths(G,mut,WT)
      monoIpath = 0
      for path in paths: 
        assert(len(path)==5)
        #if float(fithash[path[1]][condition]) <= 0.01: continue #Filter first step fitness
        monoIpath+=monoincr(path,fithash,condition,mut,WT)
      outfile.write("\t".join(map(str,[condition, mut, mutfit, len(nodes), monoIpath]))+"\n")

def main():
  WT         = 'VDGV'
  fitfile    = 'result/Mutfit'
  outfile    = 'analysis/ShortMonoPaths4ToWT'
  condition  = 'I20fit'
  fithash    = TsvWithHeader2Hash(fitfile)
  fithash    = filterfithash(fithash, condition) 
  muts       = fithash.keys()
  outfile    = open(outfile,'w')
  outfile.write("Condition\tMut\tMutfit\tNodes\tMonoIPath\n")
  pathwayanalysis(fithash,muts,WT,condition,outfile)
  outfile.close()

if __name__ == '__main__':
  main()
