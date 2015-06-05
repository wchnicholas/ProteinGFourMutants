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

def removenodes(nodes, fithash, cutoff, condition, mut, WT):
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

def pathwayanalysis(fithash,muts,WT,condition,cutoff,detailcutoffs):
  print 'Pathway analysis started. Cutoff = %f Condition = %s' % (cutoff, condition)
  if cutoff in detailcutoffs:  detailfile = open('analysis/ShortPaths4mutsCutoff'+str(int(1/cutoff))+'fold','w')
  if cutoff in detailcutoffs:  detailfile.write("Condition\tCutoff\tMut\tMutfit\tNodes\tPaths\tStuckPath\tMonoIPath\n")
  num_mut  = 0
  num_node = []
  num_path = []
  monoIpaths = []
  stuckpaths = []
  for mut in muts:
    HD = int(fithash[mut]['HD'])
    mutfit = float(fithash[mut][condition])
    if mutfit < 1: continue
    if HD == 4 and mutfit >= cutoff:
      nodes = generatenodes(mut,WT,fithash)
      if 'NA' in nodes: continue
      nodes = removenodes(nodes, fithash, cutoff, condition, mut, WT)
      G     = buildgraph(nodes)
      if nx.has_path(G,WT,mut): 
        paths = all_shortest_paths(G,WT,mut)
        stuckpath = 0
        countpath = 0
        monoIpath = 0
        for path in paths: 
          assert(len(path)==5)
          countpath+=1
          stuckpath+=stucking(path,fithash,condition,mut,WT)
          monoIpath+=monoincr(path,fithash,condition,mut,WT)
      else: countpath = 0; stuckpath = 0; monoIpath = 0; stuckpath = 0
      num_mut += 1
      num_node.append(len(nodes))
      num_path.append(countpath)
      monoIpaths.append(monoIpath)
      stuckpaths.append(stuckpath)
      if cutoff in detailcutoffs: detailfile.write("\t".join(map(str,[condition, cutoff, mut, mutfit, len(nodes), countpath, stuckpath, monoIpath]))+"\n")
  if cutoff in detailcutoffs: detailfile.close()
  return num_mut, num_node, num_path, monoIpaths, stuckpaths
   

def main():
  WT         = 'VDGV'
  fitfile    = 'result/Mutfit'
  outfile    = 'analysis/ShortPaths4muts'
  #conditions = ['I10fit','I20fit','I90fit']
  conditions = ['I20fit']
  cutoffs    = [0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1]
  detailcutoffs  = [0.1,1]
  fithash    = TsvWithHeader2Hash(fitfile)
  fithash    = filterfithash(fithash, conditions[0])  
  muts       = fithash.keys()
  header     = "\t".join(['Conditions','Cutoff','NumMut','MeanNumNode','StdNumNode','MeanNumPath','StdNumPath','MeanMonoIPath','StdMonoIPath','MeanStuckPath','StdStuckPath'])
  outfile    = open(outfile,'w')
  outfile.write(header+"\n")
  for condition in conditions:
    for cutoff in cutoffs:
      num_mut, num_node, num_path, monoIpaths, stuckpaths = pathwayanalysis(fithash,muts,WT,condition,cutoff,detailcutoffs)
      outfile.write("\t".join(map(str,[condition, cutoff, num_mut, np.mean(num_node), np.std(num_node), np.mean(num_path), np.std(num_path),np.mean(monoIpaths),np.std(monoIpaths),np.mean(stuckpaths),np.std(stuckpaths)]))+"\n")
  outfile.close()

if __name__ == '__main__':
  main()
