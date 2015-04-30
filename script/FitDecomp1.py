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

def pathwayanalysis(fithash,muts,WT,condition):
  pathinfohash = {}
  for mut in muts:
    HD = int(fithash[mut]['HD'])
    if HD == 4:
      nodes = generatenodes(mut,WT,fithash)
      if 'NA' in nodes: continue
      if WT in nodes: nodes.remove(WT)
      fits = [float(fithash[node][condition]) for node in nodes]
      pathinfohash[mut] = {}
      pathinfohash[mut]['minfit']  = min(fits)
      pathinfohash[mut]['maxfit']  = max(fits)
      pathinfohash[mut]['ragfit']  = max(fits)-min(fits)
  return pathinfohash
       
def compileout(infile, outfile, pathinfohash):
  infile  = open(infile,'r')
  outfile = open(outfile,'w')
  for line in infile.xreadlines():
    line = line.rstrip()
    if 'genotype_all' in line: outfile.write(line+"\t"+"\t".join(['FitMax','FitMin','FitRange'])+"\n")
    else: 
      mut    = line.rsplit("\t")[0]
      fitmax = pathinfohash[mut]['maxfit']
      fitmin = pathinfohash[mut]['minfit']
      fitrag = pathinfohash[mut]['ragfit']
      outfile.write(line+"\t"+"\t".join(map(str,[fitmax,fitmin,fitrag]))+"\n")
  infile.close()
  outfile.close()

def main():
  WT         = 'VDGV'
  fitfile    = 'result/Mutfit'
  infile     = 'analysis/FitnessDecompose'
  outfile    = 'analysis/FitnessDecomposeFit'
  condition  = 'I20fit'
  fithash    = TsvWithHeader2Hash(fitfile)
  fithash    = filterfithash(fithash, condition)  
  muts       = fithash.keys()
  header     = "\t".join(['Conditions','Cutoff'])
  pathinfohash = pathwayanalysis(fithash,muts,WT,condition)
  compileout(infile, outfile, pathinfohash)

if __name__ == '__main__':
  main()
