#!/usr/bin/python
import os
import sys
import operator
import networkx as nx
import numpy as np
from itertools import imap

def floor(fit):
  if fit == 'NA': return 'NA'
  elif float(fit) < 0.01: return 0.01
  else: return float(fit)

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

def stephashing(stepprobhash,step,mut,WT,fithash,condition,fixtype):
  n1 =  mut[0]+step[1]+step[2]+step[3]
  n2 = step[0]+ mut[1]+step[2]+step[3]
  n3 = step[0]+step[1]+ mut[2]+step[3]
  n4 = step[0]+step[1]+step[2]+ mut[3]
  nlist = list(set([n1,n2,n3,n4]))
  step_fit = floor(float(fithash[step][condition]))
  nsels    = {}
  while step in nlist: nlist.remove(step)
  for n in nlist: 
    n_sel = floor(float(fithash[n][condition]))-step_fit
    if n_sel < 0: nsels[n] = 0
    else: 
      if fixtype == 'equal': nsels[n] = 1
      if fixtype == 'prop':  nsels[n] = n_sel
  if sum(nsels.values()) != 0: 
    sumsel = sum(nsels.values())
    for n in nlist: stepprobhash[step+'->'+n] = float(nsels[n])/float(sumsel)
  else: 
    for n in nlist: stepprobhash[step+'->'+n] = float(0)
  return stepprobhash

def calprobpath(path,stepprobhash):
  P = float(1)
  for n in range(1,len(path)):
    step = path[n-1]+'->'+path[n]
    P = P*stepprobhash[step]
  return P

def localmaxing(G,WT,dist,fithash,condition):
  localmaxs = []
  for node in G.nodes(): 
    fit = float(fithash[node][condition])
    if hamming(WT,node) >= dist and fit > 1: 
      neighfits = [] 
      [neighfits.append(float(fithash[neigh][condition])) for neigh in G.neighbors(node)]
      if max(neighfits) < fit: localmaxs.append(node)
  return localmaxs

def monoincr(path,fithash,condition,mut,WT):
  for n in range(1,len(path)):
    stepPrev = path[n-1]
    stepCurr = path[n]
    if float(fithash[stepPrev][condition]) > float(fithash[stepCurr][condition]): return 0
  return 1

def monolocalmaxing(localmax,WT,fithash,condition,localmonopaths,G):
  paths = all_shortest_paths(G, WT, localmax)
  localmonopaths[localmax] = []
  [localmonopaths[localmax].append(path) for path in paths if monoincr(path,fithash,condition,localmax,WT) == 1]
  return localmonopaths

def giniindexing(pathprobs):
  giniindex = 0
  pathprobs = [float(0)]*(24-len(pathprobs))+sorted(pathprobs)
  accprobs  = [sum(pathprobs[0:i+1]) for i in range(0,len(pathprobs))]
  for n in range(len(accprobs)):
    if n == 0: giniindex += float(accprobs[n])/2
    else: giniindex += float(accprobs[n]+accprobs[n-1])/2
  maxarea   = float(len(accprobs))/2
  giniindex = (maxarea - float(giniindex))/maxarea
  return giniindex, accprobs
  
def pathwayparam(fithash,muts,WT,condition,outfile,HDdist,fixtype):
  countmut = 0
  countmutwithmonopath = 0
  for mut in muts:
    HD = int(fithash[mut]['HD'])
    mutfit = float(fithash[mut][condition])
    if HD == HDdist:
      countmut += 1
      if countmut%10000 == 0: print 'Processed %d HD = %d mutants' % (countmut, HDdist)
      nodes = generatenodes(mut,WT,fithash)
      if 'NA' in nodes: continue
      G     = buildgraph(nodes)
      localmaxs      = localmaxing(G,WT,1,fithash,condition)
      if len(localmaxs) == 0: continue
      stepprobhash   = {}
      localmonopaths = {}
      localmonomaxs  = []
      pathprobhash   = {}
      nodeprobhash   = {}
      nodeQ = 0
      pathQ = 0
      for localmax in localmaxs: localmonopaths = monolocalmaxing(localmax,WT,fithash,condition,localmonopaths,G)
      [localmonomaxs.append(localmonomax) for localmonomax in localmonopaths.keys() if len(localmonopaths[localmonomax]) != 0]
      if len(localmonomaxs) != 1 or hamming(localmonomaxs[0],WT) != 4: continue
      for step in G.nodes(): stepprobhash = stephashing(stepprobhash,step,mut,WT,fithash,condition,fixtype)
      for localmonomax in localmonomaxs:
        monopaths = localmonopaths[localmonomax]
        nodeprobhash[localmonomax] = 0
        for path in monopaths: 
          pathprob = calprobpath(path,stepprobhash)
          pathprobhash['->'.join(path)] = pathprob
          nodeprobhash[localmonomax] += pathprob
      pathprobnormalscale = sum(pathprobhash.values())
      nodeprobnormalscale = sum(nodeprobhash.values())
      for path in pathprobhash.keys(): pathprobhash[path] = float(pathprobhash[path])/float(pathprobnormalscale)
      for path in pathprobhash.keys(): pathQ+=float(pathprobhash[path])**2 
      for node in nodeprobhash.keys(): nodeprobhash[node] = float(nodeprobhash[node])/float(nodeprobnormalscale)
      for node in nodeprobhash.keys(): nodeQ+=float(nodeprobhash[node])**2
      assert(len(localmonomaxs) == 1 and hamming(localmonomaxs[0],WT) == 4)
      ginimax, maxaccprobs = giniindexing([1])
      giniindex, accprobs  = giniindexing(pathprobhash.values())
      giniindex = float(giniindex)/float(ginimax)
      accprobs  = "\t".join(map(str,accprobs))
      outfile.write("\t".join(map(str,[mut,HD,pathQ,nodeQ,len(localmonomaxs),len(pathprobhash.keys()),giniindex,accprobs]))+"\n")
      if len(localmonomaxs) > 0: countmutwithmonopath+=1
  print 'For HD = %d, there are %d beneficial variants with a monotonic paths' % (HDdist, countmutwithmonopath)
      #outfile.write(mut+"\t"+str(giniindex)+"\t"+"\t".join(map(str,accprob_paths))+"\n")
      
def main():
  WT         = 'VDGV'
  fitfile    = 'result/Mutfit'
  outfile    = 'analysis/PathwayParamResult'
  condition  = 'I20fit'
  fixtype    = 'equal' #Either 'prop' or 'equal'
  fithash    = TsvWithHeader2Hash(fitfile)
  fithash    = filterfithash(fithash, condition)  
  muts       = fithash.keys()
  header     = "\t".join([])
  outfile    = open(outfile+'.'+fixtype,'w')
  outfile.write("\t".join(['mut','HD','pathQ','nodeQ','localmax','monopaths','giniindex',"\t".join(map(str,range(1,25)))])+"\n")
  [pathwayparam(fithash,muts,WT,condition,outfile,HDdist,fixtype) for HDdist in [4]]
  outfile.close()

if __name__ == '__main__':
  main()
