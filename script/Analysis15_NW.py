#!/usr/bin/python
import sys
import operator
import random
import networkx as nx
from math import exp
from itertools import imap

def single_source_shortest_path(G,source,cutoff=None):
    level=0                  # the current level
    nextlevel={source:1}       # list of nodes to check at next level
    paths={source:[source]}  # paths dictionary  (paths to key from source)
    if cutoff==0:
        return paths
    while nextlevel:
        thislevel=nextlevel
        nextlevel={}
        for v in thislevel:
            for w in G[v]:
                if w not in paths:
                    paths[w]=paths[v]+[w]
                    nextlevel[w]=1
        level=level+1
        if (cutoff is not None and cutoff <= level):  break
    return paths
 
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

def TsvWithHeader2Hash(fitfile,condition):
  H = {}
  infile = open(fitfile,'r')
  countline = 0
  header = []
  for line in infile.xreadlines():
    countline += 1
    line = line.rstrip().rsplit("\t")
    if countline == 1: header = line; continue
    mut = line[0]
    if '_' in str(mut): continue
    for i in range(1,len(line)):
      if header[i] == condition and str(line[i]) != 'NA':
        H[mut] = {}
        H[mut][header[i]] = line[i] 
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
    if 'genotype_missing' in line: continue
    line = line.rstrip().rsplit("\t")
    mut  = line[0]
    fit  = exp(float(line[1]))
    #if mut in fithash.keys(): print 'Variant %d is not a missing data' % mut; sys.exit()
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

def mutsampling(fithash,condition,muts,samplingrange,samplesize):
  i = 1
  samplemuts = []
  while i < len(samplingrange):
    n = 0
    for mut in muts:
      mutfit = float(fithash[mut][condition])
      if mutfit < float(samplingrange[i]) and mutfit > float(samplingrange[i-1]): 
        samplemuts.append(mut)
        n += 1
      if samplesize == n: break
    i += 1
  return samplemuts

def main():
  fitfile      = 'result/Mutfit'
  missfitfile  = 'result/regression_missing'
  outfile      = 'analysis/LocalMaxEvolvePotSampled'
  fcutoff      = -1
  condition    = 'I20fit'
  Index2pos    = {0:39,1:40,2:41,3:54}
  fithash      = TsvWithHeader2Hash(fitfile,condition)
  print "Total # of variants in the raw data: %d" % len(fithash.keys())
  fithash      = fillinmissing(fithash,missfitfile,condition) 
  muts         = fithash.keys()
  random.shuffle(muts)
  samplingrange = [0, 0.2,0.4,0.6,0.8,1.0,1.2,1.4,1.6,1.8,2.0,
                      2.2,2.4,2.6,2.8,3.0,3.2,3.4,3.6,3.8,4.0,
                      4.2,4.4,4.6,4.8,5.0,5.2,5.4,5.6,5.8,6.0]
  samplesize = 5
  samplemuts   = mutsampling(fithash,condition,muts,samplingrange,samplesize)
  print "Total number of sampled variants for path analysis: %d" % len(samplemuts)
  print "Total # of variants after fill in with regression: %d" % len(muts)
  print "# of mutant pass cutoff: %d" % len(muts)
  G  = nx.DiGraph()
  G  = buildDigraph(G,muts,fithash,condition)
  print 'Finish building digraph with %d nodes and %d edges' % (len(G.nodes()), len(G.edges()))
  outfile = open(outfile,'w')
  header  = "\t".join(['mut','fit','AllEnds','ReachEnds'])
  outfile.write(header+"\n")
  countmut = 0
  for i in samplemuts: 
    countmut += 1
    print 'Working on %d variant with fitness %f .... %s' % (countmut,float(fithash[i][condition]),i)
    mutfit  = float(fithash[i][condition]) 
    AllEnds = len([j for j in muts if float(fithash[j][condition]) > mutfit])
    ReachEnds = len(single_source_shortest_path(G,i).keys())-1
    info      = [AllEnds, ReachEnds]
    outfile.write(i+"\t"+str(mutfit)+"\t"+"\t".join(map(str,info))+"\n")
  outfile.close()

if __name__ == '__main__':
  main()
