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
    if 'genotype_' in line: continue
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

def degreeanalysis(degfile,G,fithash,condition):
  outfile = open(degfile,'w')
  header  = "\t".join(map(str,['mut','fit','indegree','outdegree']))
  outfile.write(header+"\n")
  for node in G.nodes():
    indeg  = G.in_degree(node)
    outdeg = G.out_degree(node)
    varfit = fithash[node][condition]
    outfile.write("\t".join(map(str,[node,varfit,indeg,outdeg]))+"\n")
  outfile.close()

def main():
  WT           = 'VDGV'
  fitfile      = 'result/Mutfit'
  missfitfile  = 'result/regression_missing'
  #missfitfile  = 'result/regression_all_WT'
  outfile      = 'analysis/LocalMaxEvolvePotWT'#+'_pair'
  degfile      = 'analysis/LocalMaxEvolveDegreeWT'#+'_pair'
  fcutoff      = -1
  fitfold      = float(1)
  condition    = 'I20fit'
  Index2pos    = {0:39,1:40,2:41,3:54}
  fithash      = TsvWithHeader2Hash(fitfile,condition)
  print "Total # of variants in the raw data: %d" % len(fithash.keys())
  fithash      = fillinmissing(fithash,missfitfile,condition) 
  muts         = fithash.keys()
  print "Total # of variants after fill in with regression: %d" % len(muts)
  print "# of mutant pass cutoff: %d" % len(muts)
  G  = nx.DiGraph()
  G  = buildDigraph(G,muts,fithash,condition)
  print 'Finish building digraph with %d nodes and %d edges' % (len(G.nodes()), len(G.edges()))
  print 'Analyzing degrees for each node'
  degreeanalysis(degfile,G,fithash,condition)
  outfile = open(outfile,'w')
  header  = "\t".join(['mut','HD','fit','PathLength','Path','Direct'])
  outfile.write(header+"\n")
  WTfit   = float(fithash[WT][condition]) 
  print 'Working %s with fitness %f' % (WT, WTfit)
  AllEnds   = [j for j in muts if float(fithash[j][condition]) > WTfit*fitfold]
  ReachEnds = single_source_shortest_path(G,WT)
  del ReachEnds[WT]
  for End in ReachEnds.keys():
    if float(fithash[End][condition]) <= WTfit*fitfold: del ReachEnds[End]
  print "Total Number of variants with fitness higher than %f = %d (Reachable = %d)" % (WTfit,len(AllEnds),len(ReachEnds.keys()))
  for End in AllEnds:
    if End == 'VDGV': continue
    Endfit = fithash[End][condition]
    EndHD  = hamming(WT,End)
    if ReachEnds.has_key(End): Endpl = len(ReachEnds[End])-1; Endpath = '->'.join(ReachEnds[End])
    else: Endpl = -1; Endpath = 'NA'
    if Endpl == EndHD: Direct = 'Yes'
    elif Endpl == -1: Direct = 'Inaccessible'
    elif EndHD < Endpl: Direct = 'No'
    else: print "Something is wrong"; sys.exit()
    outfile.write("\t".join(map(str,[End,EndHD,Endfit,Endpl,Endpath,Direct]))+"\n")
  outfile.close()

if __name__ == '__main__':
  main()
