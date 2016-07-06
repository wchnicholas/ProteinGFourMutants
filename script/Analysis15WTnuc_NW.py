#!/usr/bin/python
import sys
import operator
import random
import networkx as nx
import matplotlib.pyplot as plt
from math import exp
from itertools import imap

def codondistancemap():
  codondists = {}
  dnamap = {"TTT":"F", "TTC":"F", "TTA":"L", "TTG":"L",
    "TCT":"S", "TCC":"S", "TCA":"S", "TCG":"S",
    "TAT":"Y", "TAC":"Y", "TAA":"_", "TAG":"_",
    "TGT":"C", "TGC":"C", "TGA":"_", "TGG":"W",
    "CTT":"L", "CTC":"L", "CTA":"L", "CTG":"L",
    "CCT":"P", "CCC":"P", "CCA":"P", "CCG":"P",
    "CAT":"H", "CAC":"H", "CAA":"Q", "CAG":"Q",
    "CGT":"R", "CGC":"R", "CGA":"R", "CGG":"R",
    "ATT":"I", "ATC":"I", "ATA":"I", "ATG":"M",
    "ACT":"T", "ACC":"T", "ACA":"T", "ACG":"T",
    "AAT":"N", "AAC":"N", "AAA":"K", "AAG":"K",
    "AGT":"S", "AGC":"S", "AGA":"R", "AGG":"R",
    "GTT":"V", "GTC":"V", "GTA":"V", "GTG":"V",
    "GCT":"A", "GCC":"A", "GCA":"A", "GCG":"A",
    "GAT":"D", "GAC":"D", "GAA":"E", "GAG":"E",
    "GGT":"G", "GGC":"G", "GGA":"G", "GGG":"G",}
  aas = list(set(dnamap.values()))
  for aa1 in aas:
    for aa2 in aas:
      codons1 = [codon for codon in dnamap.keys() if dnamap[codon] == aa1]
      codons2 = [codon for codon in dnamap.keys() if dnamap[codon] == aa2]
      dist    = min([hamming(codon1, codon2) for codon1 in codons1 for codon2 in codons2])
      codondists[aa1+'-'+aa2] = dist
      codondists[aa2+'-'+aa1] = dist
  return codondists

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

def shortest_path(G, source=None, target=None, weight=None):
    if source is None:
        if target is None:
            ## Find paths between all pairs.
            if weight is None:
                paths=nx.all_pairs_shortest_path(G)
            else:
                paths=nx.all_pairs_dijkstra_path(G,weight=weight)
        else:
            ## Find paths from all nodes co-accessible to the target.
            with nx.utils.reversed(G):
                if weight is None:
                    paths=nx.single_source_shortest_path(G, target)
                else:
                    paths=nx.single_source_dijkstra_path(G, target, weight=weight)

                # Now flip the paths so they go from a source to the target.
                for target in paths:
                    paths[target] = list(reversed(paths[target]))

    else:
        if target is None:
            ## Find paths to all nodes accessible from the source.
            if weight is None:
                paths=nx.single_source_shortest_path(G,source)
            else:
                paths=nx.single_source_dijkstra_path(G,source,weight=weight)
        else:
            ## Find shortest source-target path.
            if weight is None:
                paths=nx.bidirectional_shortest_path(G,source,target)
            else:
                paths=nx.dijkstra_path(G,source,target,weight)

    return paths

 
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

def ShortestAAPath(pep1, pep2, AAgraph):
  return sum([len(shortest_path(AAgraph, aa1, aa2))-1 for aa1, aa2 in zip(pep1, pep2)])

def MinNucDist(pep1, pep2, codondists):
  return sum([codondists[aa1+'-'+aa2] for aa1, aa2 in zip(pep1, pep2)])

def buildaagraph(codondists):
  aas = ['E','D','R','K','H','Q','N','S','T','P','G','C','A','V','I','L','M','F','Y','W']
  AAgraph = nx.Graph()
  [AAgraph.add_node(aa) for aa in aas]
  for aa in aas: 
    variants = aas
    nucvariants = [var for var in variants if MinNucDist(var, aa, codondists) == 1]
    [AAgraph.add_edge(aa,var) for var in nucvariants]
  return AAgraph

def buildDigraph(G,muts,fithash,condition,codondists):
  for mut in muts:
    mutfit = float(fithash[mut][condition])
    G.add_node(mut)
    variants = genvarsneighbor(mut)
    nucvariants = [var for var in variants if MinNucDist(var, mut, codondists) == 1]
    [G.add_edge(mut,var) for var in nucvariants if float(fithash[var][condition]) > mutfit]
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

def AAtransitionmatrix(AAgraph, outfile):
  outfile = open(outfile,'w')
  outfile.write("\t".join(['aa1', 'aa2', 'value'])+"\n")
  for aa1 in AAgraph.nodes():
    for aa2 in AAgraph.nodes():
      if aa1 == aa2:
        outfile.write("\t".join([aa1, aa2, str(-1)])+"\n")
      elif AAgraph.has_edge(aa1,aa2):
        outfile.write("\t".join([aa1, aa2, str(1)])+"\n")
      else: 
        outfile.write("\t".join([aa1, aa2, str(0)])+"\n")
  outfile.close()
 
def main():
  WT           = 'VDGV'
  fitfile      = 'result/Mutfit'
  missfitfile  = 'result/regression_missing'
  #missfitfile  = 'result/regression_all_WT'
  outfile      = 'analysis/LocalMaxEvolvePotWTnuc'#+'_pair'
  degfile      = 'analysis/LocalMaxEvolveDegreeWTnuc'#+'_pair'
  fcutoff      = -1
  fitfold      = float(1)
  condition    = 'I20fit'
  Index2pos    = {0:39,1:40,2:41,3:54}
  codondists   = codondistancemap() 
  AAgraph      = buildaagraph(codondists)
  AAtransitionmatrix(AAgraph, 'analysis/AAtransitionmatrix')
  fithash      = TsvWithHeader2Hash(fitfile,condition)
  print "Total # of variants in the raw data: %d" % len(fithash.keys())
  fithash      = fillinmissing(fithash,missfitfile,condition) 
  muts         = fithash.keys()
  print "Total # of variants after fill in with regression: %d" % len(muts)
  print "# of mutant pass cutoff: %d" % len(muts)
  G  = nx.DiGraph()
  G  = buildDigraph(G,muts,fithash,condition,codondists)
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
    EndHD  = hamming(WT, End)
    AADist = ShortestAAPath(WT, End, AAgraph)
    if ReachEnds.has_key(End): Endpl = len(ReachEnds[End])-1; Endpath = '->'.join(ReachEnds[End])
    else: Endpl = -1; Endpath = 'NA'
    if Endpl == AADist: Direct = 'Yes'
    elif Endpl == -1: Direct = 'Inaccessible'
    elif AADist < Endpl: Direct = 'No'
    else: print "Something is wrong"; sys.exit()
    outfile.write("\t".join(map(str,[End,EndHD,Endfit,Endpl,Endpath,Direct]))+"\n")
  outfile.close()

if __name__ == '__main__':
  main()
