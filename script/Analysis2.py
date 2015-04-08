#!/usr/bim/python
import sys
import matplotlib.pyplot as plt
import networkx as nx
import operator
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

def labelnode(var,fit):
  scale = int(round((2-fit)*35+30))
  if fit >= 2: return 'grey30'
  else: return 'grey'+str(scale)

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

def findlocalmax(muts, fithash, condition):
  localmax = 0
  count    = 0
  process  = 0
  for mut in muts:
    process += 1
    if process%10000 == 0: print 'Finished %s localmax classificiation' % process
    mutfit   = float(fithash[mut][condition])
    variants = genvarsneighbor(mut)
    neighfit = []
    assert(len(variants)==76)
    for var in variants:
      if fithash.has_key(var): neighfit.append(float(fithash[var][condition]))
      else: neighfit.append('NA'); break
    if 'NA' not in neighfit:
      assert(len(neighfit)==76)
      count += 1
      if max(neighfit) < mutfit: localmax += 1
  return count, localmax

def drawgraph(G,outfile,fithash,WT,mut,condition):
  outfile=open(outfile,'w')
  outfile.write('strict digraph{'+"\n"+"\t"+'rankdir=LR'+"\n"+"\t"+'node [shape=box]'+"\n")
  for var in G.nodes():
    fit = float(fithash[var][condition])
    col = labelnode(var,fit)
    outfile.write("\t"+var+' [fillcolor='+col+', color=black, style="filled,rounded"];'+"\n")
  for E in G.edges():
    var1 = E[0]
    var2 = E[1]
    if hamming(WT,var2) != hamming(WT,var1):
      if hamming(WT,var1) < hamming(WT,var2):   s = var1; t = var2
      elif hamming(WT,var2) < hamming(WT,var1): s = var2; t = var1
      outfile.write("\t"+s+'->'+t+' [style=bold, color=black];'+"\n")
    elif hamming(mut,var2) != hamming(mut,var1):
      if hamming(mut,var1) < hamming(mut,var2):   s = var2; t = var1
      elif hamming(mut,var2) < hamming(mut,var1): s = var1; t = var2
      outfile.write("\t"+s+'->'+t+' [style=bold, color=black];'+"\n")
    elif hamming(mut,var2) == hamming(mut,var1) and hamming(WT,var2) == hamming(WT,var1): 
      outfile.write("\t"+var1+'->'+var2+' [style=bold, color=black];'+"\n")
      outfile.write("\t"+var2+'->'+var1+' [style=bold, color=black];'+"\n")
    else: print 'Error in graph construction'; sys.exit()
  outfile.write('}'+"\n")
  outfile.close()

def parsenodefile(nodefile, mut, WT, Index2pos):
  infile = open(nodefile,'r')
  paths  = {}
  countpath = 0
  for line in infile.xreadlines():
    countpath += 1
    steps = line.rstrip().rsplit(' ')
    ID    = [str(countpath)]
    for step in steps:
      for n in range(len(step)):
        if step[n] != WT[n] and step[n] != mut[n]:
          bridgemut = WT[n]+str(Index2pos[n])+step[n]
          if bridgemut not in ID: ID.append(bridgemut)
    paths['-'.join(ID)] = steps
  infile.close()
  print 'Read in a total of %d pathways' % countpath
  return paths

def analyzepaths(paths, fithash, condition):
  pathfit = {}
  for ID in paths.keys():
    minfitstep = 99
    steps = paths[ID]
    for step in steps:
      fit = float(fithash[step][condition])
      if fit < minfitstep: minfitstep = fit
    pathfit[ID] = minfitstep
  bridgemuts = []
  [bridgemuts.extend(ID.rsplit('-')[1::]) for ID in paths.keys()]
  bridgemuts = Counter(bridgemuts)
  print 'Fitness valley of each pathway'
  pathmaxfit = ''
  maxfit     = -1
  for ID in pathfit.keys():
    print ID+"\t"+str(pathfit[ID])
    if pathfit[ID] > maxfit: 
      maxfit     = pathfit[ID]
      pathmaxfit = ID+"\t"+'->'.join(paths[ID])
  print 'Common steps:'
  for bridgemut in bridgemuts.keys():
    print bridgemut+"\t"+str(bridgemuts[bridgemut])
  print "pathway with a valley with minimum fitness cost: %s" % pathmaxfit
  print 'fitness cost = %f' % maxfit
  return pathmaxfit.rsplit("\t")[1].rsplit('->')

def main():
  WT         = 'VDGV'
  mut        = 'WNWY'
  fitfile    = 'result/Mutfit'
  nodefile   = 'analysis/WNWYpathsCutoff10fold'
  outfile    = 'WNWYpaths.dot'
  fcutoff    = 0.1
  condition  = 'I20fit'
  Index2pos  = {0:39,1:40,2:41,3:54}
  fithash    = TsvWithHeader2Hash(fitfile)
  print "Total # of variants: %d" % len(fithash.keys())
  fithash    = filterfithash(fithash,condition,fcutoff)
  muts       = fithash.keys()
  print "# of mutant pass cutoff: %d" % len(muts)

  #Draw customized graphs
  paths      = parsenodefile(nodefile, mut, WT, Index2pos)
  usefulnodes= []
  [usefulnodes.extend(path) for path in paths.values()]
  usefulnodes= analyzepaths(paths, fithash, condition)
  G          = buildgraph(usefulnodes)
  paths = all_shortest_paths(G,WT,mut)
  print paths
  drawgraph(G,outfile,fithash,WT,mut,condition)
  sys.exit()

  #Find Local Max
#  Varwithallneighbors, localmax = findlocalmax(muts, fithash, condition)
#  print 'Total # of variants with all neighbors profiled: %s' % Varwithallneighbors,
#  print 'Among these, # of variants as a local max: %s' % localmax

  #Global Analysis
  G          = buildgraph(muts)
  paths = all_shortest_paths(G,WT,mut)
  for path in paths: print path
#  print '# of connected components: %d' % len(connected_components(G))
#  nx.draw_circular(G)
#  plt.savefig("TestGraph.png")

if __name__ == '__main__':
  main()
