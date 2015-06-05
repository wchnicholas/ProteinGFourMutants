#!/usr/bim/python
import os
import sys
import matplotlib.pyplot as plt
import networkx as nx
import operator
import colorsys
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

def labelnode(var,fit,highcolscale):
  high = float(highcolscale)
  low  = float(1)
  mid  = float(high-low)/2+low
  
  #USE BLUE MONOTONE GRADIENT
  if fit > high: return colorsys.rgb_to_hsv(0, 0, 1)
  elif fit > low: return colorsys.rgb_to_hsv(1-(fit-low)/(high-low),1-(fit-low)/(high-low),1)
  else: return colorsys.rgb_to_hsv(1,1,1)
  print 'Leaking'; sys.exit()

  if fit > high:  return colorsys.rgb_to_hsv(1, 0, 0)
  elif fit > mid: return colorsys.rgb_to_hsv(1,(high-fit)/(high-mid),0)
  elif fit > low: return colorsys.rgb_to_hsv(1,1,(mid-fit)/(mid-low))
  elif fit < low: return 'black'
  else: print 'Something is wrong with the coloring function'; sys.exit()

def drawgraph(nodes,outfile,scalefile,fithash,condition,highcolscale):
  normscale = 50000
  scaling   = 4
  edgecount = 0
  scalefile=open(scalefile,'w')
  scalefile.write('strict graph{'+
                  "\n"+"\t"+'graph [ dpi = 30 ]'+
                  "\n"+"\t"+'rankdir=LR'+
                  "\n"+"\t"+'node [shape=circle,label=""]'+"\n")
  for basin in [500,5000,50000]:
    size  = str(float(basin)/float(normscale)*scaling)
    basin = str(basin)
    scalefile.write("\t"+basin+' [fillcolor=white, color=black, style="filled,rounded", fixedsize=shape, width='+size+'];'+"\n")
  scalefile.write('}'+"\n")
  scalefile.close()
  outfile=open(outfile,'w')
  #outfile.write('strict graph{'+
  #              "\n"+"\t"+'graph [ dpi = 30 ]'+
  #              "\n"+"\t"+'layout=fdp'+
  #              "\n"+"\t"+'node [shape=circle,label=""]'+"\n")
  outfile.write('strict graph{'+"\n"+"\t"+'layout=fdp'+"\n"+"\t"+'node [shape=circle]'+"\n")
  for var in nodes:
    if fithash.has_key(var): fit = float(fithash[var][condition])
    else: fit = 'NA'
    col  = labelnode(var,fit,highcolscale)
    col  = ','.join(map(str,list(col)))
    size = str(float(fithash[var]['basin'])/float(normscale)*scaling)
    outfile.write("\t"+var+' [fillcolor="'+col+'", color=black, style="filled,rounded", fixedsize=shape, width='+size+'];'+"\n")
  for i in range(len(nodes)):
    for j in range(len(nodes)):
      if i > j:
        var1 = nodes[i]
        var2 = nodes[j]
        dist = hamming(var1,var2)
        if dist == 2: outfile.write("\t"+var1+'--'+var2+' [style=bold, color=black];'+"\n"); edgecount += 1
  outfile.write('}'+"\n")
  outfile.close()
  print 'Total number of edge: %d' % edgecount

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

def main(graphfile,graphscale,localmaxfile,highcolscale):
  print 'working on localmaxfile: %s' % localmaxfile
  WT           = 'VDGV'
  condition    = 'I20fit'
  fcutoff      = float(1)
  Index2pos    = {0:39,1:40,2:41,3:54}
  localmaxhash = TsvWithHeader2Hash(localmaxfile)
  print "Total # of variants as local max: %d" % len(localmaxhash.keys())
  localmaxhash = filterfithash(localmaxhash,condition,fcutoff)
  localmaxs    = localmaxhash.keys()
  print "# of localmax pass cutoff: %d" % len(localmaxs)

  #Draw graphs of local maximum
  drawgraph(localmaxs,graphfile,graphscale,localmaxhash,condition,highcolscale)
  os.system('dot -Tpng %s -o %s' % (graphfile,graphfile.replace('.dot','.png')))
  os.system('dot -Tpng %s -o %s' % (graphscale,graphscale.replace('.dot','.png')))
  #os.system('rm %s' % graphfile)
  os.system('rm %s' % graphscale)
  print 'DONE\n'

if __name__ == '__main__':
  #main('xdot/LocalMax_pair.dot','xdot/LocalScale_pair.dot','analysis/LocalMaxCompile_weight_pair1000sim',30)
  main('xdot/LocalMax.dot','xdot/LocalScale.dot','analysis/LocalMaxCompile_weight1000sim',8)
  #main('xdot/LocalMax_greedy.dot','xdot/LocalScale_greedy.dot','analysis/LocalMaxCompile_greedy',8)
  #main('xdot/LocalMax_random.dot','xdot/LocalScale_random.dot','analysis/LocalMaxCompile_random1000sim',8)
