#!/usr/bin/python
import os
import sys
import scipy 
import numpy as np
import scipy.cluster.hierarchy as hac

def hashing(infile,peakno):
  infile  = open(infile,'r')
  data = []
  IDs  = []
  for line in infile.xreadlines():
    if 'mut' in line: continue
    line = line.rstrip().rsplit("\t")
    mut  = line[0]
    des  = map(float,line[2:-1])
    des  = des[0:peakno]+[sum(des[peakno::])]
    data.append(des)
    IDs.append(mut)
  infile.close()
  data = np.array(data)
  return IDs, data
  
def clustering(IDs,data,t):
  print 'Start Clustering'
  z = hac.linkage(data)
  clusterinfo = ac.fcluster(z,20,criterion="maxclust")
  for n in range(0,len(IDs)):
    mut  = IDs[n]
    clus = clusterinfo[n]
    dat  = map(str,list(data[n]))
    print mut, clus, dat
    

def main():
  infile  = 'analysis/LocalMaxDes_weight'
  outfile = 'analysis/LocalMaxDesCluster_weight'
  peakno  = 15
  t       = 1
  IDs, data = hashing(infile,peakno)
  clustering(IDs,data,t)

if __name__ == '__main__':
  main()
