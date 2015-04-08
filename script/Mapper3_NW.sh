#!/bin/bash
for f in pep/*
do
outfile=${f/\.pep/}
outfile=${outfile/pep\//}
echo $outfile
sort $f | uniq -c | sort -nr > count/$outfile\.count
done
