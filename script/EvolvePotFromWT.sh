echo 'Analysis without codon contraint'
echo 'Hamming == 2'
awk {'if ($2==2) print'} analysis/LocalMaxEvolvePotWT | cut -f6 | sort | uniq -c
echo 'Hamming == 3'
awk {'if ($2==3) print'} analysis/LocalMaxEvolvePotWT | cut -f6 | sort | uniq -c
echo 'Hamming == 4'
awk {'if ($2==4) print'} analysis/LocalMaxEvolvePotWT | cut -f6 | sort | uniq -c

echo 'Analysis with codon constraint'
echo 'Hamming == 2'
awk {'if ($2==2) print'} analysis/LocalMaxEvolvePotWTnuc | cut -f6 | sort | uniq -c
echo 'Hamming == 3'
awk {'if ($2==3) print'} analysis/LocalMaxEvolvePotWTnuc | cut -f6 | sort | uniq -c
echo 'Hamming == 4'
awk {'if ($2==4) print'} analysis/LocalMaxEvolvePotWTnuc | cut -f6 | sort | uniq -c
