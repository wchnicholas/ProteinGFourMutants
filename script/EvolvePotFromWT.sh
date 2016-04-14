echo 'Hamming == 2'
awk {'if ($2==2) print'} analysis/LocalMaxEvolvePotWT | cut -f4 | sort | uniq -c
echo 'Hamming == 3'
awk {'if ($2==3) print'} analysis/LocalMaxEvolvePotWT | cut -f4 | sort | uniq -c
echo 'Hamming == 4'
awk {'if ($2==4) print'} analysis/LocalMaxEvolvePotWT | cut -f4 | sort | uniq -c
