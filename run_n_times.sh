#usr/bin/bash
echo "---start execution---" > out.txt;
n=0;
while [[ $n -lt $1 ]]; 
    do
    echo "iteration : $n" >> out.txt;
    ./out/openmp 16 4 >> out.txt; n=$((n+1));
    done
