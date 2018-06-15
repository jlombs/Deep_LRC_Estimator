#!/bin/bash

arg1[0]=3
arg1[1]=2
arg1[2]=1

arg2[0]=2
arg2[1]=1
arg2[2]=0

arg3[0]=2
arg3[1]=1
arg3[2]=0

NMax=100

echo "${#arg1[@]}" > "outputData${#arg1[@]}".txt

for ((i=1;i<=$NMax;i=i+10))
  do
    for ((j=0;j<"${#arg1[@]}";j++))
      do
        arg1T[$j]=$((${arg1[$j]}*i))
        arg2T[$j]=$((${arg2[$j]}*i))
        arg3T[$j]=$((${arg3[$j]}*i))
      done
    num=$(lrcalc lrcoef ${arg1T[@]} - ${arg2T[@]} - ${arg3T[@]})
    echo "$i $num" >> "outputData${#arg1[@]}".txt
  done
