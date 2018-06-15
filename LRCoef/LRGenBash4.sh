#!/bin/bash

arg1[0]=65
arg1[1]=55
arg1[2]=45
arg1[3]=35

arg2[0]=40
arg2[1]=30
arg2[2]=20
arg2[3]=10

arg3[0]=40
arg3[1]=30
arg3[2]=20
arg3[3]=10

NMax=10

echo "${#arg1[@]}" > "outputData${#arg1[@]}".txt

for ((i=1;i<=$NMax;i++))
  do
    echo "$i/$NMax"
    for ((j=0;j<"${#arg1[@]}";j++))
      do
        arg1T[$j]=$((${arg1[$j]}*i))
        arg2T[$j]=$((${arg2[$j]}*i))
        arg3T[$j]=$((${arg3[$j]}*i))
      done
    num=$(lrcalc lrcoef ${arg1T[@]} - ${arg2T[@]} - ${arg3T[@]})
    echo "$i $num" >> "outputData${#arg1[@]}".txt
  done
