#!/bin/bash

comm="python dr2_result.py"
for n in {0..8}; do
  for m in {0..9}; do
    num=$[10*$n+$m]
    if [ $num != 0 ]
    then
      input=`ls ./sessions/session152*dec_-$num.txt`
      output=~/cdn.gea.esac.esa.int/results_wb/south_hemisphere/dec-$n$m
      $comm $input $output
    fi
  done
done

