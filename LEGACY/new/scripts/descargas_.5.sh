#!/bin/bash

comm="python dr2_result.py"
for n in {0..8}; do
  for m in {0..9}; do
    num=$[10*$n+$m]
      input=`ls ./sessions/session2_*dec_-$num.5.txt`
      output=~/cdn.gea.esac.esa.int/results_wb/south_hemisphere/dec-$n$m.5
      $comm $input $output
  done
done

