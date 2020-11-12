#!/bin/bash

#
# clean up (i.e. remove) all files produced by a rs@md run
#

if [[ $# -lt 1 ]]; then
    echo "usage: $0 <#cycles to clean up>"
    exit 0
fi 

[[ $# -eq 1 ]] && nCycles=$1
echo "will remove files for $nCycles cycles ..."

rm 0.top 0-*

for cycle in `seq 1 $nCycles`; do
    [ -f ${cycle}.top ] && rm ${cycle}.top 
    [ -f ${cycle}-rs.tpr ] && rm ${cycle}-rs* 
    [ -f ${cycle}-md.tpr ] && rm ${cycle}-md* 
    [ -f rejected-${cycle}-rs.tpr ] && rm rejected-${cycle}*
    [ -f ${cycle}.products.ndx ] && rm ${cycle}.products.ndx ${cycle}.reactants.ndx
done

[ -f reactants.tpr ] && rm reactants.* products.*
[ -f reactants_solvation.tpr ] && rm reactants_solvation* products_solvation*
rm statistics.data

echo "... done"

