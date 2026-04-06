#!/bin/sh
set -x
ev=ev2
statlat=46.768978
statlon=82.300659
OMPcorrectP1=0.0
OMPcorrectS1=0.0
TravTable=ak135_P1_S1_zerodepth.txt
BaziFile=${ev}_files/${ev}_backazi_distr.txt
P1ArrFile=${ev}_files/${ev}_P1_arrtime_distr.txt
S1ArrFile=${ev}_files/${ev}_S1_arrtime_distr.txt
N=10000
outfile=${ev}_locations.txt
stats_outfile=${ev}_location_stats.txt
python  single_array_PS_locate.py \
         ${statlat} ${statlon} ${OMPcorrectP1} ${OMPcorrectS1} \
        ${TravTable} ${BaziFile} \
        ${P1ArrFile} ${S1ArrFile} \
        $N ${outfile} ${stats_outfile}
