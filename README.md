# single_array_PS_locate
A python program to estimate the location of seismic events given probability distributions for P1 and S1 arrival times and backazimuth from a seismic station.  

[![SQAaaS badge](https://github.com/EOSC-synergy/SQAaaS/raw/master/badges/badges_150x116/badge_software_bronze.png)](https://api.eu.badgr.io/public/assertions/_CpBnNcOS3GALTMMq6Mv3Q "SQAaaS bronze badge achieved")  

[![SQAaaS badge shields.io](https://img.shields.io/badge/sqaaas%20software-bronze-e6ae77)](https://api.eu.badgr.io/public/assertions/_CpBnNcOS3GALTMMq6Mv3Q "SQAaaS bronze badge achieved")  

https://sqaaas.eosc-synergy.eu/full-assessment/report/https://raw.githubusercontent.com/eosc-synergy/single_array_PS_locate.assess.sqaaas/main/.report/assessment_output.json  


The main program is called using  
```
python single_array_PS_locate.py statlat statlon OMPcorrectP1 OMPcorrectS1 TraveltimeTableFile BaziDistFile P1ArrtimeDistFile S1ArrtimeDistFile [N=10000] [outfile=single_array_PS_locate_output.txt] [stats_outfile=single_array_PS_locate_stats_output.txt]
```
and is easiest called using a shell script of the form  
```
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
```

Let's assume that we have a traveltime table file (e.g. **ak135_P1_S1_zerodepth.txt**) with a format as follows:  
```
# Columns: P1_travel_time  S1_travel_time  (S1-P1)  distance_deg  depth_km
0.000000        0.000000        0.000700        0.00    0.00
3.834300        6.427500        2.593200        0.20    0.00
7.668600        12.854900       5.186300        0.40    0.00
11.502900       19.282300       7.779400        0.60    0.00
15.337100       25.709600       10.372500       0.80    0.00
19.171300       32.136900       12.965600       1.00    0.00
23.005400       38.564000       15.558600       1.20    0.00
26.774600       44.991000       18.216400       1.40    0.00
29.525300       50.877900       21.352600       1.60    0.00
```

then, given an S-minus-P traveltime difference, we can calculate the distance from a station using:  
```
python find_dist_from_table.py ak135_P1_S1_zerodepth.txt 32.12
```
and this will calculate the distance in degrees given an S-minus-P time of 32.12 seconds ...  
```
P_travel_time = 43.075901
S_travel_time = 75.195901
S_minus_P     = 32.120000
Distance_deg  = 2.585275
```
This appears to fit well with our table.    

We now also need to define a probability distribution for backazimuth from our station.  
We do this in the form of a single file, e.g. **sig1_backazi_distr.txt**  
```
135.0  141.1
0.1
0.2
0.3
0.4
0.5
0.4
0.3
0.2
0.1
```
where the first line gives the min and max values of backazimuth and the following lines
give that values at equally spaced points in this interval, including the extreme values.
So typing  
```
python test_azi_distr.py sig1_backazi_distr.txt 10000
```
should result in a plot resembling the output in the image below  

![Output from running "python test_azi_distr.py sig1_backazi_distr.txt 10000"](azi_distr_example.png)  

Similarly, if we define a probability distribution for the first P-wave arrival in a single file,
e.g. **ev1_P1_arrtime_distr.txt**  
```
2020-06-22T09:20:01.500 2020-06-22T09:20:02.500
0.1
0.2
0.3
0.4
0.5
0.6
0.6
0.6
0.6
0.7
0.5
0.2
```
where the first line gives the first and last times of our distribution and the next lines define the
probability function for equally spaced intervals between these limits. Now, typing  
```
python test_arrtime_distr.py   ev1_P1_arrtime_distr.txt    10000
```
should give you something like the following output  

![Output from running "python test_arrtime_distr.py   ev1_P1_arrtime_distr.txt    10000"](arrtime_distr_example.png)  

