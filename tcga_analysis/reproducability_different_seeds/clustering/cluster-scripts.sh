#!/usr/bin/env bash

bsub -n 100 -W 120:00 -R "rusage[mem=2048]" < cluster-data1_1.R
bsub -n 100 -W 120:00 -R "rusage[mem=2048]" < cluster-data1_2.R
bsub -n 100 -W 120:00 -R "rusage[mem=2048]" < cluster-data1_3.R
bsub -n 100 -W 120:00 -R "rusage[mem=2048]" < cluster-data1_4.R
bsub -n 100 -W 120:00 -R "rusage[mem=2048]" < cluster-data1_5.R

bsub -n 100 -W 120:00 -R "rusage[mem=2048]" < cluster-data2_1.R
bsub -n 100 -W 120:00 -R "rusage[mem=2048]" < cluster-data2_2.R
bsub -n 100 -W 120:00 -R "rusage[mem=2048]" < cluster-data2_3.R
bsub -n 100 -W 120:00 -R "rusage[mem=2048]" < cluster-data2_4.R
bsub -n 100 -W 120:00 -R "rusage[mem=2048]" < cluster-data2_5.R

bsub -n 100 -W 120:00 -R "rusage[mem=2048]" < cluster-data3_1.R
bsub -n 100 -W 120:00 -R "rusage[mem=2048]" < cluster-data3_2.R
bsub -n 100 -W 120:00 -R "rusage[mem=2048]" < cluster-data3_3.R
bsub -n 100 -W 120:00 -R "rusage[mem=2048]" < cluster-data3_4.R
bsub -n 100 -W 120:00 -R "rusage[mem=2048]" < cluster-data3_5.R
