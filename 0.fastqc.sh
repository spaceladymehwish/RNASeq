#!/bin/bash

current_time=$(date "+%Y.%m.%d-%H.%M.%S")

source activate qc

fastqc \
 -o ./results/fastqc \
 ./raw_data/*.fq \
 > ./logs/fastqc.$current_time.out \
 2>./logs/fastqc.$current_time.err
