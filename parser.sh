#!/usr/bin/env bash

#$ -P prod
#$ -N CHIP_PARSER
#$ -A 220131_PESA_HLA
#$ -l thread=16
#$ -l h_vmem=16G
#$ -o variant_per_sample_parser.stdout
#$ -e variant_per_sample_parser.stderr

#Juan Ignacio Alvarez Arenas

. /data3/220131_PESA_HLA/env.sh

cd /data3/220919_CHIP_Annotator/_data

python3 VCF-to-JSON-Parser-variant-per-sample.py

