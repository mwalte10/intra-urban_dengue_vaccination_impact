#!/bin/csh

#$ -N basic_
#$ -t 1-20

module load bio/R/3.4.0

Rscript basic_model_no_travel.R $SGE_TASK_ID