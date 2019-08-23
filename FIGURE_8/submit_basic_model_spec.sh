#!/bin/csh

#$ -N basic_
#$ -t 1-420

module load bio/R/3.4.0

Rscript basic_model_spec.R $SGE_TASK_ID