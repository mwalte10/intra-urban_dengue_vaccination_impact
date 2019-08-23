#!/bin/csh

#$ -N val_
#$ -t 1-5

module load bio/R/3.4.0

Rscript model_validation.R $SGE_TASK_ID