#!/bin/sh

#SBATCH --account=emm3
#SBATCH --job-name=dada2
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --output=data/logs/1_dada2_%J.out
#SBATCH --error=data/logs/1_dada2_%J.err

### DADA2 pipeline ######
## ~~ Trimming, error generation and DADA2 run ~~ ##

#Trimming of the reads (and filtering with maxEE),
#generation of an error model whit the probabilities
#transitions between bases (A -> T) and dada2 algorithm.
#BEFORE continue the pipeline, PLEASE, check the error plot
#(a .pdf generated in the output dir).

## ARGS ##

#[1] /input dataset/, with the fastq.gz (the name of the samples
#    has to be present on the left and separated by an underscore.
#    (multiple underscores are OK as long as the name is on the left).
#    Example: sample120412-g20-mallorca_R1.fastq.gz
#    Sample name will be sample120412-g20-mallorca, which is good bc it will be an unique identifier
#    But maybe you want a cleaner one. Your decision!

#[2] /output dir/ Usually it should be in data/, since afterwards you will use it
#                 for statistical analysis
#                 Best option, tidy option

#[3] /name/ A common identifier to be sure that the output is the
#            correct one. You will thank us that :^)

#[4] /trim position/ The positions to trim, both at the left and the right. The numbers have
#                     to be separated by a comma. Example: 240 (forward), 210 (rev).
#                     WATCH OUT. If you trim too close to your read length, you will lose many reads.
#                     Coded right now for the paired end case. If this is not your case,
#                     say it and we will see what can we do.

#[5] /max expected error/ The tolerated threshold of maxEE. Too stringent, the reads don't pass. Too
#                         relaxed, the modeled error can be incorrect. Usually, a maxEE of 1,2 if the
#                         read quality is really good, or a 2,4.

#[6] /minimum overlap reads/ The minimum expected overlap of reads. The default is 15, but in case you fail to merge
#                            most of the reads, an option is to lower the value to see the incremental changes. Beware not to
#                            go below 8.


module load gcc
module load R

Rscript src/data/dada2/01_dada2-error-output.R \
        data/raw/bbmo/4020 \
        data/cleaned/dada2 \
        bbmo_4020 \
        230,220 \
        2,6 \
        15

Rscript src/data/dada2/01_dada2-error-output.R \
        data/raw/bbmo/5001 \
        data/cleaned/dada2 \
        bbmo_5001 \
        250,210 \
        2,5 \
        15


