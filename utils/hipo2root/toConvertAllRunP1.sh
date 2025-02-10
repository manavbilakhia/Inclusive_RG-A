#!/bin/csh -f
source /group/clas12/packages/setup.csh
module avail
module load clas12/pro
clas12root -q -b ana12GeVShortFCQA.C  --in=allRunsP1NickPart_2023.dat >> logs/logNorm_allRunsP1NickPart_2023
