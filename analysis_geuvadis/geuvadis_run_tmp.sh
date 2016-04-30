#!/bin/bash

RCODE=/home/gosia/R/drimseq_paper/analysis_geuvadis
RWD=/home/Shared/data/seq/geuvadis
ROUT=$RWD/Rout

population='YRI'

for chr in {1..2}
do 

  echo "${population}_${chr}"

  R32loc CMD BATCH --no-save --no-restore "--args rwd='$RWD' workers=1 population='${population}' chr=${chr}" $RCODE/geuvadis_drimseq_0_3_3_run_permutations.R $ROUT/geuvadis_drimseq_0_3_3_run_permutations_${population}_chr${chr}.Rout
  
  tail $ROUT/geuvadis_drimseq_0_3_3_run_permutations_${population}_chr${chr}.Rout
done




RCODE=/home/gosia/R/drimseq_paper/analysis_geuvadis
RWD=/home/Shared/data/seq/geuvadis
ROUT=$RWD/Rout

population='YRI'

for chr in {3..4}
do 

  echo "${population}_${chr}"

  R32loc CMD BATCH --no-save --no-restore "--args rwd='$RWD' workers=1 population='${population}' chr=${chr}" $RCODE/geuvadis_drimseq_0_3_3_run_permutations.R $ROUT/geuvadis_drimseq_0_3_3_run_permutations_${population}_chr${chr}.Rout
  
  tail $ROUT/geuvadis_drimseq_0_3_3_run_permutations_${population}_chr${chr}.Rout

done


RCODE=/home/gosia/R/drimseq_paper/analysis_geuvadis
RWD=/home/Shared/data/seq/geuvadis
ROUT=$RWD/Rout

population='YRI'

for chr in {5..6}
do 

  echo "${population}_${chr}"

  R32loc CMD BATCH --no-save --no-restore "--args rwd='$RWD' workers=1 population='${population}' chr=${chr}" $RCODE/geuvadis_drimseq_0_3_3_run_permutations.R $ROUT/geuvadis_drimseq_0_3_3_run_permutations_${population}_chr${chr}.Rout
  
  tail $ROUT/geuvadis_drimseq_0_3_3_run_permutations_${population}_chr${chr}.Rout
done




RCODE=/home/gosia/R/drimseq_paper/analysis_geuvadis
RWD=/home/Shared/data/seq/geuvadis
ROUT=$RWD/Rout

population='YRI'

for chr in {7..8}
do 

  echo "${population}_${chr}"

  R32loc CMD BATCH --no-save --no-restore "--args rwd='$RWD' workers=1 population='${population}' chr=${chr}" $RCODE/geuvadis_drimseq_0_3_3_run_permutations.R $ROUT/geuvadis_drimseq_0_3_3_run_permutations_${population}_chr${chr}.Rout
  
  tail $ROUT/geuvadis_drimseq_0_3_3_run_permutations_${population}_chr${chr}.Rout

done



RCODE=/home/gosia/R/drimseq_paper/analysis_geuvadis
RWD=/home/Shared/data/seq/geuvadis
ROUT=$RWD/Rout

population='YRI'

for chr in {9..10}
do 

  echo "${population}_${chr}"

  R32loc CMD BATCH --no-save --no-restore "--args rwd='$RWD' workers=1 population='${population}' chr=${chr}" $RCODE/geuvadis_drimseq_0_3_3_run_permutations.R $ROUT/geuvadis_drimseq_0_3_3_run_permutations_${population}_chr${chr}.Rout
  
  tail $ROUT/geuvadis_drimseq_0_3_3_run_permutations_${population}_chr${chr}.Rout

done



RCODE=/home/gosia/R/drimseq_paper/analysis_geuvadis
RWD=/home/Shared/data/seq/geuvadis
ROUT=$RWD/Rout

population='YRI'

for chr in {11..12}
do 

  echo "${population}_${chr}"

  R32loc CMD BATCH --no-save --no-restore "--args rwd='$RWD' workers=1 population='${population}' chr=${chr}" $RCODE/geuvadis_drimseq_0_3_3_run_permutations.R $ROUT/geuvadis_drimseq_0_3_3_run_permutations_${population}_chr${chr}.Rout
  
  tail $ROUT/geuvadis_drimseq_0_3_3_run_permutations_${population}_chr${chr}.Rout
done



RCODE=/home/gosia/R/drimseq_paper/analysis_geuvadis
RWD=/home/Shared/data/seq/geuvadis
ROUT=$RWD/Rout

population='YRI'

for chr in {13..14}
do 

  echo "${population}_${chr}"

  R32loc CMD BATCH --no-save --no-restore "--args rwd='$RWD' workers=1 population='${population}' chr=${chr}" $RCODE/geuvadis_drimseq_0_3_3_run_permutations.R $ROUT/geuvadis_drimseq_0_3_3_run_permutations_${population}_chr${chr}.Rout
  
  tail $ROUT/geuvadis_drimseq_0_3_3_run_permutations_${population}_chr${chr}.Rout
done



RCODE=/home/gosia/R/drimseq_paper/analysis_geuvadis
RWD=/home/Shared/data/seq/geuvadis
ROUT=$RWD/Rout

population='YRI'

for chr in {15..16}
do 

  echo "${population}_${chr}"

  R32loc CMD BATCH --no-save --no-restore "--args rwd='$RWD' workers=1 population='${population}' chr=${chr}" $RCODE/geuvadis_drimseq_0_3_3_run_permutations.R $ROUT/geuvadis_drimseq_0_3_3_run_permutations_${population}_chr${chr}.Rout
  
  tail $ROUT/geuvadis_drimseq_0_3_3_run_permutations_${population}_chr${chr}.Rout
done




RCODE=/home/gosia/R/drimseq_paper/analysis_geuvadis
RWD=/home/Shared/data/seq/geuvadis
ROUT=$RWD/Rout

population='YRI'

for chr in {17..18}
do 

  echo "${population}_${chr}"

  R32loc CMD BATCH --no-save --no-restore "--args rwd='$RWD' workers=1 population='${population}' chr=${chr}" $RCODE/geuvadis_drimseq_0_3_3_run_permutations.R $ROUT/geuvadis_drimseq_0_3_3_run_permutations_${population}_chr${chr}.Rout
  
  tail $ROUT/geuvadis_drimseq_0_3_3_run_permutations_${population}_chr${chr}.Rout

done



RCODE=/home/gosia/R/drimseq_paper/analysis_geuvadis
RWD=/home/Shared/data/seq/geuvadis
ROUT=$RWD/Rout

population='YRI'

for chr in {19..20}
do 

  echo "${population}_${chr}"

  R32loc CMD BATCH --no-save --no-restore "--args rwd='$RWD' workers=1 population='${population}' chr=${chr}" $RCODE/geuvadis_drimseq_0_3_3_run_permutations.R $ROUT/geuvadis_drimseq_0_3_3_run_permutations_${population}_chr${chr}.Rout

  
  tail $ROUT/geuvadis_drimseq_0_3_3_run_permutations_${population}_chr${chr}.Rout
done



RCODE=/home/gosia/R/drimseq_paper/analysis_geuvadis
RWD=/home/Shared/data/seq/geuvadis
ROUT=$RWD/Rout

population='YRI'

for chr in {21..22}
do 

  echo "${population}_${chr}"

  R32loc CMD BATCH --no-save --no-restore "--args rwd='$RWD' workers=1 population='${population}' chr=${chr}" $RCODE/geuvadis_drimseq_0_3_3_run_permutations.R $ROUT/geuvadis_drimseq_0_3_3_run_permutations_${population}_chr${chr}.Rout
  
  tail $ROUT/geuvadis_drimseq_0_3_3_run_permutations_${population}_chr${chr}.Rout
  
done


### CEU

RCODE=/home/gosia/R/drimseq_paper/analysis_geuvadis
RWD=/home/Shared/data/seq/geuvadis
ROUT=$RWD/Rout

population='CEU'

for chr in {1..2}
do 

  echo "${population}_${chr}"

  R32loc CMD BATCH --no-save --no-restore "--args rwd='$RWD' workers=1 population='${population}' chr=${chr}" $RCODE/geuvadis_drimseq_0_3_3_run_permutations.R $ROUT/geuvadis_drimseq_0_3_3_run_permutations_${population}_chr${chr}.Rout
  
  tail $ROUT/geuvadis_drimseq_0_3_3_run_permutations_${population}_chr${chr}.Rout
done




RCODE=/home/gosia/R/drimseq_paper/analysis_geuvadis
RWD=/home/Shared/data/seq/geuvadis
ROUT=$RWD/Rout

population='CEU'

for chr in {3..4}
do 

  echo "${population}_${chr}"

  R32loc CMD BATCH --no-save --no-restore "--args rwd='$RWD' workers=1 population='${population}' chr=${chr}" $RCODE/geuvadis_drimseq_0_3_3_run_permutations.R $ROUT/geuvadis_drimseq_0_3_3_run_permutations_${population}_chr${chr}.Rout
  
  tail $ROUT/geuvadis_drimseq_0_3_3_run_permutations_${population}_chr${chr}.Rout

done


RCODE=/home/gosia/R/drimseq_paper/analysis_geuvadis
RWD=/home/Shared/data/seq/geuvadis
ROUT=$RWD/Rout

population='CEU'

for chr in {5..6}
do 

  echo "${population}_${chr}"

  R32loc CMD BATCH --no-save --no-restore "--args rwd='$RWD' workers=1 population='${population}' chr=${chr}" $RCODE/geuvadis_drimseq_0_3_3_run_permutations.R $ROUT/geuvadis_drimseq_0_3_3_run_permutations_${population}_chr${chr}.Rout
  
  tail $ROUT/geuvadis_drimseq_0_3_3_run_permutations_${population}_chr${chr}.Rout
done




RCODE=/home/gosia/R/drimseq_paper/analysis_geuvadis
RWD=/home/Shared/data/seq/geuvadis
ROUT=$RWD/Rout

population='CEU'

for chr in {7..8}
do 

  echo "${population}_${chr}"

  R32loc CMD BATCH --no-save --no-restore "--args rwd='$RWD' workers=1 population='${population}' chr=${chr}" $RCODE/geuvadis_drimseq_0_3_3_run_permutations.R $ROUT/geuvadis_drimseq_0_3_3_run_permutations_${population}_chr${chr}.Rout
  
  tail $ROUT/geuvadis_drimseq_0_3_3_run_permutations_${population}_chr${chr}.Rout

done



RCODE=/home/gosia/R/drimseq_paper/analysis_geuvadis
RWD=/home/Shared/data/seq/geuvadis
ROUT=$RWD/Rout

population='CEU'

for chr in {9..10}
do 

  echo "${population}_${chr}"

  R32loc CMD BATCH --no-save --no-restore "--args rwd='$RWD' workers=1 population='${population}' chr=${chr}" $RCODE/geuvadis_drimseq_0_3_3_run_permutations.R $ROUT/geuvadis_drimseq_0_3_3_run_permutations_${population}_chr${chr}.Rout
  
  tail $ROUT/geuvadis_drimseq_0_3_3_run_permutations_${population}_chr${chr}.Rout

done



RCODE=/home/gosia/R/drimseq_paper/analysis_geuvadis
RWD=/home/Shared/data/seq/geuvadis
ROUT=$RWD/Rout

population='CEU'

for chr in {11..12}
do 

  echo "${population}_${chr}"

  R32loc CMD BATCH --no-save --no-restore "--args rwd='$RWD' workers=1 population='${population}' chr=${chr}" $RCODE/geuvadis_drimseq_0_3_3_run_permutations.R $ROUT/geuvadis_drimseq_0_3_3_run_permutations_${population}_chr${chr}.Rout
  
  tail $ROUT/geuvadis_drimseq_0_3_3_run_permutations_${population}_chr${chr}.Rout
done



RCODE=/home/gosia/R/drimseq_paper/analysis_geuvadis
RWD=/home/Shared/data/seq/geuvadis
ROUT=$RWD/Rout

population='CEU'

for chr in {13..14}
do 

  echo "${population}_${chr}"

  R32loc CMD BATCH --no-save --no-restore "--args rwd='$RWD' workers=1 population='${population}' chr=${chr}" $RCODE/geuvadis_drimseq_0_3_3_run_permutations.R $ROUT/geuvadis_drimseq_0_3_3_run_permutations_${population}_chr${chr}.Rout
  
  tail $ROUT/geuvadis_drimseq_0_3_3_run_permutations_${population}_chr${chr}.Rout
done



RCODE=/home/gosia/R/drimseq_paper/analysis_geuvadis
RWD=/home/Shared/data/seq/geuvadis
ROUT=$RWD/Rout

population='CEU'

for chr in {15..16}
do 

  echo "${population}_${chr}"

  R32loc CMD BATCH --no-save --no-restore "--args rwd='$RWD' workers=1 population='${population}' chr=${chr}" $RCODE/geuvadis_drimseq_0_3_3_run_permutations.R $ROUT/geuvadis_drimseq_0_3_3_run_permutations_${population}_chr${chr}.Rout
  
  tail $ROUT/geuvadis_drimseq_0_3_3_run_permutations_${population}_chr${chr}.Rout
done




RCODE=/home/gosia/R/drimseq_paper/analysis_geuvadis
RWD=/home/Shared/data/seq/geuvadis
ROUT=$RWD/Rout

population='CEU'

for chr in {17..18}
do 

  echo "${population}_${chr}"

  R32loc CMD BATCH --no-save --no-restore "--args rwd='$RWD' workers=1 population='${population}' chr=${chr}" $RCODE/geuvadis_drimseq_0_3_3_run_permutations.R $ROUT/geuvadis_drimseq_0_3_3_run_permutations_${population}_chr${chr}.Rout
  
  tail $ROUT/geuvadis_drimseq_0_3_3_run_permutations_${population}_chr${chr}.Rout

done



RCODE=/home/gosia/R/drimseq_paper/analysis_geuvadis
RWD=/home/Shared/data/seq/geuvadis
ROUT=$RWD/Rout

population='CEU'

for chr in {19..20}
do 

  echo "${population}_${chr}"

  R32loc CMD BATCH --no-save --no-restore "--args rwd='$RWD' workers=1 population='${population}' chr=${chr}" $RCODE/geuvadis_drimseq_0_3_3_run_permutations.R $ROUT/geuvadis_drimseq_0_3_3_run_permutations_${population}_chr${chr}.Rout

  
  tail $ROUT/geuvadis_drimseq_0_3_3_run_permutations_${population}_chr${chr}.Rout
done



RCODE=/home/gosia/R/drimseq_paper/analysis_geuvadis
RWD=/home/Shared/data/seq/geuvadis
ROUT=$RWD/Rout

population='CEU'

for chr in {21..22}
do 

  echo "${population}_${chr}"

  R32loc CMD BATCH --no-save --no-restore "--args rwd='$RWD' workers=1 population='${population}' chr=${chr}" $RCODE/geuvadis_drimseq_0_3_3_run_permutations.R $ROUT/geuvadis_drimseq_0_3_3_run_permutations_${population}_chr${chr}.Rout
  
  tail $ROUT/geuvadis_drimseq_0_3_3_run_permutations_${population}_chr${chr}.Rout
  
done



###################################################
## Run sQTLSeekeR analysis
###################################################

RCODE=/home/gosia/R/drimseq_paper/analysis_geuvadis
RWD=/home/Shared/data/seq/geuvadis
ROUT=$RWD/Rout

population='YRI'

  R31 CMD BATCH --no-save --no-restore "--args rwd='$RWD' workers=1 population='${population}'" $RCODE/geuvadis_sqtlseeker_2_1.R $ROUT/geuvadis_sqtlseeker_2_1_run_${population}.Rout
  
  tail $ROUT/geuvadis_sqtlseeker_2_1_run_${population}.Rout



















