#!/bin/bash


odir=$1
infile=$2
group=$3
fcfile=$4
corranaly="Analys $5 "
companaly="Analys $6 "
pval=$7
fdr=$8


#based on /data/rodriguesrr/Koltsova/analysis/Jan2022/TransNet/analys/merged/corr/il17-il22-double/p1_cp1_cfdr1/per_analysis/create-nets-bash-cmd.sh

# usage: 
# cd /data/rodriguesrr/Koltsova/analysis/July2022/analys/corr/merged/CrePosKO_CreNegWT/p1_cp1_cfdr1/per_analysis
# bash /data/rodriguesrr/scripts/bash/create-nets-bash-cmd.sh net /data/rodriguesrr/Koltsova/analysis/July2022/analys/corr/merged/CrePosKO_CreNegWT/p1_cp1_cfdr1/merged_sp_FolChMedian_merged-parallel-output.csv Cre /data/rodriguesrr/Koltsova/analysis/July2022/analys/comp/mw_comp-output.csv 1 1 0.15 1

OUTFILE="run-create-nets-cmds.txt"



for consis in $(ls *.txt)
do
  echo $consis
  lines=`wc -l $consis | awk '{ print $1 }'`
  if [ $lines -gt 0 ]
  then
    mkdir -p $odir/$consis
    #cmd="Rscript /data/rodriguesrr/scripts/R/create-network_bw.R --file $infile --group $group --foldchange $fcfile --consistent $consis --indivPvalCutoff 1 --output $odir/$consis/net --analysiscorr '$corranaly' --analysisfc '$companaly' --combPvalCutoff $pval --combFDRCutoff $fdr"
    
    # with  --numbDataFromFCfile
    cmd="Rscript /data/rodriguesrr/scripts/R/create-network_bw.R --file $infile --group $group --foldchange $fcfile --consistent $consis --indivPvalCutoff 1 --output $odir/$consis/net --analysiscorr '$corranaly' --analysisfc '$companaly' --combPvalCutoff $pval --combFDRCutoff $fdr --numbDataFromFCfile"
    echo $cmd >> $OUTFILE
    #$cmd # running it here directly throws errors. its better to run the file conatining all commands
  else
    echo "Skipping $consis because it is empty"
  fi
done
