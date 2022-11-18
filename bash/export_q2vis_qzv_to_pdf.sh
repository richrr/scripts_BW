#!/bin/bash


# cd /data/MicrobiomeCore/analysis/GTAM45_20220714

# usage: bash /data/rodriguesrr/scripts/bash/export_q2vis_qzv_to_pdf.sh pdfdump Group > ResultsDump.csv


# every csv can be read by R and written to table in pdf?


idir=$(pwd)
#FILES=$(find $idir)

output=$1
groupCols=$2  # comma separated list of columns to be added later


#https://unix.stackexchange.com/questions/156534/bash-script-error-with-strings-with-paths-that-have-spaces-and-wildcards
#FILES=$(/bin/ls -d "$idir"/*/q2_visualizations_*/)
#https://stackoverflow.com/questions/246215/how-can-i-generate-a-list-of-files-with-their-absolute-path-in-linux
FILES=$(find "$PWD")

#alphametrics=("chao" "simpson_e" "obs_otus" "shannon" "faith_pd") 

module unload qiime
module load qiime/2-2020.2
module load biom-format

ml R/4.1


# https://stackoverflow.com/questions/6723426/looping-over-arrays-printing-both-index-and-value
for f in ${FILES[@]}
do
  #echo "${f}"
  # https://stackoverflow.com/questions/229551/how-to-check-if-a-string-contains-a-substring-in-bash

  if [[ $f =~ "core_metrics_" ]]; then  #[[ $f =~ "core_metrics_" && $f =~ ".qzv" ]]
   #echo "exporting $f"
   outf=${f/"$idir"/""}


   # beta: 
     # plot for each metric (using phyloseq, rarefied or relativized table?)
     # beta_sginificance (export qzv: and use *pairwise.csv result files; parse the index.html for the global anosim/permanova p-value; not using png/pdf plots)
     #qiime tools export --input-path /data/MicrobiomeCore/analysis/GTAM45_20220714/D5/q2_visualizations_202208082210/beta/beta_sig_Group_20220808221428/bray_curtis_anosim_bsig.qzv --output-path pdfdump//D5/q2_visualizations_202208082210/beta/beta_sig_Group_20220808221428/bray_curtis_anosim_bsig.qzv
     #qiime tools export --input-path /data/MicrobiomeCore/analysis/GTAM45_20220714/D5/q2_visualizations_202208082210/beta/beta_sig_Group_20220808221428/bray_curtis_permanova_bsig.qzv --output-path pdfdump//D5/q2_visualizations_202208082210/beta/beta_sig_Group_20220808221428/bray_curtis_permanova_bsig.qzv

     ### comment this; since this cannot get anosim and permanova results ###
     ### which first need to be exported to be able to parse from index.html ###
     #if [[ $f =~ "_pcoa_results.qza" ]] ; then
      #  echo "File: $f"
      #  cmd1="Rscript /data/rodriguesrr/scripts/R/qiime2R_parse_qza.R $f $groupCols $output"
      #  echo $cmd1
        
      #  break
     #fi

     if [[ $f =~ "_bsig.qzv" ]] ; then
        #echo "File: $f"
        cmd1="qiime tools export --input-path ${f} --output-path $output/${outf}"
        #echo $cmd1
        $cmd1

        #cmd2="grep -A30 'method name' $output/${outf}/index.html"
        cmd2="Rscript /data/rodriguesrr/scripts/R/html_table_to_csv.R $output/${outf}/index.html"
        #echo $cmd2
        #TestName=$($cmd2)
        #echo $TestName
        $cmd2
        echo
        cat "$output/${outf}/index.html_df.csv"
        echo

        echo "pairwise beta_significance results:"
         if [[ $f =~ "permanova" ]] ; then
           cat "$output/${outf}/permanova-pairwise.csv"
         elif [[ $f =~ "anosim" ]] ; then
           cat "$output/${outf}/anosim-pairwise.csv"
         fi
        echo

        cmd3="Rscript /data/rodriguesrr/scripts/R/qiime2R_parse_qza_dump.R $f $idir/ $groupCols $output $output/${outf}/index.html_df.csv"
        #echo $cmd3
        $cmd3
        echo

        #break

     fi

   # alpha: (export qzv and use csv result files; and parse the jsonp to get global pvalue and either the jsonp or metadata.tsv to make boxplot if needed)
     # boxplot for each metric (skip)
     # alpha_significance.qzv
     #qiime tools export --input-path /gpfs/gsfs6/users/MicrobiomeCore/analysis/GTAM45_20220714/D1/core_metrics_20220808229//chao_alpha_significance.qzv --output-path pdfdump/qzv

     if [[ $f =~ "_alpha_significance.qzv" ]] ; then
        #echo "File: $f"
        cmd1="qiime tools export --input-path ${f} --output-path $output/${outf}"
        #echo $cmd1
        $cmd1

        # https://stackoverflow.com/questions/1955505/parsing-json-with-unix-tools
        cmd2="grep -Po '"H":.*?,' $output/${outf}/column-$groupCols.jsonp"
        #echo $cmd2
        Hresult=$(grep -Po '"H":.*?,' "$output/${outf}/column-$groupCols.jsonp")
        #echo $Hresult


        cmd3="grep -Po '"p":.*?,' $output/${outf}/column-$groupCols.jsonp"
        #echo $cmd3
        Presult=$(grep -Po '"p":.*?,' "$output/${outf}/column-$groupCols.jsonp")
        #echo $Presult

        echo "Global alpha_significance result:, $Hresult $Presult"
        echo
        echo "kruskal-wallis-pairwise alpha_significance results:"
        cat "$output/${outf}/kruskal-wallis-pairwise-$groupCols.csv"
        echo
     fi
   #done

  #break
  fi

done





# taxa-bar-plots (csv files), alpha rarefaction (csv and jsonp) skip



# get the csv files and plot as needed?
#qiime tools export --input-path /data/MicrobiomeCore/analysis/GTAM45_20220714/D5/q2_visualizations_202208082210/rarefied_taxa-bar-plots.qzv --output-path pdfdump//D5/q2_visualizations_202208082210/rarefied_taxa-bar-plots.qzv


#qiime tools export --input-path /data/MicrobiomeCore/analysis/GTAM45_20220714/D5/q2_visualizations_202208082210/alpha/alpha_rarefaction.qzv --output-path pdfdump//D5/q2_visualizations_202208082210/alpha/alpha_rarefaction.qzv
