library(readxl)
library(stringr)

args = commandArgs(trailingOnly=TRUE)

# needs (JAMS) excel file with sheets for abundance, metadata, and lkt names AND comp analys file from Transnet



# abundance

# metadata file

# lkt name

# Transnet comps file


# usage #
# cd /data/Goldszmid_Lab/analysis/collab_evangelos/JAMS/JAMSbeta/manual_plots/TransNet/analys/lefse
# give full path as argument
# Rscript /data/rodriguesrr/scripts/R/prep_JAMS_for_lefse_w_comp_file.R ./Evangelos_Metadata_and_Relabund_PPM_light_2022-02-04_.xlsx /data/Goldszmid_Lab/analysis/collab_evangelos/JAMS/JAMSbeta/manual_plots/TransNet/data/analys_file_comp_3_18_2022.xlsx



#TransNet/analys/lefse 


indir = getwd()
#--- 
# list the excel jam file
#--- 
all_files = args[1]
all_files



#--- 
# comp_file
#--- 
comp_file_df = read_excel(args[2], na = "N_A") # tibble
#head(comp_file_df)

# remove commented analysis and renumber analysis
comp_file_df = data.frame(comp_file_df)
comp_file_df = comp_file_df[which(!grepl("#", comp_file_df[,1])), ]
rownames(comp_file_df) = NULL
head(comp_file_df)






BIGLIST = c()

for (f in all_files){
  print(f)
  #Sheets = excel_sheets(f)
  #print(Sheets)


  #--- 
  # create output directory with all needed files and for results
  #---
  outdir = str_split(f, "/")[[1]]
  outdir = outdir[length(outdir)]
  outdir = str_replace(outdir, ".xlsx", "")
  print(outdir)
  dir.create(outdir)
  
  
  #--- 
  # metadata
  #--- 
  Metadata_df = read_excel(f, sheet = "Metadata", na = "N_A") # tibble
  Metadata_df = data.frame(Metadata_df)

  
  #--- 
  # LKT abundance
  #--- 
  LKT_PPM_df = read_excel(f, sheet = "LKT_PPM", na = "N_A") # tibble
  LKT_PPM_df = data.frame(LKT_PPM_df)
  colnames(LKT_PPM_df)[1] = c("#OTU ID")
  #print(head(LKT_PPM_df))
  
  
  #--- 
  # lkt names
  #--- 
  LKT_featuretable_df = read_excel(f, sheet = "LKT_featuretable", na = "N_A") # tibble
  LKT_featuretable_df = data.frame(LKT_featuretable_df)
  #print(head(LKT_featuretable_df))
  
  taxonomy = paste(LKT_featuretable_df[, "Kingdom"], LKT_featuretable_df[, "Phylum"],  
                LKT_featuretable_df[, "Class"], LKT_featuretable_df[, "Order"],
                LKT_featuretable_df[, "Family"], LKT_featuretable_df[, "Genus"], 
                LKT_featuretable_df[, "Species"], sep="; ")
  #print(head(taxonomy))
  
  name_info = cbind(LKT_featuretable_df[, "row.names"], taxonomy)
  colnames(name_info) = c("#OTU ID", "taxonomy")
  print(head(name_info))
  #print(dim(name_info))
  
  
  #--- 
  # otu table for prep for lefse command
  #--- 
  otuTable = merge(LKT_PPM_df, name_info, by="#OTU ID", all.x=T)
  print(head(otuTable))
  print(dim(LKT_PPM_df))
  print(dim(otuTable))
  
  
  
  # ----------------
  # subset as per comps file
  # ----------------
  
  for (numb in 1:nrow(comp_file_df)){
      comp = comp_file_df[numb,]
      print(comp)
      ColID = comp[, "ColumID"]
      GrpA = comp[, "Group_A"]
      GrpB = comp[, "Group_B"]
      
      subdir = paste0(c("Analy", numb, comp), collapse="__")
      odir = paste(outdir, subdir, sep="/")
      print(odir)
      dir.create(file.path(outdir,subdir), recursive = TRUE)


      samples_to_keep = Metadata_df[which(Metadata_df[, ColID] == GrpA | Metadata_df[, ColID] == GrpB), 1]
      print(samples_to_keep)
      
      # subset the map file and the abundance file
      keep_map_df = Metadata_df[which(Metadata_df[,1] %in% samples_to_keep), c("row.names", ColID)]
      print(head(keep_map_df))
      write.table(keep_map_df, paste0(odir, "/map.txt"), row.names=F, quote=F, sep="\t")

      keep_abund_df = otuTable[, c("#OTU ID", samples_to_keep, "taxonomy")]
      print(head(keep_abund_df))
      
      ofname = paste0(odir, "/", subdir, ".tmp.tsv")
      write.table(keep_abund_df, ofname, row.names=F, quote=F, sep="\t")
      
      outfname = paste0(odir, "/", subdir, ".tsv")
      comand = paste("echo '# Constructed from biom file' >", outfname)
      print(comand)
      system(comand)
      
      comand2 = paste("cat", ofname, ">>", outfname)
      print(comand2)
      system(comand2)
      

  
      #--- 
      # commands for running lefse pipeline
      #--- 
      cmd1 = paste("cd", paste0(indir, "/", odir))
      cmd2 = paste("Rscript /data/rodriguesrr/scripts/R/prep_qiime2_for_lefse.R", paste0(subdir, ".tsv"), "map.txt 1")
      cmd3 = paste("chmod +x", paste0(subdir, ".tsv.lefse.sh"))
      cmd4 = paste0("./", paste0(subdir, ".tsv.lefse.sh"))
      
      BIGLIST = c(BIGLIST, cmd1, cmd2, cmd3, cmd4)
  
  }
}

print(BIGLIST)
write(BIGLIST, "launch.cmds.for.lefse.sh", sep="\n")
