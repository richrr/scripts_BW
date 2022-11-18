
#ml R/4.1
# usage: Rscript /data/rodriguesrr/scripts/R/html_table_to_csv.R index.html

args = commandArgs(trailingOnly=TRUE)

library(rvest)

infile = args[1]

# https://stackoverflow.com/questions/32400916/convert-html-tables-to-r-data-frame
df = as.data.frame( read_html(infile) %>%
  html_element("table") %>%
  html_table() )

write.csv(df , file = paste0(infile, "_df.csv"), row.names=F, quote=F)

#colnames(df)[2]
#tolower(df[1,2])
