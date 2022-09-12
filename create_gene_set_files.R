#Setting working directory
setwd('ENTER DIRECTORY PATH')

for (set in 1:64) {
  
  cmd_str <- paste0("plink --bfile target_train.QC --make-set Set_", as.character(set), ".txt --make-set-collapse-group --write-set --out snplist_set_", as.character(set))
  system(paste (cmd_str ,sep=" "))
  
}

#Processing set files
for (set in 1:64){
  set_list <- fread(paste0('snplist_set_',as.character(set),'.set'), header=TRUE)
  
  #Rename column
  colnames(set_list) <- 'rsid'
  
  #Remove last 2 redundant rows
  set_list <- head(set_list,-2)
  
  file_name <- paste0('snp_set_',as.character(set),'.txt')
  write.table(set_list, file_name, col.names=F, quote=F, row.names=F)
  
}

