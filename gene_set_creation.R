#Seting working directory
setwd('ENTER DIRECTORY PATH')

#Loading required packages
library(data.table)

#-------------------------------------------
# FUMA Gene Set Generation
#-------------------------------------------

#Download Yengo Height GWAS sumstats: https://portals.broadinstitute.org/collaboration/giant/images/6/63/Meta-analysis_Wood_et_al%2BUKBiobank_2018.txt.gz
#Upload to FUMA using default settings
#Download FUMA MAGMA results
#Use list of sig sets reported by MAGMA to identify genes within each sig set

#Base folder for storing MAGMA results
base_folder <- 'FUMA_job...'
magma <- readLines(paste0(base_folder,'/magma.gsa.sets.genes.out'))
sets <- magma[grepl('VARIABLE = ',magma)]
for(i in 1:length(sets)){
  set_i<-magma[grepl(paste0('^_SET',i,'_'), magma)]
  set_i<-gsub('\\s+',' ',set_i)
  set_i<-do.call(rbind, strsplit(set_i,' '))
  set_i<-data.frame(set_i)
  names(set_i)<-set_i[1,]
  set_i<-set_i[-1,]
  set_file<-data.frame(CHR=paste0('chr',set_i$CHR),
                    START=set_i$START,
                    STOP=set_i$STOP,
					ID=set_i$GENE,
					Set=gsub(' .*','',gsub('.*= ','',sets[i])))
  write.table(set_file, paste0(base_folder, '/Set_',i,'.txt'), col.names=F, row.names=F, quote=F)
}
