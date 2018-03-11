# create date: Mar 1 2018
# purpose: collapse left and right visual field
# calculate mean error for collapsed conditions


#--------------------------------------------------------------------------------------------

library(tidyverse)

# all subjects were used
subno<-c(1:19,24:36,38,40,41,43:46)
subs<-paste0("Sub",subno)


# "s": input subs; "e": 1, SD method, 2, MAD method; "p": number of SD or number of MAD 
s <- subs
e <- 1
p <- 3


for (subNo in seq_along(s)){
  
  session1<-sprintf('%s_Session1_rawinfo.txt',s[subNo])
  session2<-sprintf('%s_Session2_rawinfo.txt',s[subNo])
  
  files<-list(session1,session2)
  rawdata<-files %>% map(function(x) read.table(x, sep="\t", 
                                                colClasses = c("factor","factor","factor","factor","factor","factor","numeric",
                                                               "numeric","numeric","numeric","numeric","numeric","numeric"))) %>%
    reduce(rbind)
  
  colnames(rawdata)<-c("subject","session","block","type","loc","val","response","oriLeft","oriRight","RT","ISI","error90","error180")
  
  if (e == 1){
    # to remove those with longer response time and higher errors (pSD away from the mean and shorter than 0.1s)
    rawdata$abs90<-abs(rawdata$error90)
    
    attach(rawdata)
    
    aver_RT<-mean(RT)
    sd_RT<-sd(RT)
    outlier_RT<-aver_RT+p*sd_RT
    detach(rawdata)
    data_RT<-rawdata %>% filter(RT<outlier_RT & RT>0.1)
    
    attach(data_RT)
    aver_error<-mean(abs90)
    sd_error<-sd(abs90)
    outlier_error<-aver_error+p*sd_error
    detach(data_RT)
    AftRemOutli<-data_RT %>% filter(abs90<outlier_error)
  } else if (e == 2) {
    # remove trials using MAD method
    rawdata$abs90<-abs(rawdata$error90)
    attach(rawdata)
    med_RT<-median(RT)
    mad_RT<-mad(RT)
    outlier_RT<-med_RT+p*mad_RT
    detach(rawdata)
    data_RT<-rawdata %>% filter(RT<outlier_RT & RT>0.1)
    
    attach(data_RT)
    med_error<-median(abs90)
    mad_error<-mad(abs90)
    outlier_error<-med_error+p*mad_error
    
    detach(data_RT)
    AftRemOutli<-data_RT %>% filter(abs90<outlier_error)
  }
  
  
  ## con100
  con100 <- AftRemOutli %>% filter(val==-1) %>% group_by(ISI) %>% summarise(m_error=mean(abs90))
 
  write.table(con100,sprintf('Con100_sub%d.txt',subno[subNo]),sep="\t",row.names=FALSE,col.names = FALSE)
  
  
  ## con50_valid
  con50_val <- AftRemOutli %>% filter(type==2,val==1) %>% group_by(ISI) %>% summarise(m_error=mean(abs90))
  
  write.table(con50_val,sprintf('Con50_valid_sub%d.txt',subno[subNo]),sep="\t",row.names=FALSE,col.names = FALSE)
  
  
  ## con50_invalid
  con50_inval <- AftRemOutli %>% filter(type==2,val==0) %>% group_by(ISI) %>% summarise(m_error=mean(abs90))

  write.table(con50_inval,sprintf('Con50_invalid_sub%d.txt',subno[subNo]),sep="\t",row.names=FALSE,col.names = FALSE)
}

#----------------------------------------------------------------------------------------------------
# move files
library(fs)
file_move(dir_ls(".",glob="*sub*.txt"),collapseLR)

#----------------------------------------------------------------------------------------------------
# calculating mean error for each condition
# for con100

subno<-c(1:19,24:36,38,40,41,43:46)
allsubs<-paste0("sub",subno)

condition_list <- c("Con100_","Con50_valid_","Con50_invalid_")
con_mean <- data.frame(con100 = rep(0,length(subno)),
                       con50_valid = rep(0,length(subno)),
                       con50_invalid = rep(0,length(subno)))

for(con in seq_along(condition_list)){
  allfiles <- as.list(paste0(paste0(condition_list[con],allsubs),'.txt'))
  con_mean[,con] <- allfiles %>% map(read.table) %>% map(colMeans) %>% map('V2') %>% reduce(rbind)
}
con_mean$subno <- subno
write.csv(con_mean,'CollapsedLR_mean.csv',col.names = TRUE,row.names = FALSE)
