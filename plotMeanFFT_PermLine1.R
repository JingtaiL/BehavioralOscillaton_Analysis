# update : 03-20-2018
# this script plot mean observed FFT and the 95% threshold line

library(tidyverse)
library(cowplot)


# set important parameters

rm(list = ls())

fs <- 25
nfft <- 32*1 # refer to 'do_spectral.m'
pvalue <- 0.05

filepattern <- c("LP_Smoothed_")

#--------------------New Method to plot data--------------------------------------------
rm(list = ls())

fs <- 25
nfft <- 32*1 # refer to 'do_spectral.m'
pvalue <- 0.05

# calcualte confidence interval for each condition
sub_all <- c(1:19,24:36,38,40,41,43:46)
outlier <- c(5,11,18)
sub_valid <- setdiff(sub_all,outlier)
sub_prefix <- paste0('sub',sub_valid)

# A function to generate plot
plot_meanFFT <- function(cond){
  setwd("C:/Users/psw/Dropbox/2nd project/Gabor_WM_v4/realResults/DataByCondition/collapseLR/collapsedLR_3madwithin_nosmooth/perm500_nfft1_nosmooth_newperm_newp_movave6_3outliers_ver1/individual_fft")
  
  # prepare data for mean FFT (reference: https://stackoverflow.com/questions/32669473/plotting-the-means-with-confidence-intervals-with-ggplot)
  constantfft<-sprintf('_LP_Smoothed_%s_ind_ampli.txt',cond)
  indfft_file <- paste0(sub_prefix,constantfft)
  
  all_load <- indfft_file %>% map(read.table) %>% reduce(rbind)
  all_ampli_hz <- data.frame(Hz = c((0:(nfft/2-1))*fs/nfft), Amplitude = all_load$V1)
  
  # calculated the threshold line (the cutoff of 5% in permutated FFTs)
  setwd("C:/Users/psw/Dropbox/2nd project/Gabor_WM_v4/realResults/DataByCondition/collapseLR/collapsedLR_3madwithin_nosmooth/perm500_nfft1_nosmooth_newperm_newp_movave6_3outliers_ver1/permutated_meansub_con")
  
  permu_data <- read.table(sprintf('./LP_Smoothed_%s_mean_ampli_rand.txt',cond),header=FALSE)
  
  perc5_threshold <- vector()
  for (i in 1:16){
    select_column <- paste0('V',i)
    percent5 <- permu_data %>% select(select_column) %>% top_n(50) %>% top_n(-1)
    if (length(percent5[[1]]) > 1) {
      percent5 <- percent5[[1]][1]
    }
    perc5_threshold[i] <- percent5
  }
  indcon_perm_threshold <- perc5_threshold %>% reduce(cbind)
  thres_df <- data.frame(Hz = c((0:(nfft/2-1))*fs/nfft), Amplitude = as.vector(indcon_perm_threshold))
  
  #  'random_p' folder
  # find out the exact Hz and insert a vertical line
  setwd("C:/Users/psw/Dropbox/2nd project/Gabor_WM_v4/realResults/DataByCondition/collapseLR/collapsedLR_3madwithin_nosmooth/perm500_nfft1_nosmooth_newperm_newp_movave6_3outliers_ver1/random_p")

  con_p <- read.table(sprintf('LP_Smoothed_%s_pvalue_uncor.txt',cond))
  con_p_index <- which(con_p<pvalue) # at 3.90625Hz for con100
  
  
 plot_beauty <- ggplot(all_ampli_hz,aes(x=Hz,y=Amplitude))+
    stat_summary(geom="ribbon", fun.data=mean_cl_normal, 
                 fun.args=list(conf.int=0.95), fill="lightblue", alpha = 0.4)+
    stat_summary(geom="line", fun.y=mean, size=1.5,color="#3D79F37F",alpha=0.7)+
    xlab("Frequency (Hz)")+ylab("Spectral Amplitude")+
    theme_classic()+scale_x_continuous(breaks=seq(0,15,2))+
    theme(axis.title.x=element_text(size=16),
          axis.title.y=element_text(size=16),
          axis.text.x = element_text(size=16),
          axis.text.y = element_text(size=16))+
    geom_line(data=thres_df,aes(x=Hz,y=Amplitude),
              size=0.8,color="#E6352F7F",alpha=0.5,linetype = 2)+
    geom_vline(xintercept = thres_df$Hz[con_p_index],linetype = 3)
 
 return(plot_beauty)
}

# Con100
plot_meanFFT('Con100')

# Con50valid
plot_meanFFT('Con50_valid')

# con50invalid
plot_meanFFT('Con50_invalid')
#------------------------------------------------------------------------
  


#------------------OLD method (without margin)---------------------------------------------

#### plot mean FFT results for all conditions
# in ‘meansub_fft_con' folder
# prepare data
constantfft<-c('Con100_allsub_ampli_mean.txt','Con50_invalid_allsub_ampli_mean.txt',
               'Con50_valid_allsub_ampli_mean.txt')

fftfiles<-paste0(filepattern, constantfft)

all_con<- fftfiles %>% map(read.table) %>% reduce(rbind)

trans_all_con<-data.frame(t(all_con))

colnames(trans_all_con)<-c("con100","50invalid","50valid")

trans_all_con$Hz<-c((0:(nfft/2-1))*fs/nfft)
rownames(trans_all_con)<-NULL

long_all_con<- trans_all_con %>% gather(Conditions,Amplitude,"con100":"50valid")
long_all_con$Conditions<-factor(long_all_con$Conditions, levels=(c("con100","50invalid","50valid")))


# prepare data for each condition
# in 'permutated_meansub_con' folder
# con100 -------------------------------------------------
con100<-long_all_con %>% filter(Conditions=="con100")
data_100 <- read.table('./LP_Smoothed_Con100_mean_ampli_rand.txt',header=FALSE)

perc5_threshold_100 <- vector()
for (i in 1:16){
  select_column <- paste0('V',i)
  percent5 <- data_100 %>% select(select_column) %>% top_n(50) %>% top_n(-1)
  perc5_threshold_100[i] <- percent5
}
indcon_perm_threshold <- perc5_threshold_100 %>% reduce(cbind)


# ' random_p' folder
# find out the exact Hz and insert a vertical line
con100_p <- read.table('LP_Smoothed_Con100_pvalue_uncor.txt')
con100_p_index <- which(con100_p<0.05) # at 3.90625Hz


#plot average FFT results + threshold line + significant frequency line

con100$shuff_amp <- as.vector(indcon_perm_threshold)
plot_100<-ggplot(con100,aes(x=Hz,y=Amplitude))+geom_line(size=1.5,color="#3D79F37F",alpha=0.7)+
  xlab("Frequency (Hz)")+ylab("Spectral Amplitude")+
  theme_classic()+scale_x_continuous(breaks=seq(0,15,2))+
  theme(axis.title.x=element_text(size=16),
        axis.title.y=element_text(size=16),
        axis.text.x = element_text(size=16),
        axis.text.y = element_text(size=16))
plot_100+geom_line(data=con100,aes(x=Hz,y=shuff_amp),
                   size=0.8,color="#E6352F7F",alpha=0.5,linetype = 2)+
  geom_vline(xintercept = con100$Hz[con100_p_index],linetype = 3)


#--------------------------------------------------------------------
# prepare data for each condition
# in 'permutated_meansub_con' folder
# con50 valid -------------------------------------------------
con50_valid<-long_all_con %>% filter(Conditions=="50valid")
setwd("~/Dropbox/2nd project/Gabor_WM_v4/realResults/DataByCondition/collapseLR/nosmooth/perm1000_nfft1_allsub_nosmooth_newperm_newp_movave6/permutated_meansub_con")
data_50valid <- read.table('./LP_Smoothed_Con50_valid_mean_ampli_rand.txt',header=FALSE)

perc5_threshold_50valid <- vector()
for (i in 1:16){
  select_column <- paste0('V',i)
  percent5 <- data_50valid %>% select(select_column) %>% top_n(50) %>% top_n(-1)
  if (length(percent5[[1]]) > 1) {
    percent5 <- percent5[[1]][1]
  }
  perc5_threshold_50valid[i] <- percent5
}
indcon_perm_threshold <- perc5_threshold_50valid %>% reduce(cbind)


# ' random_p' folder
# find out the exact Hz and insert a vertical line
setwd("~/Dropbox/2nd project/Gabor_WM_v4/realResults/DataByCondition/collapseLR/nosmooth/perm1000_nfft1_allsub_nosmooth_newperm_newp_movave6/random_p")
con50valid_p <- read.table('LP_Smoothed_Con50_valid_pvalue_uncor.txt')
con50valid_p_index <- which(con50valid_p<0.05) # at 10.15625Hz in this case


#plot average FFT results + threshold line + significant frequency line

con50_valid$shuff_amp <- as.vector(indcon_perm_threshold)

plot_50valid<-ggplot(con50_valid,aes(x=Hz,y=Amplitude))+geom_line(size=1.5,color="#3D79F37F",alpha=0.7)+
  xlab("Frequency (Hz)")+ylab("Spectral Amplitude")+
  theme_classic()+scale_x_continuous(breaks=seq(0,15,2))+
  theme(axis.title.x=element_text(size=16),
        axis.title.y=element_text(size=16),
        axis.text.x = element_text(size=16),
        axis.text.y = element_text(size=16))
plot_50valid+geom_line(data=con50_valid,aes(x=Hz,y=shuff_amp),
                   size=0.8,color="#E6352F7F",alpha=0.5,linetype = 2)+
  geom_vline(xintercept = con50_valid$Hz[con50valid_p_index],linetype = 3)



# not finished ---------




