
# working folder: 'permutated_indsub_con' or in 'regression_analysis' %>% 'shuffleFFT'
# calculate the mean amplitude in theta/delta/alpha band averaged across all iterations for all subjects

# 4- 8Hz (6-11)

fs <- 25
nfft <- 32*1 # refer to 'do_spectral.m'
Hz <- c((0:(nfft/2-1))*fs/nfft)


library(tidyverse)
# -------------------------------------------------------------------------
allsub <- c(1:19,24:36,38,40,41,43:46)
outlier <- c(5,11,18)
sub_valid <- setdiff(allsub,outlier)
valid_sub <- sub_valid # get only valid subjects

subNo <- paste0("sub",valid_sub)
prefix <- "_LP_Smoothed_"
suffix <- "_indshuf_ampli.txt"

#offset_all <- read.csv('./CollapsedLR_mean_nooutlier.csv',header = TRUE)

# function starts here
cond_df<-function(cond,band){
  # if (cond=='Con100'){
  #   offset_index <- 1
  # } else if (cond=='Con50_valid'){
  #   offset_index <- 2
  # } else if (cond=="Con50_invalid"){
  #   offset_index <- 3
  # }
  
  ind_files <- as.list(paste0(subNo,prefix,cond,suffix))
  
  Meanampli <- ind_files %>% map(read.table) %>% map(function(x) x[,band]) %>% 
    map(colMeans) %>% map(mean) %>% reduce(rbind)
  
  #offset_cond <- offset_all[,offset_index]
  
  results_cond <- data_frame(subj = seq_along(valid_sub),
                             ampli = as.vector(Meanampli))#,
                             #offset = offset_cond)
  return(results_cond)
}

# --- do some changes here ------
output_suffix <- 'onefreq'
#----------------------------

# -------------------------- THETA ------------------------------------------------
band <- 10:12 # 7.03 - 8.59hz
# con100
cond <- "Con100"
cond100_df_theta <- cond_df(cond,band)
write.csv(cond100_df_theta,sprintf('cond100_meanShuffle_theta_%s.csv',output_suffix),row.names = FALSE, col.names = TRUE)

# con50_valid
cond <- "Con50_valid"
cond50_valid_df_theta <- cond_df(cond,band)
write.csv(cond50_valid_df_theta,sprintf('cond50_valid_meanShuffle_theta_%s.csv',output_suffix),row.names = FALSE, col.names = TRUE)

# con50_invalid
cond <- "Con50_invalid"
cond50_invalid_df_theta <- cond_df(cond,band)
write.csv(cond50_invalid_df_theta,sprintf('cond50_invalid_meanShuffle_theta_%s.csv',output_suffix),row.names = FALSE, col.names = TRUE)


# --------------------------Delta ------------------------------------------------
# con100
cond <- "Con100"
cond100_df_delta <- cond_df(cond,1:7)
write.csv(cond100_df_delta,'cond100_meanShuffle_delta_narrow.csv',row.names = FALSE, col.names = TRUE)

# con50_valid
cond <- "Con50_valid"
cond50_valid_df_delta <- cond_df(cond,1:7)
write.csv(cond50_valid_df_delta,'cond50_valid_meanShuffle_delta_narrow.csv',row.names = FALSE, col.names = TRUE)

# con50_invalid
cond <- "Con50_invalid"
cond50_invalid_df_delta <- cond_df(cond,1:7)
write.csv(cond50_invalid_df_delta,'cond50_invalid_meanShuffle_delta_narrow.csv',row.names = FALSE, col.names = TRUE)


# Alpha ------------------------------------------------
# con100
cond <- "Con100"
cond100_df_alpha <- cond_df(cond,12:16)
write.csv(cond100_df_alpha,'cond100_meanShuffle_alpha_narrow.csv',row.names = FALSE, col.names = TRUE)

# con50_valid
cond <- "Con50_valid"
cond50_valid_df_alpha <- cond_df(cond,12:16)
write.csv(cond50_valid_df_alpha,'cond50_valid_meanShuffle_alpha_narrow.csv',row.names = FALSE, col.names = TRUE)

# con50_invalid
cond <- "Con50_invalid"
cond50_invalid_df_alpha <- cond_df(cond,12:16)
write.csv(cond50_invalid_df_alpha,'cond50_invalid_meanShuffle_alpha_narrow.csv',row.names = FALSE, col.names = TRUE)

