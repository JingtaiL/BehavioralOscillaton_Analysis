%% update date: 10/6/17
% This program lists every step, including R processing steps


%%

close all
clear all
% in office
main_dir=pwd; % path to main folder containing all data, functions, toolboxes, scripts...
addpath(genpath(main_dir))
% in office
dat_dir=[main_dir '\DataByCondition']; % group data by condition was done by "AccuracyResultReport.Rmd"
% at home
dat_dir=[main_dir '/DataByCondition'];

outlier = [5 11 17 18];
%% SG smooth
cd ../DataByCondition/fs7ord3
order=3;
framelen=7;
prefix='SG'; %generated files will be named 'prefix_Smooth_Conxxx'
do_SGsmooth(outlier,order,framelen,prefix)

%% low-pass filter
cd ../DataByCondition/low_pass-8Hz
%outlier = [];
cutoff=12; % cut-off frequency
prefix='LP'; %generated files will be named 'prefix_Smooth_Conxxx'
order = 1;
do_lowpass(outlier,cutoff,order,prefix)

%% band-pass filter
cd ../DataByCondition/low_pass-8Hz
%outlier = [];
cutoff=[1 10]; % cut-off frequency
prefix='BP'; %generated files will be named 'prefix_Smooth_Conxxx'
order = 2;
do_bandpass(outlier,cutoff,order,prefix)

%% do spectral analysis
% still in smooth data folder
% if the prefix is SG, don't forget to check the `do_spectral` file to change
% name
nfft_times = 2;
lp_cutoff = 12;
lp_order = 2;
windows = 1;
outputfile = sprintf('perm1000_nfft%s_2mad_noscale_LPsmooth_newperm_newp',num2str(nfft_times));
do_spectral2(outlier,nfft_times,1000,outputfile,windows,lp_cutoff,lp_order) % newp accounts for p=0

parfol = 'phase';
do_phase(outlier,nfft_times,parfol)
%% phase concentration in 100L/R condition
% wd: phase
con = 1; % left
con = 2; % right
index = 14; %3.3854Hz
index = 20; % 4.9479Hz

[ind_phase,mean_phase,p] = phase_concentration(outlier,con,index)

p_all = [];
for i = 1:48
    [ind_phase,mean_phase,p] = phase_concentration(outlier,con,i);
    p_all(i) = p;
end
[freq_Axis',p_all']
    
rad = circ_ang2rad(ind_phase);
figure(5)
circ_plot(rad','hist',[],20,true,true,'linewidth',2,'color','r')

%% extract significant frequencies and plot 
% this step used 'find_sigHz_plot.Rmd', to use this script properly,
% remember to put both script and processing files in the same folder. 

%% plot raw against smoothed
% this step used 'Raw_Smooth_plot.Rmd'

%% doing circular statistics
% working directory is the smooth data folder
%con = 1; % left VF
con = 2; % right VF
times =3;
smooth = 'lp';
outlier = [5 11 17 18];

% run different frequency bands in for-loop
freqs = [0:0.2604:8];
Allfreqs_ind_phaseDiff = zeros(length(freqs),35);
Allfreqs_mean_phaseDiff = zeros(length(freqs),1);
Allfreqs_p = zeros(length(freqs),1);

for f = 1:length(freqs)
    Allfreqs_phaseDiff(f,:) = do_calc_phase_diff(outlier,con,times,smooth,freqs(f));
    deg = Allfreqs_phaseDiff(f,:)';
    rad = circ_ang2rad(deg); % convert to radians
    mean_rad = circ_mean(rad); % calculate the mean
    mean_deg = circ_rad2ang(mean_rad); % convert back to degrees
    if mean_deg < 0
        mean_deg = circ_rad2ang(2*pi+mean_rad);
    else
        mean_deg = mean_deg;
    end
    p_uniform = circ_rtest(rad);
    
    Allfreqs_mean_phaseDiff(f) = mean_deg;
    Allfreqs_p(f) = p_uniform;
end

list_allfreqs = [freqs;Allfreqs_mean_phaseDiff';Allfreqs_p'];
save_allfreqs=sprintf('right_allfreqs_phasediff.txt');
openfile=fopen(save_allfreqs,'w');
fprintf(openfile,'%f %f\n',list_allfreqs);
fclose(openfile);

%%
targetHz = [2.34]; % left
targetHz = [5.469 5.989]; % left
targetHz = [2.34 5.989];
targetHz = [2.86 3.64];
targetHz = [0.5 4];
targetHz = [2 4];
targetHz = [4 8];
% good_subjects = [1 2 7 8 13 16 19 25 28 30 33 36 38 40 43 44 45];
% bad_subjects = [3 4 6 9 10 12 14 24 26 27 29 31 32 34 35 41 46];
% outlier = [outlier,bad_subjects];

%%
AllSubPhase = do_calc_phase_diff(outlier,con,times,smooth,targetHz);

deg=AllSubPhase';

rad = circ_ang2rad(deg); % convert to radians
mean_rad = circ_mean(rad);
circ_rad2ang(mean_rad)
% circ_rad2ang(2*pi+mean_rad) % mean degree, if the mean angle is negative,
% then to transform it to positive, use (2*pi+mean_rad)

figure(3)
%circ_plot(rad,'pretty','bo',true,'linewidth',2,'color','r')
circ_plot(rad,'hist',[],20,true,true,'linewidth',2,'color','r')
% Rayleigh test
p_uniform = circ_rtest(rad)

%%
% individual results
sub_all = [1:19,24:36,38,40,41,43:46];
sub_valid = setdiff(sub_all,outlier);
[deg,sub_valid']
%%
% % V test
% p_180 = circ_vtest(rad,circ_ang2rad(180));
% p_0 = circ_vtest(rad,circ_ang2rad(0));

% compute confidence intervals (recommended to be used to test if the
% hypothesized orientation is signficant)
interval=circ_confmean(rad,0.05)
fprintf('Mean, low/up 95 perc. CI:\t\t\t%.2f \t%.2f\n', circ_rad2ang([mean_rad-interval mean_rad+interval]))

% One sample test for the mean angle
p = circ_mtest(rad,circ_ang2rad(180)); % expect p=1 
p = circ_mtest(rad,circ_ang2rad(0)); % expect p=0
% the requirements for confidence level methods is the data
% should be uniformaly distributed.

% One-Sample test for the mean angle.
%     H0: the population has mean dir.
%     HA: the population has not mean dir.
%  h:     0 if H0 can not be rejected, 1 otherwise


% test = [45 -50 60 -70 -60 90 -180 -145]';
% test = [0 -90 -30 -60]';
% rad_test = circ_ang2rad(test);
% mean_rad_test = circ_mean(rad_test);
% circ_rad2ang(mean_rad_test+2*pi)
% figure(5)
% circ_plot(rad_test,'pretty','bo',true,'linewidth',2,'color','r')


%%
% %% to change frequency band, refer to this:
% fs=25;
% t=0.3:1/fs:(1.5-1/fs);
% times=3;
% nfft=(2.^(nextpow2(length(t))))*times; % times of next power of 2
% freqAxis1=(0:nfft/2-1)*fs/nfft;
% index=[1:length(freqAxis1)];
% refer=[freqAxis1',index']

% note: use different smoothing parameters will not affect the circular
% results (doubtable)