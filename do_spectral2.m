function do_spectral2(outlier,times,iter,parfol,windows,cutoff,order)
%% update date: 02/12/2018
% change the way to do permutation test, firstly, shuffle the raw data 1000
% times, then do the filtering and other preprocessing steps (Ronconi & Melcher et al., 2017).
% * if get rid of detrending (worse, should add back).
%% update date: 12/13/2016
% detrend the smoothed data with a second order polynomial just like
% Fiebelkorn did (2013)

%% update date: 08/31/2016
% This script is written to run randomization test on individual condition
% (no comparision between conditions), and generate output files.
%
% To use it, provide five parameters:
% "outliers": outlier index, like[5,11,17,18]
% "times": how many times you want to set for next power of 2, e.g. pow2=3
% "iter": how many iterations you want to run, e.g. iter=1000
% "parfol': which parent folder you want the generated files to be moved
% into, e.g. parfol=['perm_1000_nfft3']
% "windows": if using windows system then assigns 1, Mac assigns 0

%%
% clear all;
% clc;
% times=3;
% iter=1000;
% parfol=['Perm1000_nfft3_3MAD_notscaled_smooth'];

% iterate across multiple subjects and multiple conditions

sub_all = [1:19,24:36,38,40,41,43:46];
sub_valid = setdiff(sub_all,outlier);
sub_x = {};
for i=1:length(sub_valid)
    sub_x(i)=cellstr(sprintf('sub%s',num2str(sub_valid(i))));
end
sub = sub_x';
   
   
con_x_raw=['Con50L_valid  ';...
       'Con50L_invalid';...
       'Con50R_valid  ';...
       'Con50R_invalid';...
       'Con100L       ';...
       'Con100R       '];

% con_x=['SG_Smoothed_Con50L_valid  ';...
%        'SG_Smoothed_Con50L_invalid';...
%        'SG_Smoothed_Con50R_valid  ';...
%        'SG_Smoothed_Con50R_invalid';...
%        'SG_Smoothed_Con100L       ';...
%        'SG_Smoothed_Con100R       '];
   
con_x=['LP_Smoothed_Con50L_valid  ';...
       'LP_Smoothed_Con50L_invalid';...
       'LP_Smoothed_Con50R_valid  ';...
       'LP_Smoothed_Con50R_invalid';...
       'LP_Smoothed_Con100L       ';...
       'LP_Smoothed_Con100R       '];
   
con=cellstr(con_x);
con_raw = cellstr(con_x_raw);

%% the main program

%format long
for c = 1: length(con)
    conname=sprintf('%s_ampli_allsub',char(con(c))); % con(c) is a cell, so it's needed to be conveted to character
    conname=[]; % row: subjects; col: amplitudes of frequency bins
    damp_allsub=sprintf('%s_allsub_rand_ampli',char(con(c)));
    damp_allsub=[]; % has three dimensions: (col: amplitudes of frequency bins);
                    % (row: iterations);(third dimension: subjects)
    for s = 1: length(sub)
        filename=sprintf('%s_%s.txt',char(con(c)),char(sub(s)));
        data=dlmread(filename);
        error=data(:,2); % only error data will be used
        % do spectral analysis on raw data
        fs=25;
        t=0.3:1/fs:(1.5-1/fs);%30 points,next power of 2=96 %%%%%%!!!!!!!
        nfft=(2.^(nextpow2(length(t))))*times; % times of next power of 2
        zeronum=nfft-length(t); % to make sure the nfft is the result of (next power of 2)*times
        %error_dc=error;
        error_dc=detrendnonlin(error); % detrend with 2nd order polynomials  
        x=error_dc'.*hanning(length(error_dc))'; % hanning window
        x1=[x, zeros(1,zeronum)]; % zero-padding (before and after the signal, like Song, 2014)
        spectrum=fft(x1,nfft);
        realampli=abs(spectrum); %no normalization
       % realampli=realampli.^2; % to calculate power!!
        realampli_half=realampli(1:nfft/2); %1:nfft/2 to match the size of freq_Axis
        % freq_Axis=(0:nfft/2-1)*fs/nfft; %from 0 to nyquist frequency(not included)
      
        conname=[conname;realampli_half];               
        
        % save individual data
        % this is the OBSERVED amplitude of each condition of each subject
        save_ampli=sprintf('%s_%s_ind_ampli.txt',char(sub(s)),char(con(c)));
        fid_ampli=fopen(save_ampli,'w');
        fprintf(fid_ampli,'%f\n',realampli_half);
        fclose(fid_ampli);

        % this step is used to SHUFFLE raw data
        reps=iter; % iteration steps
        n=length(error);
        dram=zeros(n,reps); %dram stands for data randomization
                            % row: error points
                            % col: reps
        for i=1:reps
            filename_raw=sprintf('%s_%s.txt',char(con_raw(c)),char(sub(s)));
            data_raw=dlmread(filename_raw);
            error_raw=data_raw(:,2); % only error data will be used
            temp=randsample(error_raw,n); % resample without replacement
            dram(1:n,i)=temp; 
        end

        % this step is used to do spectral analysis on the SHUFFLED data for each
        % subject
        dampname=sprintf('%s_%s_rand_ampli',char(sub(s)),char(con(c)));
        dampname=zeros(reps,nfft/2); % randomization for each subject
                                     % col: amplitude for frequency bin
                                     % row: iteration
        %damp=zeros(reps,nfft/2);
        for i=1:reps
            % firstly low-passed data
            error_raw1 = dram(:,i);
            lp_filtered = filter1('lp',error_raw1,'fs',25,'fc',cutoff,'order',order); % in vertical order 
            %dram_dc = lp_filtered;
            dram_dc=detrendnonlin(lp_filtered);
            dram_hann=dram_dc'.*hanning(length(dram_dc))';
            dram_hann_zero=[dram_hann, zeros(1,zeronum)];
            dram_fft=fft(dram_hann_zero,nfft);
            dram_fft_h=dram_fft(1:nfft/2);
            dram_fft_mag=abs(dram_fft_h);
           % dram_fft_mag=dram_fft_mag.^2; % calcualte power!!
            dampname(i,:)=dram_fft_mag; % to store amplitude information for each frequency bin    
        end
        damp_allsub(:,:,s)=dampname; % compile all subjects for this condition
    end
    
        % observed data: save OBSERVED amplitude for all subject at one condition
        save_ampli_allsub_con=sprintf('%s_allsub_ampli.txt',char(con(c)));
        dlmwrite(save_ampli_allsub_con,conname,'delimiter','\t') 
        
        % observed data: average across subjects
        save_ampli_allsub_con_mean=sprintf('%s_allsub_ampli_mean.txt',char(con(c)));
        dlmwrite(save_ampli_allsub_con_mean,mean(conname),'delimiter','\t')
        
        % get amplitude for each frequency bin, repeat "reps" times,
        % average across subjects
        aver_allsub_rand_ampli=mean(damp_allsub,3); % average subjects
        save_ampli_allsub_random=sprintf('%s_mean_ampli_rand.txt',char(con(c)));
        dlmwrite(save_ampli_allsub_random,aver_allsub_rand_ampli,'delimiter','\t') 
        
        %calculate p values
        Pval=sprintf('%s_pvalue',char(con(c)));
        Pval=[];
        for b=1:size(aver_allsub_rand_ampli,2)
            refdis=aver_allsub_rand_ampli(:,b);
            observed_allbin=mean(conname);
            observed_indbin=observed_allbin(b);
            % here I used a biased method which was introduced by Smyth &
            % Phipson (2011, Permutation p-values should never be zero)
            Pval(b)=(sum(abs(refdis)>=abs(observed_indbin))+1)/(size(aver_allsub_rand_ampli,1)+1); 
        end
        save_pval=sprintf('%s_pvalue_uncor.txt',char(con(c)));
        dlmwrite(save_pval,Pval,'delimiter','\t','precision','%.8f')
     
end

if windows == 1
    random_p=sprintf('%s\\random_p',parfol); % windows
    meansub_fft_cond=sprintf('%s\\meansub_fft_cond',parfol);   
    individual_fft=sprintf('%s\\individual_fft',parfol);  
    allsub_fft_cond=sprintf('%s\\allsub_fft_cond',parfol);    
    permutated_meansub_con=sprintf('%s\\permutated_meansub_con',parfol);
else
    random_p=sprintf('%s//random_p',parfol); % mac
    meansub_fft_cond=sprintf('%s//meansub_fft_cond',parfol);
    individual_fft=sprintf('%s//individual_fft',parfol);
    allsub_fft_cond=sprintf('%s//allsub_fft_cond',parfol);
    permutated_meansub_con=sprintf('%s//permutated_meansub_con',parfol);
end
    
mkdir(random_p);
mkdir(meansub_fft_cond);
mkdir(individual_fft);
mkdir(allsub_fft_cond);
mkdir(permutated_meansub_con);
           
fprintf('--------------------- Ok,I am moving files -----------------------\n')
movefile('./*uncor.txt',sprintf('%s/random_p',parfol));
movefile('./*ind_ampli.txt',sprintf('%s/individual_fft',parfol));
movefile('./*mean.txt',sprintf('%s/meansub_fft_cond',parfol));
movefile('./*allsub_ampli.txt',sprintf('%s/allsub_fft_cond',parfol));
movefile('./*rand.txt',sprintf('%s/permutated_meansub_con',parfol));

fprintf('---------------------- I am done! Good luck with the results! ------------------------\n') 
end

