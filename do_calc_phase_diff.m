
function AllSubPhase = do_calc_phase_diff(outlier,con,times,smooth,targetHz)
%% date: 10/5/17
% this script can be used to calcualte phase difference between two
% conditions
% if con=1, left visual field; if con=2, right visual field
% outlier: excluded subject index, like [5,11,17,18]
% times: determine nfft size
% smooth: indicate smooth method: 'lp' or 'sg'
% targetHz: get target Hz from 'find_sigHz_plot.Rmd'. can be one value or a
% band
%%
% get valid subjects
sub_all = [1:19,24:36,38,40,41,43:46];
sub_valid = setdiff(sub_all,outlier);
sub_x = {};
for i=1:length(sub_valid)
    sub_x(i)=cellstr(sprintf('sub%s',num2str(sub_valid(i))));
end
sub = sub_x';
% set compare conditions
if con == 1 && strcmp('sg',smooth)
    con_x=['SG_Smoothed_Con50L_valid  ';'SG_Smoothed_Con50L_invalid'];
elseif con == 1 && strcmp('lp',smooth)
    con_x=['LP_Smoothed_Con50L_valid  ';'LP_Smoothed_Con50L_invalid'];
elseif con == 2 && strcmp('sg',smooth)
    con_x=['SG_Smoothed_Con50R_valid  ';'SG_Smoothed_Con50R_invalid'];
elseif con == 2 && strcmp('lp',smooth)
    con_x=['LP_Smoothed_Con50R_valid  ';'LP_Smoothed_Con50R_invalid'];
end

con_x=cellstr(con_x);

AllSubPhase=[];

for s=1:length(sub)
    valid = sprintf('%s_%s.txt',char(con_x(1)),char(sub(s)));
    invalid = sprintf('%s_%s.txt',char(con_x(2)),char(sub(s)));

    % do spectral analysis on observed valid data
    data_val = dlmread(valid);
    data_inval = dlmread(invalid);
    error_val = data_val(:,2); % only error data will be used % verticle data
    error_inval = data_inval(:,2);

    % set parameters    
    fs=25;
    t=0.3:1/fs:1.5;%30 points,next power of 2=96 (whether using 0.3:1/fs:(1.5-1/fs) or not, the nfft will be the same)
    % but since the freqAxis1 included the nyquist frequency, to match the
    % size, here we should use 1.5, instead of 1.5-1/fs
    nfft=(2.^(nextpow2(length(t))))*times; % times of next power of 2
    zeronum=nfft-length(t); % to make sure the nfft is the result of (next power of 2)*times

    dc_val=detrendnonlin(error_val); 
    dc_inval=detrendnonlin(error_inval); % detrend to get rid of the peak frequency at 0Hz

    hann_val=dc_val'.*hanning(length(dc_val))'; % hanning window
    zero_hann_val=[hann_val, zeros(1,zeronum)]; % zero-padding

    hann_inval=dc_inval'.*hanning(length(dc_inval))'; % hanning window
    zero_hann_inval=[hann_inval, zeros(1,zeronum)]; % zero-padding

    val_fft=fft(zero_hann_val,nfft);
    inval_fft=fft(zero_hann_inval,nfft);
    freqAxis1=(0:nfft/2)*fs/nfft; % it is fine to include the nyquist frequency, because this frquency
                                  % is not of interest, only need to make
                                  % sure the valid and invalid are matched.


    phase_diff=zeros(1,nfft/2+1);
    for i=1:length(freqAxis1)
        phase_diff(i)=angle(conj(val_fft(i))*inval_fft(i)); % key part
        phase_diff(i)=phase_diff(i)*180/pi;
    end

    format short
    % extract the results
    index=[1:length(freqAxis1)]';
    conj_multi=[freqAxis1',phase_diff',index];
    
    % find the the phase difference at target Hz
    if length(targetHz) == 1
        [cc iindex]=min(abs(freqAxis1-targetHz));
        fprintf('The target Hz is %.4fHz\n',freqAxis1(iindex));
        freqband=conj_multi(iindex,1:2); % get phase difference at target Hz
    elseif length(targetHz) == 2
        [cc1 iindex1]=min(abs(freqAxis1-targetHz(1)));
        [cc2 iindex2]=min(abs(freqAxis1-targetHz(2)));
        fprintf('The target Hz is from %.4fHz to %.4fHz\n',freqAxis1(iindex1),freqAxis1(iindex2));
        freqband=conj_multi(iindex1:iindex2,1:2); % get phase difference at target Hz        
    end

    MeanPhase=circ_rad2ang(circ_mean(circ_ang2rad(freqband(:,2)))); % calculate the average phase difference in the band
    AllSubPhase=[AllSubPhase,MeanPhase]; % store each phase difference 

end

% % read out results
% if con == 1 
%     fid=fopen('Con50L_peakF_PhaseDiff.txt','w');
% elseif con == 2
%     fid=fopen('Con50R_peakF_PhaseDiff.txt','w');
% end
% fprintf(fid,'%f\n',AllSubPhase);
% fclose(fid);

end

% %% testing
% % signal parameters
% fs = 25;
% f0 = 3;
% f1 = 5;
% T = 1.1;
% 
% % preparation of the time vector
% N = round(T*fs);
% t = (0:N-1)/fs;
% 
% % generation of the signal
% x = sin(2*pi*f0*t) + sin(2*pi*f1*t+pi/3) + 0.02*randn(1, N);
% y = 0.5*sign(sin(2*pi*f0*t + pi/4)) + sin(2*pi*f1*t+pi/6) + 0.02*randn(1, N);
%         
% nfft=(2.^(nextpow2(length(t))))*3; % times of next power of 2
% zeronum=nfft-length(t); % to make sure the nfft is the result of (next power of 2)*times
% 
% dc_val=detrendnonlin(x); 
% dc_inval=detrendnonlin(y); % detrend to get rid of the peak frequency at 0Hz
% 
% hann_val=dc_val'.*hanning(length(dc_val))'; % hanning window
% zero_hann_val=[hann_val, zeros(1,zeronum)]; % zero-padding
% 
% hann_inval=dc_inval'.*hanning(length(dc_inval))'; % hanning window
% zero_hann_inval=[hann_inval, zeros(1,zeronum)]; % zero-padding
% 
% val_fft=fft(zero_hann_val,nfft);
% inval_fft=fft(zero_hann_inval,nfft);
% freqAxis1=(0:nfft/2)*fs/nfft; % it is fine to include the nyquist frequency, because this frquency
%                               % is not of interest, only need to make
%                               % sure the valid and invalid are matched.
% 
% 
% phase_diff=zeros(1,nfft/2+1);
% for i=1:length(freqAxis1)
%     phase_diff(i)=angle(conj(val_fft(i))*inval_fft(i)); % key part
%     phase_diff(i)=phase_diff(i)*180/pi;
% end
% 
% 
% index=[1:length(freqAxis1)]';
% conj_multi=[freqAxis1',phase_diff',index];
