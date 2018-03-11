
%% Update date: Aug 17,2016
% update date: Aug 18, 2016
% some changes are made to adapt to new experiment:
% 1: 60-40 condition is removed
% 2: probe.error90 and probe.error180 are added

%% This script was written to read out necessary data to a text file
% note: This script is relied on MGL, so it can only be run on a Mac
% system.
% You need to provide -directory : e.g. './Sub1_S2'
%                     -subj      : e.g. 6
%                     -session   : e.g. 1
%                     -outputdir : e.g. './Sub1_S2
% Example: One_ExtractVariables('../Gabor_WM_v3_R/data/Sub3_S1',3,1,'../Gabor_WM_v3_R/data/Sub3_S1')
%%
function One_ExtractVariables(directory,subj,session,outputdir)
    SubjNo=subj;
    SessNum=session;
    dataname=sprintf('Sub%i_S%i',SubjNo,SessNum);
    dataname=[]; %each folder will have an individual name in this case: to overcome overwrite?

    %%%%%%%%%%%%%%%%%%%%%% Main Program %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    myfiles=dir(fullfile(directory,'*.mat'));
    blockNum=0;
    for k=1:length(myfiles)
        filename=myfiles(k).name;
        fullfilename=fullfile(directory,filename);
        session=load(fullfilename);

        exp=getTaskParameters(session.myscreen,session.task);
        if strcmp(session.myscreen.type,'100-0') 
            Stimutype = 1;                       
            ISIs=session.stimulus.allSOAs;
            ntrial=exp.nTrials;
        elseif strcmp(session.myscreen.type,'50-50')
            Stimutype = 2;
            ISIs=session.stimulus.allSOAs;
            ntrial=exp.nTrials;
        end

        ISI_first7 = [];
        for ii = 1:ntrial % ii is trial number
            for jj = 1:7 %jj=segment number -1
                ISI_first7(ii,jj) = exp.trials(ii).segtime(jj+1) - exp.trials(ii).segtime(jj);
            end
        end
        RT=ISI_first7(:,6); % response time % not equals to exp.responseTime

        % extract error data
        error90=session.stimulus.probe.error90;
        error180=session.stimulus.probe.error180;
        % extract ISI data
        time=ISIs(exp.parameter.soaIdx);
        te=[time',error90',error180']; %combine time and error
        % extract condition informaiton
        if strcmp(session.myscreen.type,'100-0')
            val=repmat(-1,exp.nTrials,1); % assign an unexisted -1 to 100-0 condition
        else
            val=(exp.parameter.valid)';
        end

        loc=(exp.parameter.retroCuedLoc)';

        % subject information
        sub=repmat(SubjNo,exp.nTrials,1);
        sess=repmat(SessNum,exp.nTrials,1);
        blockNum=blockNum+1;
        block=repmat(blockNum,exp.nTrials,1);
        type=repmat(Stimutype,exp.nTrials,1);
        response=(session.stimulus.probe.setting)';
        orileft=(session.stimulus.oriLeft.adjust)';
        oriright=(session.stimulus.oriRight.adjust)';
        % all informaiton
        BLOCK{blockNum}=[sub,sess,block,type,loc,val,response,orileft,oriright,RT,te];
        dataname=[dataname;BLOCK{blockNum}]; 
    end

    % use dlmwrite to write matrix to txt file.
    writedata=sprintf('Sub%i_Session%i_rawinfo.txt',SubjNo,SessNum);
    outputwrite=sprintf('%s/%s',outputdir,writedata);
    dlmwrite(outputwrite,dataname,'delimiter','\t')
end






