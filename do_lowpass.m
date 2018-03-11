% update date: 10/3/2017
% usage: it can be used to smooth data, use it on raw data 
% outlier: provide outlier index
% order: provide polynomial order
% framelen: provide frame length
% prefix: provide prefix of the output file, will be
% 'prefix_Smoothed_SGF_con_sub'


function do_lowpass(outlier,cutoff,order,prefix)

% iterate across multiple subjects and multiple conditions
    sub_all = [1:19,24:36,38,40,41,43:46];
    sub_valid = setdiff(sub_all,outlier);
    sub_x = {};
    for i=1:length(sub_valid)
        sub_x(i)=cellstr(sprintf('sub%s',num2str(sub_valid(i))));
    end
    sub = sub_x';
    con_x=['Con50L_valid  ';...
           'Con50L_invalid';...
           'Con50R_valid  ';...
           'Con50R_invalid';...
           'Con100L       ';...
           'Con100R       '];

    con=cellstr(con_x);

    % smooth (with low-pass filter: http://www.mathworks.com/matlabcentral/fileexchange/53534-filter1)
    for c = 1: length(con)
        for s = 1: length(sub)
            filename=sprintf('%s_%s.txt',char(con(c)),char(sub(s)));
            data=dlmread(filename);
            error=data(:,2); % only error data will be used
            smoothed=filter1('lp',error,'fs',25,'fc',cutoff,'order',order);       
            time=data(:,1);
            smooth=[time';smoothed'];
            readoutfile=sprintf('%s_Smoothed_%s_%s.txt',prefix,char(con(c)),char(sub(s)));
            fid=fopen(readoutfile,'w');
            fprintf(fid,'%f %f\n',smooth);
            fclose(fid);
        end
    end
end

