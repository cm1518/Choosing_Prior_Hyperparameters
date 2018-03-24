%% DoConditionalMoments_JBES.m
% -------------------------------------------------------------------------
% This file calculates evolving long-run means (conditional on time t) for 
% the TVP-VAR with estimated hyperparameters. Local stability conditions
% areimposed for the long-run mean to be defined.
% -------------------------------------------------------------------------
clear('all'); clc; close('all')

do_stability = 1;

WhichVersion = 5;   % 5: Euroarea
                    % 9: UK
                    
prior_sel   = 1;    % 1: fixed
                    % 2: uniform
                    % 3: inv-gamma
                    % 4: half-cauchy
                    % 5: half-t

sampling_sel= 0;    % 0: RW-MH, (Benchmark choice!)
                    % 1: I-MH with uniform proposal 
                    % 2: I-MH w/ half-cauchy proposal

data_sel    = WhichVersion; 


main_path = pwd;
if ispc
    data_path = [main_path '/Data/'];
else
    data_path = [main_path '/Data/'];
end
cd(main_path)

switch WhichVersion
    case 5
        savefile        = 'Euroarea';
        if ispc
            save_path = [main_path '\Results\' savefile '\'];
        else
            save_path = [main_path '/Results/' savefile '/'];
        end
        cd(save_path)
        eval(sprintf('load %sFinalMat_%s_Prior_%d_SamplingScheme_%d_StableDraws',save_path,savefile,prior_sel,sampling_sel))
        load([data_path 'euroarea_data'])
        Y               = data;
        options_.tau    = 40;
        
    case 9
        savefile        = 'UK';
        if ispc
            save_path = [main_path '\Results\' savefile '\'];
        else
            save_path = [main_path '/Results/' savefile '/'];
        end
        cd(save_path)
        eval(sprintf('load %sFinalMat_%s_Prior_%d_SamplingScheme_%d_StableDraws',save_path,savefile,prior_sel,sampling_sel))
        load([data_path 'uk_data'])
        Y               = data;
        options_.tau    = 24;
end


number_of_draws = size(Bt_post,3);  % number of draws stored

M               = size(Sigt_post,1);
t               = size(Sigt_post,2);
p               = options_.p;

%% Priliminaries
tic;
counter         = zeros(t,number_of_draws);
mean_tmp        = nan(M,t,number_of_draws);
counter_i       = 0;

for j = 1 : number_of_draws
    if mod(j,100)==0,
        fprintf('Draw:\t%7.0f\n',j)
    end
    Btdraw=Bt_post(:,:,j);
    
    biga = zeros(M*p,M*p);
    for jj = 1:p-1
        biga(jj*M+1:M*(jj+1),M*(jj-1)+1:jj*M) = eye(M);
    end
    
    
    for i=1:t
    
        bbtemp = Btdraw(M+1:end,i);  % get the draw of B(t) at time i=1,...,T  (exclude intercept)

        splace = 0;
        for ii = 1:p
            for iii = 1:M
                biga(iii,(ii-1)*M+1:ii*M) = bbtemp(splace+1:splace+M,1)';
                splace = splace + M;
            end
        end

        ctemp1 = biga;
        if max(abs(eig(ctemp1)))>=1;
            counter(i,j)=counter(i,j)+1;
        end
    end
    if ~sum(counter(:,j)) % local stability conditions satisfied
        counter_i = counter_i + 1;
        for i = 1 : t;
            tmp=[];
            bbtemp = Btdraw(M+1:end,i);  % get the draw of B(t) at time i=1,...,T  (exclude intercept)
            splace = 0;
            for ii = 1:p
                for iii = 1:M
                    biga(iii,(ii-1)*M+1:ii*M) = bbtemp(splace+1:splace+M,1)';
                    splace = splace + M;
                end
            end
            intercepts  = Btdraw(1:M,i);
            ctemp1 = biga;
            tmp = [eye(size(ctemp1,1)) - ctemp1] \ [intercepts;zeros(M,1)];
            mean_tmp(:,i,counter_i) = tmp(1:M,1);
        end
    end
end;

mean_tmp = mean_tmp(:,:,1:counter_i);
var_name = {'Real GDP growth','Inflation','Short Rate'};

clear -regexp _post$;

save_str = sprintf('FinalConmeansMat_JBES%s_Prior_%d_SamplingScheme_%d',savefile,prior_sel,sampling_sel);

eval(sprintf('save %s.mat',save_str))

cd ../..
