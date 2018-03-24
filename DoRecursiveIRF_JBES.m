%% DoRecursiveIRF_JBES.m
% -------------------------------------------------------------------------
% This file calculates recursive impulse response functions for the 
% TVP-VAR with estimated hyperparameters. 
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
main_path   = pwd;
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
            save_path = [main_path '\Results\Euroarea\'];
        else
            save_path = [main_path '/Results/Euroarea/'];
        end
        cd(save_path)
        if ~do_stability
            eval(sprintf('load %sFinalMat_%s_Prior_%d_SamplingScheme_%d',save_path,savefile,prior_sel,sampling_sel))
        else
            eval(sprintf('load FinalMat_%s_Prior_%d_SamplingScheme_%d_StableDraws',savefile,prior_sel,sampling_sel))
        end
    case 9
        savefile        = 'UK';
        if ispc
            save_path = [main_path '\Results\UK\'];
        else
            save_path = [main_path '/Results/UK/'];
        end
        cd(save_path)
        if ~do_stability
            eval(sprintf('load FinalMat_%s_Prior_%d_SamplingScheme_%d',savefile,prior_sel,sampling_sel))
        else
            eval(sprintf('load FinalMat_%s_Prior_%d_SamplingScheme_%d_StableDraws',savefile,prior_sel,sampling_sel))
        end
end

number_of_draws = size(Bt_post,3);  % number of draws stored
M               = size(Sigt_post,1);
t               = size(Sigt_post,2);
p               = options_.p;
DO_SAVE         = 1;
nhorz           = 21;
h               = nhorz;
tube.recnirf    = nan(M,t,nhorz,number_of_draws);
irf.nirf_tmp    = [];
draw_index      = 1:number_of_draws;%randsample(nrep,number_of_draws);

tic;
for l = 1 : number_of_draws
    if mod(l,50)==0
        fprintf('Draw: %8.0f\t(%4.2f\t minutes.)\n',l,toc/60)
    end
    ctemp1  = zeros(M,M*p);
    Btdraw  = Bt_post(:,:,draw_index(l));
    %%  Create the VAR covariance matrix H(t). It holds that: A(t) x H(t) x A(t)' = SIGMA(t) x SIGMA(t) '
    capAt = zeros(M*t,M);
    for i = 1:t
        capatemp = eye(M);
        aatemp = At_post(:,i,draw_index(l));
        ic=1;
        for j = 2:M
            capatemp(j,1:j-1) = aatemp(ic:ic+j-2,1)';
            ic = ic + j - 1;
        end
        capAt((i-1)*M+1:i*M,:) = capatemp;
    end
    sigtemp = eye(M);
    sigt = zeros(M*t,M);
    for i = 1:t
        for j = 1:M
            sigtemp(j,j) = exp(0.5*Sigt_post(j,i,draw_index(l)));
        end
        sigt((i-1)*M+1:i*M,:) = sigtemp;
    end
    Ht      = zeros(M*t,M);
    Htsd    = zeros(M*t,M);
    counter = 0;
    J       = [eye(M) zeros(M,M*(p-1))];
    for i = 1 : t
        inva    = inv(capAt((i-1)*M+1:i*M,:));
        stem    = sigt((i-1)*M+1:i*M,:);
        Hsd     = capAt((i-1)*M+1:i*M,:)\stem;
        Hdraw   = Hsd*Hsd';
        A1      = reshape(Btdraw(M+1:M+M^2,i),M,M)';
        A2      = reshape(Btdraw(M+M^2+1:end,i),M,M)';
        ctemp1  = [A1 A2; eye(M*(p-1)) zeros(M*(p-1),M)];
        % =================================================================
        result.recnirf  = [];%nan(M,M,nhorz,nalpha);
        tmprec      = 0;
        tmpnrec     = 0;
        impf_tmp4   = zeros(M,M,h);
        for k = 1 : nhorz
            phi                 = J*ctemp1^(k-1)*J';
            impf_tmp4(:,:,k)    = phi*inva;             % Recursive IRF (normed to 1 %)
        end; % end horizon
        tube.recnirf(:,i,:,l)     = squeeze(impf_tmp4(:,end,:));
        
    end % end time t
end % end draw l

if DO_SAVE
    clear -regexp _post$;
    tube.recnirf      = squeeze(prctile(tube.recnirf,[5:5:95],4));
    if options_.do_stability
        save_str = sprintf('FinalIRFMat_%s_Prior_%d_SamplingScheme_%d_StableDraws',savefile,prior_sel,sampling_sel);
    else
        save_str = sprintf('FinalIRFMat_%s_Prior_%d_SamplingScheme_%d',savefile,prior_sel,sampling_sel);
    end
    eval(sprintf('save %s.mat -v7.3',save_str))
end

cd ../..

