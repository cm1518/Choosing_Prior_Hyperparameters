%% This file replicates figure 6 of the paper. 
% -------------------------------------------------------------------------
% Figure 6:     Comparison of posterior median of evolving reduced form
%               parametres based on fixed vs. estimated hyperparameters for 
%               the country of choice [in ll. 15-16] over time.
%
% -------------------------------------------------------------------------
close('all'); clc

EdgeColor1  = [0 0 0]+.2;  PatchColor1 = [0 0 0]+.2;
EdgeColor2  = [0 0 0]+.6;  PatchColor2 = [0 0 0]+.7;


do_stability    = 1;
WhichVersion    = 5;    % 5: Euroarea
                        % 9: UK
data_sel        = WhichVersion;

prior_sel       = 1;    % 1: fixed
                        % 2: uniform
                        % 3: inv-gamma
                        % 4: half-cauchy
                        % 5: half-t

spl_index       = [1:3:9;2:3:9;3:3:9]';

sampling_sel    = 0;    % 0: RW-MH, (Benchmark choice!)
                        % 1: I-MH with uniform proposal
                        % 2: I-MH w/ half-cauchy proposal

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
        if ~do_stability
            eval(sprintf('load %sFinalMat_%s_Prior_%d_SamplingScheme_%d',save_path,savefile,1,sampling_sel))
            At_post_fix = At_post;
            Bt_post_fix = Bt_post;
            sig_post_fix = sig_post;
            eval(sprintf('load %sFinalMat_%s_Prior_%d_SamplingScheme_%d',save_path,savefile,3,sampling_sel))
        else
            eval(sprintf('load %sFinalMat_%s_Prior_%d_SamplingScheme_%d_StableDraws',save_path,savefile,1,sampling_sel))
            At_post_fix = At_post;
            Bt_post_fix = Bt_post;
            sig_post_fix = sig_post;
            eval(sprintf('load %sFinalMat_%s_Prior_%d_SamplingScheme_%d_StableDraws',save_path,savefile,3,sampling_sel))
        end
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
        if ~do_stability
            eval(sprintf('load %sFinalMat_%s_Prior_%d_SamplingScheme_%d',save_path,savefile,1,sampling_sel))
            At_post_fix = At_post;
            Bt_post_fix = Bt_post;
            sig_post_fix = sig_post;
            eval(sprintf('load %sFinalMat_%s_Prior_%d_SamplingScheme_%d',save_path,savefile,3,sampling_sel))
        else
            eval(sprintf('load %sFinalMat_%s_Prior_%d_SamplingScheme_%d_StableDraws',save_path,savefile,1,sampling_sel))
            At_post_fix = At_post;
            Bt_post_fix = Bt_post;
            sig_post_fix = sig_post;
            eval(sprintf('load %sFinalMat_%s_Prior_%d_SamplingScheme_%d_StableDraws',save_path,savefile,3,sampling_sel))
        end
        load([data_path 'uk_data'])
        Y               = data;
        options_.tau    = 24;
end

number_of_draws = size(Bt_post,3);  % number of draws stored

M               = size(Sigt_post,1);
t               = size(Sigt_post,2);
p               = options_.p;
tau             = options_.tau;

var_name = {'Real GDP growth','Inflation','Short Rate'};
yearlab     = yearlab(options_.tau+p+1:end);

figure
hold('on')
for iii = 1 : M
    subplot(3,3,spl_index(iii,1)),hold('on')
    plot(yearlab,squeeze(prctile(At_post_fix(iii,:,:),50,3)),'Color',PatchColor1,'LineWidth',2)
    plot(yearlab,squeeze(prctile(At_post(iii,:,:),50,3)),'Color',PatchColor2,'LineWidth',2)
    axis tight,box on
    if iii==1,title('a_t'),end
end
for iii = 1 : M
    subplot(3,3,spl_index(iii,2)),hold('on')
    plot(yearlab,squeeze(prctile(sig_post_fix(:,iii,:),50,3)),'Color',PatchColor1,'LineWidth',2)
    plot(yearlab,squeeze(prctile(sig_post(:,iii,:),50,3)),'Color',PatchColor2,'LineWidth',2)
    axis tight,box on
    if iii==1,title('h_t'),end
end
subplot(3,3,[spl_index(:,3)]'),hold('on')
plot(yearlab,squeeze(prctile(Bt_post_fix,50,3)),'Color',PatchColor1,'LineWidth',2)
plot(yearlab,squeeze(prctile(Bt_post,50,3)),'Color',PatchColor2,'LineWidth',2)
axis tight,box on
title('b_t')

if ~do_stability
    print(gcf, '-depsc', '-r300',sprintf('EvolvingStates_%s_Comp_SamplingScheme_%dBW.eps',savefile,sampling_sel))
else
    print(gcf, '-depsc', '-r300',sprintf('EvolvingStates_%s_Comp_SamplingScheme_%d_StableDrawsBW.eps',savefile,sampling_sel))
end
close all

cd ../..

