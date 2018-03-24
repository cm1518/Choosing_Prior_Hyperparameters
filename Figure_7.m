%% This file replicates figure 4 and 5 of the paper. 
% -------------------------------------------------------------------------
% Figure 7:     Comparison of posterior long-run mean/infinite horizon 
%               forecast based on fixed vs. estimated hyperparameters for 
%               the country of choice [in ll. 15-16] over time.
%
% -------------------------------------------------------------------------
close('all'); clc

do_stability    = 1;
sampling_sel    = 0;    % 0: RW-MH, (Benchmark choice!)
                        % 1: I-MH with uniform proposal
                        % 2: I-MH w/ half-cauchy proposal
prior_sel       = 3;
% savefile       = 'Euroarea';
savefile        = 'UK';

main_path = pwd;
if ispc
    save_path = [main_path '\Results\' savefile '\'];
else
    save_path = [main_path '/Results/' savefile '/'];
end
cd(save_path)

% STABLE CASE / Conditional means
eval(sprintf('load FinalConmeansMat_JBES%s_Prior_%d_SamplingScheme_%d',savefile,1,sampling_sel))
mean_tmp_fix = mean_tmp;
eval(sprintf('load FinalConmeansMat_JBES%s_Prior_%d_SamplingScheme_%d',savefile,3,sampling_sel))


switch savefile
    case 'Euroarea'
        load([data_path 'euroarea_data'])
    case 'UK'
        load([data_path 'uk_data'])
end

var_name    = {'Real GDP growth','Inflation','Short Rate'};
yearlab     = yearlab(options_.tau+p+1:end);
EdgeColor1  = [0 0 0]+.1;  PatchColor1 = [0 0 0]+.1;
EdgeColor2  = [0 0 0]+.6;  PatchColor2 = [0 0 0]+.6;

figure
for iii = 1 : M
    subplot(1,M,iii),hold on    
    plot(yearlab,zeros(1,t),'color',[0 0 0],'LineWidth',.5)    
    tmp_time    = [yearlab flip(yearlab,2)];
    tmp_fix     = [squeeze(prctile(mean_tmp_fix(iii,:,:),16,3)) flip(squeeze(prctile(mean_tmp_fix(iii,:,:),84,3)),2)];
    tmp_kappa   = [squeeze(prctile(mean_tmp(iii,:,:),16,3))     flip(squeeze(prctile(mean_tmp(iii,:,:),84,3)),2)];
    plot(yearlab,squeeze(prctile(mean_tmp_fix(iii,:,:),[16 84 ],3)),'Color',EdgeColor1,'LineWidth',1);    
    patch(tmp_time,tmp_kappa,PatchColor2,'EdgeColor','None');
    title(var_name(iii))
    grid,axis tight,box on,alpha(.5)
end
print(gcf, '-depsc', '-r300',sprintf('LongRunMean_%s_PriorComp_SamplingScheme_%d_BW.eps',savefile,sampling_sel))
close all

cd ../..