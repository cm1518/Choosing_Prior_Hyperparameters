%% This file replicates figure 4 and 5 of the paper. 
% -------------------------------------------------------------------------
% Figure 4/5:   Comparison of posterior IRF for based on fixed vs. estimated 
%               hyperparameters for the country of choice [in ll. 16-18] at 
%               selected points in time.
%
% -------------------------------------------------------------------------
close('all'); clc

do_stability    = 1;
sampling_sel    = 0;    % 0: RW-MH, (Benchmark choice!)
                        % 1: I-MH with uniform proposal
                        % 2: I-MH w/ half-cauchy proposal
prior_sel       = 3;

% savefile = 'Euroarea';
savefile = 'UK';

save_path = [main_path '/Results/' savefile '/'];
cd(save_path)

if ~do_stability % UNSTABLE CASE / Forecast means
    eval(sprintf('load FinalIRFMat_%s_Prior_%d_SamplingScheme_%d',savefile,1,sampling_sel))
    mean_tmp_fix = tube.recnirf;clear tube
    eval(sprintf('load FinalIRFMat_%s_Prior_%d_SamplingScheme_%d',savefile,3,sampling_sel))
else % STABLE CASE / Conditional means
    eval(sprintf('load FinalIRFMat_%s_Prior_%d_SamplingScheme_%d_StableDraws',savefile,1,sampling_sel))
    mean_tmp_fix = tube.recnirf;clear tube
    eval(sprintf('load FinalIRFMat_%s_Prior_%d_SamplingScheme_%d_StableDraws',savefile,3,sampling_sel))
end

switch savefile
    case 'Euroarea'
        load([data_path 'euroarea_data'])
    case 'UK'
        load([data_path 'uk_data'])
end

yearlab     = yearlab(options_.tau+p+1:end);
var_name    = {'Real GDP growth','Inflation','Short Rate'};
varname     = var_name;

%% 1. Plot settings
year_indx   = [find(yearlab==1987),find(yearlab==1990),find(yearlab==1995),...
    find(yearlab==2000),find(yearlab==2005),find(yearlab==2010)];
year_indx2  = [find(yearlab==1987),find(yearlab==1990),find(yearlab==2000),find(yearlab==2010)];
RGB         = [1 1 1];
DO_SAVE     = 1;
EdgeColor1  = [0 0 0]+.2;  PatchColor1 = [0 0 0]+.2;
EdgeColor2  = [0 0 0]+.6;  PatchColor2 = [0 0 0]+.7;
xmin        = [];
xmax        = [];
ylim_tmp    = [];
clear gg
horz_end = h;
for nvar = 1; %[1 2 3];
    ylim_tmp    = [];
    y_lim = [];
    figure('Name',sprintf('%s Time IRF to 1Std MP Shock',char(varname(nvar))));orient('Landscape')
    for ii = 1 : length(year_indx2)
        gg(ii)=subplot(2,2,ii);
        hold('on');grid('on');
        % Zero line
        plot(1:horz_end,zeros(1,horz_end),'Color',[0 0 0],'LineStyle','-','LineWidth',.5)
%         % fixed hyperparameters
        plot(1:1:horz_end,squeeze(mean_tmp_fix(nvar,year_indx2(ii),1:horz_end,[3 17])), 'Color',EdgeColor1,'LineStyle','-','LineWidth',1)
        %%%plot(1:1:horz_end,squeeze(mean_tmp_fix(nvar,year_indx2(ii),1:horz_end,10)),     'Color',EdgeColor1,'LineStyle','-','LineWidth',2)
        % Estimated hyperparameters
        for jj = 3%3:7 % 1:9
            tmp = [1:1:horz_end horz_end:-1:1;squeeze(tube.recnirf(nvar,year_indx2(ii),1:horz_end,jj))' flipdim(squeeze(tube.recnirf(nvar,year_indx2(ii),1:horz_end,end+1-jj)),1)']';
            patch(tmp(:,1),tmp(:,2),PatchColor2,'EdgeColor',EdgeColor2)
        end
        %%%plot(1:1:horz_end,squeeze(tube.recnirf(nvar,year_indx2(ii),1:horz_end,10)),'LineStyle','-','Color',EdgeColor2,'LineWidth',2)
        
        % if ii < 3
        %     title(sprintf('%s IRF in %d:Q4',char(varname(nvar)),yearlab(year_indx2(ii))-1),'FontSize',12);%,'FontWeight','Bold')
        % end
        title(sprintf('%s IRF in %d:Q4',char(varname(nvar)),yearlab(year_indx2(ii))-1),'FontSize',12);%,'FontWeight','Bold')
        
        axis('tight')
        set(gca,'XTick',[1 5:4:21],'XTicklabel',[0 1:5])
        if any(ii==[1:2])
            set(gca,'XTicklabel',[]),
        else
            xlabel('Horizon (Years)')
        end
        if any(ii==[1 3])
            ylabel(sprintf('Percent'))
            % ylabel(sprintf('%s Reponse (Percent)',char(varname(nvar))))
        end
        alpha(.5);box('on'),hold('off')
        set(findall(gcf, '-property', 'FontSize'), 'FontSize', 12, 'fontWeight', 'normal','FontName','Times New Roman')
        ylim_tmp = [ min([ylim_tmp ylim]') max([ylim_tmp ylim]')];
    end
    % set(gg,'XLim',[1 horz_end],'FontSize',12)
    set(gg,'YLim',[min(ylim_tmp(:)) max(ylim_tmp(:))],'XLim',[1 horz_end],'FontSize',12)
    % set(gg,'YLim',[floor(min(ylim_tmp(:))/10)*10 ceil(max(ylim_tmp(:))/10)*10],'XLim',[1 horz_end],'FontSize',12)
    if DO_SAVE
        if ~do_stability % UNSTABLE CASE / Forecast means
            print(gcf, '-depsc', '-r300',sprintf('TimeSelectedIRF_%s_PriorComp_SamplingScheme_%d_%s_BW.eps',savefile,sampling_sel,char(varname(nvar))))
        else
            print(gcf, '-depsc', '-r300',sprintf('TimeSelectedIRF_%s_PriorComp_SamplingScheme_%d_StableDraws_%s_BW.eps',savefile,sampling_sel,char(varname(nvar))))
        end
    end
end
close('all')

cd ../..
