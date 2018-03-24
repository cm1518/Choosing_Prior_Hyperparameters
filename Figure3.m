%% This file replicates figure 3 of the paper. 
% -------------------------------------------------------------------------
% Figure 3:     Comparison of prior vs posterior distribution of 
%               hyperparameters  while comparing to the fixed values chosen 
%               in Primiceri (2005,REStud). The first row of subplots will 
%               report results for the Euroarea data while the second row 
%               reports results for the UK data.
% -------------------------------------------------------------------------
clc;

do_stability    = 1;

if ispc
    str1 = [ main_path '\Results\Euroarea'];
    str2 = [ main_path '\Results\UK'];    
else
    str1 = [ main_path '/Results/Euroarea'];
    str2 = [ main_path '/Results/UK'];
end

if ~do_stability
    cd(str1);
    try RWMH.Kappa1.Prior(3).Data(1).Posterior=load('FinalMat_Euroarea_Prior_3_SamplingScheme_0kappas_posterior'); end

    cd(str2);
    try RWMH.Kappa1.Prior(3).Data(2).Posterior=load('FinalMat_UK_Prior_3_SamplingScheme_0kappas_posterior'); end
else
    cd(str1);
    try RWMH.Kappa1.Prior(3).Data(1).Posterior=load('FinalMat_Euroarea_Prior_3_SamplingScheme_0_StableDrawskappas_posterior'); end
    
    cd(str2);
    try RWMH.Kappa1.Prior(3).Data(2).Posterior=load('FinalMat_UK_Prior_3_SamplingScheme_0_StableDrawskappas_posterior'); end                                                     
end

% =========================================================================
ytit    = {'Euroarea','UK'};
tit     = {'\kappa_{\Omega_{b,constant}}','\kappa_{\Omega_{b,dynamic}}','\kappa_{\Omega_{h}}','\kappa_{\Omega_{a}}'};

color_tmp2  = [0 0 0]+.8;
color_tmp1  = [0 0 0]+.1;
color_tmp3  = [.4 .4 .4];
color_tmp4  = [.4 .4 .4];

fix_hyp     = [.01 .01 .01 .1];

LineWidth   = 1.5;
LineWidth2  = 1.0;
LineColor   = [0 0 0];
LineStyle   = '-';
FontSize    = 12;

beg_draw    = options_.nrep/2;% 
DO_SAVE     = 1;

% =========================================================================
% Specification: RW-MH, MultiKappa
% =========================================================================

df2     = 2;
scale3  = .1; 
xxx     = linspace(sqrt(eps),1,500)';
tmp     = exp(linv_gam_pdf(xxx,df2/2,scale3*df2/2));

xmin =[];
xmax =[];

xmin_mat= nan(2,4);
xmax_mat= nan(2,4);
ymax_mat= nan(2,4);
max_y   = 100;

figure;orient('Landscape');
plt_index = 0; gg = []; xmin=nan(2,4); xmax=nan(2,4);
for ii = 1 : 2
    if ii == 1 % Euroarea;
        try
            tmp_data = RWMH.Kappa1.Prior(3).Data(ii).Posterior;
            succes = 1;
        catch
            succes = 0;
        end
    elseif ii == 2 % UK;
        try
            tmp_data = RWMH.Kappa1.Prior(3).Data(ii).Posterior;
            succes = 1;
        catch
            succes = 0;
        end
    end
    
    if succes
        gg(plt_index+1)=subplot(2,4,plt_index+1); hold('on');
        [ff,xx] = ksdensity(tmp_data.k_Q_post(beg_draw:end,1)); 
        area(xx,ff,'FaceColor',color_tmp1,'EdgeColor','None'); alpha(.5)
        x_hndl = plot(xxx,tmp,'Color',color_tmp2,'LineStyle',LineStyle,'LineWidth',LineWidth+.5);
        if max(xx) < .1,x = [min(xx) .1]; y = max(ff);else,x = [min(xx) max(xx)];y = max(ff);end
        if any(plt_index==[0:4:15]);ylabel(sprintf('%s',char(ytit(ii))),'FontSize',FontSize,'FontWeight','Bold');end
        if plt_index<4;title(sprintf('%s',char(tit(plt_index+1))),'FontSize',FontSize,'FontWeight','Bold');end
        xmin_mat(ii,1) = min([x(:);fix_hyp(1)]); xmax_mat(ii,1) = max([x(:);fix_hyp(1)]); ymax_mat(ii,1) = y;
        max_y = max([ff(:);tmp(:)]);
        hh = line([fix_hyp(1);fix_hyp(1)],[0;max_y]);set(hh,'Color',LineColor);set(hh,'LineWidth',LineWidth2,'LineStyle',LineStyle)

        axis tight
        
        gg(plt_index+2)=subplot(2,4,plt_index+2);hold('on');
        [ff,xx] = ksdensity(tmp_data.k_Q_post(beg_draw:end,2)); 
        area(xx,ff,'FaceColor',color_tmp1,'EdgeColor','None'); alpha(.5)
        x_hndl = plot(xxx,tmp,'Color',color_tmp2,'LineStyle',LineStyle,'LineWidth',LineWidth+.5);
        if max(xx) < .1,x = [min(xx) .1]; y = max(ff);else,x = [min(xx) max(xx)];y = max(ff);end
        if plt_index<4;title(sprintf('%s',char(tit(plt_index+2))),'FontSize',FontSize,'FontWeight','Bold');end
        xmin_mat(ii,2) = min([x(:);fix_hyp(2)]); xmax_mat(ii,2) = max([x(:);fix_hyp(2)]); ymax_mat(ii,2) = y;
        max_y = max([ff(:);tmp(:)]);
        hh = line([fix_hyp(2);fix_hyp(2)],[0;max_y]);set(hh,'Color',LineColor);set(hh,'LineWidth',LineWidth2,'LineStyle',LineStyle)
        axis tight
        
        gg(plt_index+3)=subplot(2,4,plt_index+3);hold('on');
        [ff,xx] = ksdensity(tmp_data.k_W_post); 
        area(xx,ff,'FaceColor',color_tmp1,'EdgeColor','None'); alpha(.5)
        x_hndl = plot(xxx,tmp,'Color',color_tmp2,'LineStyle',LineStyle,'LineWidth',LineWidth+.5);
        if max(xx) < .1,x = [min(xx) .1]; y = max(ff);else,x = [min(xx) max(xx)];y = max(ff);end
        if plt_index<4;title(sprintf('%s',char(tit(plt_index+3))),'FontSize',FontSize,'FontWeight','Bold');end
        xmin_mat(ii,3) = min([x(:);fix_hyp(3)]); xmax_mat(ii,3) = max([x(:);fix_hyp(3)]); ymax_mat(ii,3) = y;
        max_y = max([ff(:);tmp(:)]);
        hh = line([fix_hyp(2);fix_hyp(2)],[0;max_y]);set(hh,'Color',LineColor);set(hh,'LineWidth',LineWidth2,'LineStyle',LineStyle)
        axis tight
        
        gg(plt_index+4)=subplot(2,4,plt_index+4);hold('on');
        [ff,xx] = ksdensity(tmp_data.k_S_post); 
        area(xx,ff,'FaceColor',color_tmp1,'EdgeColor','None'); alpha(.5)
        x_hndl = plot(xxx,tmp,'Color',color_tmp2,'LineStyle',LineStyle,'LineWidth',LineWidth+.5);
        if max(xx) < .1,x = [min(xx) .1]; y = max(ff);else,x = [min(xx) max(xx)];y = max(ff);end
        if plt_index<4;title(sprintf('%s',char(tit(plt_index+4))),'FontSize',FontSize,'FontWeight','Bold');end
        xmin_mat(ii,4) = min([x(:);fix_hyp(4)]); xmax_mat(ii,4) = max([x(:);fix_hyp(4)]); ymax_mat(ii,4) = y;
        max_y = max([ff(:);tmp(:)]);
        hh = line([fix_hyp(4);fix_hyp(4)],[0;max_y]);set(hh,'Color',LineColor);set(hh,'LineWidth',LineWidth2,'LineStyle',LineStyle)
        axis tight
        
        xmin(ii,1) = min([min(x(:));xmin(ii,1)]);
        xmin(ii,2) = min([min(x(:));xmin(ii,2)]);
        xmin(ii,3) = min([min(x(:));xmin(ii,3)]);
        xmin(ii,4) = min([min(x(:));xmin(ii,4)]);
        xmax(ii,1) = max([max(x(:));xmax(ii,1)]);
        xmax(ii,2) = max([max(x(:));xmax(ii,2)]);
        xmax(ii,3) = max([max(x(:));xmax(ii,3)]);
        xmax(ii,4) = max([max(x(:));xmax(ii,4)]);
        
    end
    plt_index=plt_index+4;
end

set(gg,'Xgrid','on','Ygrid','on','box','on')
set(gg,'FontSize',FontSize-6);
set(gg,'LineWidth',.5);
for ii = 1 : length(gg)
    set(get(gg(ii),'Title'),'FontSize',14,'FontWeight','Normal')
end

cd ..

for ii = 1 : 4
    set(gg([ii:4:8]),'XLim',[0 max(xmax_mat(:,ii)) ])
end

sub_index = 0;
for ii = 1 : 2
    for jj =  1 : 4
        sub_index = sub_index + 1;
        set(gg(sub_index),'YLim',[0 ymax_mat(ii,jj)])
        
    end
end

if DO_SAVE
    if ~do_stability
        print(gcf, '-depsc', '-r290',sprintf('EuroareaUKMultiKappaRWMH_RevisionAxisAdjustedBW.eps'))
    else
        print(gcf, '-depsc', '-r290',sprintf('EuroareaUKMultiKappaRWMH_StableDraws_RevisionAxisAdjustedBW.eps'))
    end

end
close('all')

cd ..
