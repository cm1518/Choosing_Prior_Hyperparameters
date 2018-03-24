function [ k_S,accept_flag ] = draw_k_S(S,k_Sold, S_prmean,S_prvar,tuning_para,options_,M,p)
% =========================================================================
%   Input: 
%       Sdraw       : old draw iW of S
%       S_prmean    : iW Mean given old draw of S
%       S_prvar     : iW Variance given old draw of S
%       k_Smax      : max value tuning parameter can take (discipline on degree of time variation)
% 
%   Output: 
%       k_S         : new draw of degree of time variation S
% =========================================================================

%% Metropolis step for k_S
% -------------------------------------------------------------------------
% Get product of likelihood of S/ sum of log likelihood given the current Ks. This is just the prior of S(given Ks)

lpriordens_old = priordens(k_Sold, options_);

for i=1:M-1
    [lik_old,loglik_old_tmp(i)]    = invwishpdf(S{i},S_prmean{i},S_prvar(i));
end;

loglik_old=sum(loglik_old_tmp)+ lpriordens_old;


% Draw a new kS
if options_.sampling==0  % random walk MH
    k_Snew = k_Sold+randn(1,1)*tuning_para; 
elseif options_.sampling==1    % independent MH with uniform proposal
    k_Snew = rand(1,1)*tuning_para;
elseif options_.sampling==2    % independent MH with half-cauchy proposal
    k_Snew = rhalfcauchy(1,tuning_para);
end;

% decide whether to keep the new draw based on prior support
if k_Snew<options_.lower_bound || k_Snew>options_.upper_bound
    accept_flag=0;
    k_S = k_Sold;
    return
end


% Calculate the likelihood of W given the new Kw
% S_prmean_new = (k_Snew/k_Sold)^2*S_prmean; % k_W changes the mean, so

S_prmean_new = cell(p,1);

for i=1:M-1
    S_prmean_new{i} = (k_Snew/k_Sold)^2*S_prmean{i}; % k_S changes the mean, so
end;

% divide out the old value and multiply in the new.
% [lik_new,loglik_new] = invwishpdf(W,W_prmean_new,W_prvar);

lpriordens_new = priordens(k_Snew, options_);

for i=1:M-1
    [lik_new,loglik_new_tmp(i)]    = invwishpdf(S{i},S_prmean_new{i},S_prvar(i));
end;

loglik_new=sum(loglik_new_tmp)+ lpriordens_new;



% Calculate the acceptance ratio; Since the prior on k_S is uniform, I have omitted it
logalpha = loglik_new-loglik_old;

% if using independent MH with half-cauchy proposal, need to add proposal densities. 
% For RW-MH and independent MH with uniform proposal, there is no need for this step
if options_.sampling==2 
    logalpha=logalpha+dhalfcauchy(k_Sold,tuning_para)- dhalfcauchy(k_Snew,tuning_para); % adding proposal densities for independent MH (half-cauchy)...
end;

% Accept according to Metropolis algorithm
if(log(rand)<logalpha) && k_Snew>options_.lower_bound && k_Snew<options_.upper_bound
%if(rand<alpha) 
    k_S = k_Snew;
    accept_flag=1;
else
    k_S = k_Sold;
    accept_flag=0;
end



