function [ k_W,accept_flag ] = draw_k_W(W,k_Wold, W_prmean,W_prvar,tuning_para,options_,M,p);
% =========================================================================
%   Input: 
%       Qdraw       : old draw iW of W
%       W_prmean    : iW Mean given old draw of W
%       W_prvar     : iW Variance given old draw of W
%       k_Wmax      : max value tuning parameter can take (discipline on degree of time variation)
% 
%   Output: 
%       k_W         : new draw of degree of time variation W
% =========================================================================

%% Metropolis step for k_W
% -------------------------------------------------------------------------
% Get likelihood of W given the current Kw. This is just the prior of W(given Kw)

[lik_old,loglik_old]    = invwishpdf(W,W_prmean,W_prvar);
lpriordens_old          = priordens(k_Wold, options_);
loglik_old              = loglik_old+lpriordens_old; % posterior kernel of the old draw

% Draw a new Kw
if options_.sampling==0         % random walk MH
    k_Wnew = k_Wold+randn(1,1)*tuning_para; 
elseif options_.sampling==1     % independent MH with uniform proposal
    k_Wnew = rand(1,1)*tuning_para;
elseif options_.sampling==2     % independent MH with half-cauchy proposal
    k_Wnew = rhalfcauchy(1,tuning_para);
end

% decide whether to keep the new draw based on prior support
if k_Wnew<options_.lower_bound || k_Wnew>options_.upper_bound
    accept_flag = 0;
    k_W         = k_Wold;
    return
end

% Calculate the likelihood of W given the new Kw
W_prmean_new = (k_Wnew/k_Wold)^2*W_prmean; % k_W changes the mean, so

% divide out the old value and multiply in the new.
[lik_new,loglik_new] = invwishpdf(W,W_prmean_new,W_prvar);

lpriordens_new = priordens(k_Wnew, options_); %evaluate prior density with the new kW draw  

loglik_new=loglik_new+lpriordens_new;
% Calculate the acceptance ratio; Since the prior on k_W is uniform, I have omitted it
logalpha = loglik_new-loglik_old;

% if using independent MH with half-cauchy proposal, need to add proposal densities. 
% For RW-MH and independent MH with uniform proposal, there is no need for this step
if options_.sampling==2 
    logalpha=logalpha+dhalfcauchy(k_Wold,tuning_para)- dhalfcauchy(k_Wnew,tuning_para); % adding proposal densities for independent MH (half-cauchy)...
end;

% Accept according to Metropolis algorithm
if(log(rand)<logalpha) 
    k_W = k_Wnew;
    accept_flag=1;
else
    k_W = k_Wold;
    accept_flag=0;
end



