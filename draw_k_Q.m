function [k_Q,accept_flag] = draw_k_Q(Q,k_Qold, Q_prmean,Q_prvar,tuning_para,options_,M,p)
% =========================================================================
%   Input: 
%       Q       : old draw iW of Q
%       Q_prmean    : iW Mean given old draw of Q
%       Q_prvar     : iW Variance given old draw of Q
%  
%   Output: 
%       k_Q         : new draw of degree of time variation Q
% =========================================================================
% 22.05.2015: This function implements multivariate scaling for kQ only.
% kQs are broken up into two blocks: constants and VAR parameters. Priors
% for kQs are i.i.d

% Metropolis step for K_q
% -------------------------------------------------------------------------
% Get likelihood of Q given the current kQ. This is just the prior of Q(given kQ)
[lik_old,loglik_old] = invwishpdf(Q,Q_prmean,Q_prvar);


if options_.mvs_sampling==1 % multivariate sampling of kQ?
    nkQ             = 2;
    % extracting scaling parameters k_Q_const and k_Q_coeff and evaluate prior density with the old kQ draw
    k_Qold_tmp      = [k_Qold(1,1);k_Qold(M+1,M+1)];
    lpriordens_old  = priordens(k_Qold_tmp, options_);
else
    nkQ                 = 1;
    k_Qold_tmp          = k_Qold;
    lpriordens_old(1,1) = priordens(k_Qold_tmp, options_);
end;
    
loglik_old = loglik_old + sum(lpriordens_old); % posterior kernel of the old draw

% Draw a new Kq
if options_.sampling==0                 % random walk MH
    k_Qnew = k_Qold_tmp+randn(nkQ,1)*tuning_para; 
elseif options_.sampling==1             % independent MH with uniform proposal
    k_Qnew = rand(nkQ,1)*tuning_para;
elseif options_.sampling==2             % independent MH with half-cauchy proposal
    k_Qnew = rhalfcauchy(nkQ,tuning_para);
end

% decide whether to keep the new draw based on prior support
if any(k_Qnew<options_.lower_bound) || any(k_Qnew>options_.upper_bound)
    accept_flag=0;
    k_Q = k_Qold;
    return
end


% Calculate the likelihood of Q given the new Kq
if options_.mvs_sampling==1 % multivariate sampling of kQ?
    lpriordens_new      = priordens(k_Qnew, options_); %evaluate prior density with the new kQ draw    
    k_factor            = k_Qnew./k_Qold_tmp;% divide out the old value and multiply in the new.
    tmp                 = diag([k_factor(1,1)*ones(M,1);k_factor(2,1)*ones(p*(M^2),1)]);
    Q_prmean_new        = tmp*Q_prmean*tmp; % k_Q changes the mean, so
else
    lpriordens_new(1,1) = priordens(k_Qnew(1,1), options_); %evaluate prior density with the new kQ draw
    Q_prmean_new        = (k_Qnew/k_Qold)^2*Q_prmean; % k_Q changes the mean, so
end;


[lik_new,loglik_new] = invwishpdf(Q,Q_prmean_new,Q_prvar);


loglik_new = loglik_new+sum(lpriordens_new); % posterior kernel of the new draw

% Calculate the acceptance ratio; Since the prior on k_Q is uniform, I have omitted it
logalpha = loglik_new-loglik_old;

% if using independent MH with half-cauchy proposal, need to add proposal densities. 
% For RW-MH and independent MH with uniform proposal, there is no need for this step
if options_.sampling==2 
    logalpha = logalpha+sum(dhalfcauchy(k_Qold_tmp,tuning_para))- sum(dhalfcauchy(k_Qnew,tuning_para)); % adding proposal densities for independent MH (half-cauchy)...
end;


%alpha = min(1,exp(logalpha));

% Accept according to Metropolis algorithm
if(log(rand)<logalpha) 
    if options_.mvs_sampling==1 % multivariate sampling of kQ
        k_Q = diag([k_Qnew(1,1)*ones(M,1);k_Qnew(2,1)*ones(p*(M^2),1)]);
    else
        k_Q=k_Qnew;
    end;
    accept_flag=1;
else
    accept_flag=0;
    k_Q = k_Qold;
end


