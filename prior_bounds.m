function [lb,ub]=prior_bounds(options_global,options_kappa)
% determine endogenously the lower and upperbound of the prior density for
% given parameterization. modelled after the dynare function "prior_bounds"

prior_trunc = options_global.prior_trunc;

switch options_kappa.prior
    
    case 'fixed'
        lb=NaN;
        ub=NaN;
    
    case 'inv-gamma' % inverse gamma (scaled-inverse-chi2 specification)
        a=options_kappa.df/2;
        b=options_kappa.scale*a;
        lb=options_global.hard_lb;
        ub=1/gaminv(1-(1-prior_trunc),a,b);
        %ub=b/(a+1)*500;
    case 'half-cauchy'
        lb=options_global.hard_lb;
        ub=qhalfcauchy(1-prior_trunc,options_kappa.scale);
    case 'half-t' 
        % using half-cauchy specification due to the difficulty 
        % of evaluating half-t quantile function
        lb=options_global.hard_lb;
        ub=qhalfcauchy(1-prior_trunc,options_kappa.scale);
    case 'uniform'
        lb=options_kappa.lower_bound;
        ub=options_kappa.upper_bound;
        
end;