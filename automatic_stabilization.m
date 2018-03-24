function [tuning_para_kQ,tuning_para_kW,tuning_para_kS]=automatic_stabilization(tuning_para_kQ,tuning_para_kW,tuning_para_kS,jsux_kQ,jsux_kW,jsux_kS,jrep,options_,options_kQ,options_kW,options_kS,M)
% For RW-MH sampling:
% adjusting the std of the proposal increment innovation. if acceptance
% ratio is too low (cfactor<1), decrease the std, increse otherwise
%
% For independent MH sampling with uniform proposal:
% adjusting the upper bound of the uniform(0,tuning_para). if acceptance
% ratio is too low (cfactor<1), decrease the upper bound, increse otherwise
%
% For independent MH sampling with half-cauchy proposal
% adjusting the scale of the half-cauchy distribution. if acceptance
% ratio is too low (cfactor<1), decrease the scale, increse otherwise
%
% the tuning algorithm is controlled with a relaxation parameter
% 0<relax<1, adjustment is "faster" for higher relaxation parameter

% Tuning phase 
if ~strcmp(options_kQ.prior,'fixed')
    test1 = jsux_kQ/jrep;
    %disp(['Average Acceptance Ratio kQ - Burn In: ' num2str(test1)])
    if test1<1e-6
        error('Average acceptance ratio for kQ during the tuning phase is too low')
    end;
    cfactor1= test1/options_.AcceptanceTarget;
    tuning_para_kQ =  (1-options_.relax)*tuning_para_kQ+options_.relax*(tuning_para_kQ.*cfactor1); 
end;

if ~strcmp(options_kW.prior,'fixed')
    test2 = jsux_kW/jrep;
    %disp(['Average Acceptance Ratio kW - Burn In: ' num2str(test2)])
    if test2<1e-6
        error('Average acceptance ratio for kW during the tuning phase is too low')
    end;
    cfactor2= test2/options_.AcceptanceTarget;    
    tuning_para_kW =  (1-options_.relax)*tuning_para_kW+options_.relax*(tuning_para_kW.*cfactor2);  
end;


if ~strcmp(options_kS.prior,'fixed') && M > 1 % AR or VAR?   
    test3 = jsux_kS/jrep;
    %disp(['Average Acceptance Ratio kS - Burn In: ' num2str(test3)])
    if test3<1e-6
        error('Average acceptance ratio for kS during the tuning phase is too low')
    end
    cfactor3= test3/options_.AcceptanceTarget;
    tuning_para_kS =  (1-options_.relax)*tuning_para_kS+options_.relax*(tuning_para_kS.*cfactor3);   
end;



