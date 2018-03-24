%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Main script to replicate empirical results in
% 
%       "Choosing Prior Hyperparameters: With Applications To Time-Varying 
%       Parameter Models" by Amir-Ahmadi, Matthes and Wang (2018,JBES)
% 
% -------------------------------------------------------------------------
% Some notational translation 
% -------------------------------------------------------------------------
% Paper notation: \kappa_{{\Omega_{b}},\kappa_{{\Omega_{h}}, \kappa_{{\Omega_{a}}
% Respective code notation: k_Q, k_W, k_S
% 
% -------------------------------------------------------------------------
% Acknowledgements
% -------------------------------------------------------------------------
% These codes build on TVP-VAR codes kindly provided by Gary Koop and 
% Dimitris Korobilis on their respective websites. 
% 
% -------------------------------------------------------------------------
% Contact:
% -------------------------------------------------------------------------
% Pooyan Amir-Ahmadi	
% University of Illinois at Urbana-Champaign
% Email: pooyan@illinois.edu
% 
% Christian Matthes
% Federal Reserve Bank of Richmond
% Email: christian.matthes@rich.frb.org
% 
% Mu-Chun Wang
% Hamburg University
% Email: Mu-Chun.Wang@uni-hamburg.de
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear('all'); close('all'); clc

%% Key Choice here 
% =========================================================================
% I.    Data selection
% II.   Sampling selection
% III.  Prior selection
% =========================================================================
data_sel        = 9;    % 5: euroarea_data.mat    
                        % 9: uk_data.mat
                        
sampling_sel    = 0;    % 0: RW-MH, (Benchmark choice!)
                        % 1: I-MH with uniform proposal 
                        % 2: I-MH w/ half-cauchy proposal

prior_sel       = 1;    % 1: fixed (as in Priiceri (2005, REStud))
                        % 2: uniform prior
                        % 3: inv-gamma prior
                        % 4: half-cauchy prior
                        % 5: half-t prior
                        
switch prior_sel
    case 1
        prior_str = 'fixed';
    case 2
        prior_str = 'uniform';
    case 3
        prior_str = 'inv-gamma';
    case 4
        prior_str = 'half-cauchy';
    case 5
        prior_str = 'half-t';    
end

%% Organizing path and lebelling
main_path   = pwd;

if ispc
    data_path   = [main_path '\Data\'];
else
    data_path   = [main_path '/Data/'];
end

switch data_sel
    case 5
        savefile        = 'Euroarea';
        load([data_path 'euroarea_data'])
        Y               = data;
        options_.tau    = 40;
    case 9
        savefile        = 'UK';
        load([data_path 'uk_data'])
        Y               = data;
        options_.tau    = 24;
end

if ispc
    save_path   = [main_path '\Results\'  savefile '\' ];
    figure_path = [main_path '\Results\'  savefile '\' ];
else
    save_path   = [main_path '/Results/'  savefile '/' ];
    figure_path = [main_path '/Results/'  savefile '/' ];
end
mkdir(save_path);

cd(main_path)

%% global estimation parameters
options_.nrep                   =  2000;        % Number of replications
options_.nburn                  = 48000;        % Number of burn-in replications (automatic tuning only during burn-in phase)
options_.nsave                  =  2000;        % Number of draws saved 
options_.nadjust                =   500;        % tuning phase adjustment sample

options_.acc_tuning             =     1;        % 1: with automatic stabilization phase (default), 0: without automatic stabilization phase 
options_.Qprior_spec            =     1;        % 0: Primiceri setup 1: Minimal setup for Q 
options_.p                      =     2;        % VAR lags
% options_.tau                    =    95;        % Size of training sample 
options_.sampling_from_old_init =     0;        % 1: if the sampling should be continued from a previous MCMC sampling
options_.AcceptanceTarget       =   1/3;        % target average acceptance rate
options_.relax                  =   0.8;        % relaxation parameter, 0<relax<1, adjustment is "faster" for higher relaxation parameter
options_.do_save                =     1;        % 0: no saving of results   
options_.do_stability           =     1;        % 1: imposing stationarity condtion
options_.hard_lb                =   eps;        % hard coded lower bound for scale parameters
options_.prior_trunc            =   eps;        % prior truncation value
options_.const                  =     1;        % 1: with constants

%% Automatic Stabilization doesn't work with I-MH
if sampling_sel ~= 0 % 0: RW-MH,1: I-MH U-prop,2: I-MH w/ half-cauchy proposal
    options_.acc_tuning = 0;        % 1: with automatic stabilization phase (default), 0: without automatic stabilization phase
end

%% kQ parameters
% choice of prior distribution: 'fixed, 'uniform', 'inv-gamma', 'half-cauchy', 'half-t'
options_kQ.prior            = prior_str;
options_kQ.scale            = 0.1;          % scale parameter for 'inv-gamma', 'half-cauchy', 'half-t'
options_kQ.df               = 2;            % degree of freedom for 'inv-gamma','half-t'
options_kQ.lower_bound      = 1e-10;        % lower bound for the 'uniform'
options_kQ.upper_bound      = 1;            % upper bound for the 'uniform'
options_kQ.init             = 0.1;          % initial starting value
options_kQ.tuning_para_init = 0.001;        % initial starting value of the tuning parameter,
options_kQ.sampling         = sampling_sel; % sampling method: 0 random walk MH, 1 independent MH with uniform proposal, 2 independent MH with half-cauchy proposal
options_kQ.mvs_sampling     = 1;            % 1 for multivariate sampling of kQ (separate kQs for constants and VAR coefficients each)

%% kW parameters
% choice of prior distribution: 'fixed, 'uniform', 'inv-gamma', 'half-cauchy', 'half-t'
options_kW.prior            = prior_str;
options_kW.scale            = 0.1;          % scale parameter for 'inv-gamma', 'half-cauchy', 'half-t'
options_kW.df               = 2;            % degree of freedom for 'inv-gamma','half-t'
options_kW.lower_bound      = 1e-10;        % lower bound for the 'uniform'
options_kW.upper_bound      = 1;            % upper bound for the 'uniform'
options_kW.init             = 0.1;          % initial starting value
options_kW.tuning_para_init = 0.001;        % initial starting value of the tuning parameter
options_kW.sampling         = sampling_sel; % sampling method: 0 random walk MH, 1 independent MH with uniform proposal, 2 independent MH with half-cauchy proposal

%% kS parameters
% choice of prior distribution: 'fixed, 'uniform', 'inv-gamma', 'half-cauchy', 'half-t'
options_kS.prior            = prior_str;
options_kS.scale            = 0.1;          % scale parameter for 'inv-gamma', 'half-cauchy', 'half-t'
options_kS.df               = 2;            % degree of freedom for 'inv-gamma','half-t'
options_kS.lower_bound      = 1e-10;        % lower bound for the 'uniform'
options_kS.upper_bound      = 1;            % upper bound for the 'uniform'
options_kS.init             = 0.1;          % initial starting value
options_kS.tuning_para_init = 0.001;        % initial starting value of the tuning parameter
options_kS.sampling         = sampling_sel; % sampling method: 0 random walk MH, 1 independent MH with uniform proposal, 2 independent MH with half-cauchy proposal

if options_.do_stability
    save_str = sprintf('FinalMat_%s_Prior_%d_SamplingScheme_%d_StableDraws',savefile,prior_sel,sampling_sel);
else
    save_str = sprintf('FinalMat_%s_Prior_%d_SamplingScheme_%d',savefile,prior_sel,sampling_sel);
end

%% Run main code for TVP-VAR with hyperparameter estimation
estimate_TVVAR_func4
