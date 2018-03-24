%% ------------------------------------------------------------------------
% TVP-VAR-SV function
% -------------------------------------------------------------------------
% options_: structures containing all estimation specifications, see note
%
% -------------------------------------------------------------------------

%% Part I. Load data and set specifications
% =========================================================================
t           = size(Y,1);            % t is the time-series observations of Y
M           = size(Y,2);            % M is the dimensionality of Y
p           = options_.p;           % # lags
numa        = M*(M-1)/2;            % non-0,non-1 elements of A_t
tau         = options_.tau;         % Size of training sample 
Qprior_spec = options_.Qprior_spec; % 0: Primiceri setup 1: Minimal setup

% VAR EQUATION
% -------------------------------------------------------------------------
ylag        = mlag2(Y,p);               % lagged Y is [T x M]. ylag is [T x (Mp)]
ylag        = ylag(p+tau+1:t,:);        % Form RHS matrix X_t = [1 y_t-1 y_t-2 ... y_t-k] for t=1:T
if options_.const == 1                  % including constants
    K           = M + p*(M^2);          % K is the number of elements in the state vector
    Z           = zeros((t-tau-p)*M,K); % Create Z_t matrix.
    for i = 1:t-tau-p
        ztemp = eye(M);
        for j = 1:p
            xtemp = ylag(i,(j-1)*M+1:j*M);
            xtemp = kron(eye(M),xtemp);
            ztemp = [ztemp xtemp];
        end
        Z((i-1)*M+1:i*M,:) = ztemp;
    end
else
    K           = p*(M^2);                  % K is the number of elements in the state vector
    Z           = zeros((t-tau-p)*M,K);     % Create Z_t matrix.
    for i = 1:t-tau-p
        ztemp =[];
        for j = 1:p
            xtemp = ylag(i,(j-1)*M+1:j*M);
            xtemp = kron(eye(M),xtemp);
            ztemp = [ztemp xtemp];
        end
        Z((i-1)*M+1:i*M,:) = ztemp;
    end
    options_kQ.mvs_sampling=0;              % if no constants involved, multivariate kQ sampling is not allowed
end

y   = Y(tau+p+1:t,:)';
yFF = [Y mlag(Y,p-1)];
yFF = yFF(end,:)';
t   = size(y,2);                            % t is now smaller

%% Specify posterior sampler
% -------------------------------------------------------------------------
nrep        = options_.nrep;    % Number of replications
nburn       = options_.nburn;   % Number of burn-in replications (automatic tuning only during burn-in phase)
nsave       = options_.nsave;   % Number of draws saved in each file
nadjust     = options_.nadjust; % tuning phase adjustment sample
it_print    = 5000;             % Print in the screen every "it_print"-th iteration

try options_.forehor % BAD SOLUTION, modify!
catch 
    options_.forehor=0; 
end

draw_index      = sort(randsample(nrep,nsave))+nburn;   % randomly saves nsave number of draws after burn-in
tuning_para_kQ  = options_kQ.tuning_para_init;          % initializing tuning parameters for the automatic stabilization algorithm
tuning_para_kW  = options_kW.tuning_para_init;
tuning_para_kS  = options_kS.tuning_para_init;


%% PRIORS
% -------------------------------------------------------------------------
% specifying lower and upper bound of the prior density
[options_kQ.lower_bound,options_kQ.upper_bound] = prior_bounds(options_,options_kQ);
[options_kW.lower_bound,options_kW.upper_bound] = prior_bounds(options_,options_kW);
[options_kS.lower_bound,options_kS.upper_bound] = prior_bounds(options_,options_kS);

% Training sample prior
if M==1;
    [B_OLS,VB_OLS,sigma_OLS] = ts_prior_ar(Y,tau,p,options_.const); % AR model
else
    [B_OLS,VB_OLS,A_OLS,sigma_OLS,VA_OLS] = ts_prior(Y,tau,M,p,options_.const);
end;

sizeW   = M;                % Size of matrix W
sizeS   = 1:M;              % Size of matrix S
sizeB   = K;                % Size of matrxi B

% Prior means and variances (_prmean / _prvar)
% -------------------------------------------------------------------------
% Kalman filter initial conditions for B(t), A(t) and (log) SIGMA(t).
B_0_prmean      = B_OLS;                    % B_0 ~ N(B_OLS, 4Var(B_OLS))
B_0_prvar       = 4*VB_OLS;
if M>1 % AR or VAR?
    A_0_prmean      = A_OLS;                % A_0 ~ N(A_OLS, 4Var(A_OLS))
    A_0_prvar       = 4*VA_OLS;
end;
sigma_prmean    = sigma_OLS;                % log(sigma_0) ~ N(log(sigma_OLS),I_n)
sigma_prvar     = 4*eye(M);
% Q ~ IW(k2_Q*size(tau)*Var(B_OLS),size(stau))

if options_kQ.mvs_sampling==1; % sampling scale for constants and VAR-coefficients separately
    % building the diag(kQ) matrix
    k_Q = diag([options_kQ.init*ones(M,1);options_kQ.init*ones(p*(M^2),1)]);
else
    k_Q =options_kQ.init;
end;

k_W =options_kW.init;
k_S =options_kS.init;

if Qprior_spec==0
        Q_prmean        = k_Q*tau*VB_OLS*k_Q;         
        Q_prvar         = tau;
else
        Q_prmean        = k_Q*(sizeB+1)*VB_OLS*k_Q;        
        Q_prvar         = sizeB+1;                         
end;

% W ~ IW(k2_W*(1+dim(W))*I_n,(1+dim(W)))
W_prmean        = ((options_kW.init)^2)*(1 + sizeW)*eye(M); % W_prmean = 0.00001*eye(M);
W_prvar         = 1 + sizeW;                    % W_prvar = tau;
% S ~ IW(k2_S*(1+dimension(S)*Var(A_OLS),(1+dimension(S)))
if M>1 % AR or VAR?
    S_prmean        = cell(M-1,1);
    S_prvar         = zeros(M-1,1);
    ind = 1;
    for ii = 2:M % S is block diagonal as in Primiceri (2005)
        S_prmean{ii-1}  = ((options_kS.init)^2)*(1 + sizeS(ii-1))*VA_OLS(((ii-1)+(ii-3)*(ii-2)/2):ind,((ii-1)+(ii-3)*(ii-2)/2):ind);
        S_prvar(ii-1)   = 1 + sizeS(ii-1);
        ind             = ind + ii;
    end
end;

% Parameters of the 7 component mixture approximation to a log(chi^2) density:
q_s     = [  0.00730;  0.10556;  0.00002; 0.04395; 0.34001; 0.24566;  0.25750]; % probabilities
m_s     = [-10.12999; -3.97281; -8.56686; 2.77786; 0.61942; 1.79518; -1.08819]; % means
u2_s    = [  5.79596;  2.61369;  5.17950; 0.16735; 0.64009; 0.34023;  1.26261]; % variances

%% INITIALIZE MATRICES:
% -------------------------------------------------------------------------
% Specify covariance matrices for measurement and state equations
consQ   = 0.0001;
consS   = 0.0001;
consH   = 0.01;
consW   = 0.0001;

Ht          = kron(ones(t,1),consH*eye(M));         % Initialize Htdraw, a draw from the VAR covariance matrix
Htchol      = kron(ones(t,1),sqrt(consH)*eye(M));   % Cholesky of Htdraw defined above
Qdraw       = consQ*eye(K);                         % Initialize Qdraw, a draw from the covariance matrix Q
if M>1 % AR or VAR?
    Atdraw      = zeros(numa,t);                    % Initialize Atdraw, a draw of the non 0 or 1 elements of A(t)    
    Sdraw       = consS*eye(numa);                  % Initialize Sdraw, a draw from the covariance matrix S
    Sblockdraw  = cell(M-1,1);                      % ...and then get the blocks of this matrix (see Primiceri)
    ijc = 1;
    for jj=2:M
        Sblockdraw{jj-1} = Sdraw(((jj-1)+(jj-3)*(jj-2)/2):ijc,((jj-1)+(jj-3)*(jj-2)/2):ijc);
        ijc = ijc + jj;
    end
else
    Atdraw = 1;
end;
Wdraw       = consW*eye(M);                 % Initialize Wdraw, a draw from the covariance matrix W
Btdraw      = zeros(K,t);                   % Initialize Btdraw, a draw of the mean VAR coefficients, B(t)
Sigtdraw    = zeros(M,t);                   % Initialize Sigtdraw, a draw of the log-diagonal of SIGMA(t)
sigt        = kron(ones(t,1),0.01*eye(M));  % Matrix of the exponent of Sigtdraws (SIGMA(t))
statedraw   = 5*ones(t,M);                  % initialize the draw of the indicator variable (of 7-component mixture of Normals approximation)
Zs          = kron(ones(t,1),eye(M));
prw         = zeros(numel(q_s),1);
forecast    = zeros(M,options_.forehor,nsave);

%% if sampling is continued from a previous simulation, please provide 'postdraws_old_init'
if options_.sampling_from_old_init==1
    load([save_path 'postdraws_old_init']);
    options_.acc_tuning=0; % switching automatic stabilization OFF!!!
end;

%% Storage matrices for posteriors and stuff
% -------------------------------------------------------------------------
Bt_post     = zeros(K,t,nsave);             % regression coefficients B(t)
Sigt_post   = zeros(M,t,nsave);             % diagonal std matrix SIGMA(t)
Q_post      = zeros(K,K,nsave);             % covariance matrix Q of B(t)
if M>1 % AR or VAR?
    At_post     = zeros(numa,t,nsave);      % lower triangular matrix A(t)
    ikc = 1;
    for kk = 2:M
        Sdraw(((kk-1)+(kk-3)*(kk-2)/2):ikc,((kk-1)+(kk-3)*(kk-2)/2):ikc)=Sblockdraw{kk-1};
        ikc = ikc + kk;
    end
    S_post      = zeros(numa,numa,nsave);   % covariance matrix Q of B(t)
    cor_post    = zeros(t,numa,nsave);      % off-diagonal elements of the VAR cov matrix
end;
W_post      = zeros(M,M,nsave);             % covariance matrix Q of B(t)
sig_post    = zeros(t,M,nsave);             % diagonal of the VAR covariance matrix

if options_kQ.mvs_sampling==1;
    k_Q_post     = zeros(nrep+nburn,2);
else
    k_Q_post     = zeros(nrep+nburn,1);
end;
    k_W_post     = zeros(nrep+nburn,1);    
    k_S_post     = zeros(nrep+nburn,1);

%% II. Start posterior sampler
% =========================================================================
tic; 
disp('Number of iterations');

iirep   = 1;
jrep    = 0;

jsux_kQ = 0; kQ_acc_flag = 0; isux_kQ = 0;
jsux_kW = 0; kW_acc_flag = 0; isux_kW = 0;
jsux_kS = 0; kS_acc_flag = 0; isux_kS = 0;

average_acc = [];

for irep = 1 : nrep + nburn    % GIBBS iterations starts here
    
    if mod(irep,it_print) == 0;
        fprintf('Draw:\t%8.0f (took %6.2f minutes)\n',irep,toc/60);
        disp(['Average Acceptance Ratio kQ - Current:' num2str(isux_kQ/irep)])
        disp(['Average Acceptance Ratio kW - Current:' num2str(isux_kW/irep)])
        if M>1
            disp(['Average Acceptance Ratio kS - Current:' num2str(isux_kS/irep)])
        end;
    end
    
    %% p(B_t,Q|Y,...): sample coeff. states and respective residual covariance
    % ---------------------------------------------------------------------
    draw_beta_corrected
    
    %% p(A_t,S|Y,...): sample cov. states and respective residual covariance
    % ---------------------------------------------------------------------
    draw_alpha_corrected
       
    %% Correction sampling for covariance states and SV
    % ---------------------------------------------------------------------
    [statedraw,yss,capAt]           = draw_sTcomp_corrected(Atdraw,Sigtdraw,yhat,m_s,u2_s,q_s,M,t);
    [Sigtdraw,log_lik3,sigt,Wdraw]  = draw_sigma_correctedSep(statedraw,Wdraw,yss,Zs,m_s,u2_s,M,t,sigma_prmean,sigma_prvar,W_prmean,W_prvar);
    
    %% Create the VAR covariance matrix H(t). It holds that: A(t) x H(t) x A(t)' = SIGMA(t) x SIGMA(t) '
    % ---------------------------------------------------------------------
    Ht      = zeros(M*t,M);
    Htsd    = zeros(M*t,M);
    for i = 1 : t
        stem                    = sigt((i-1)*M+1:i*M,:);
        Hsd                     = capAt((i-1)*M+1:i*M,:)\stem;
        Hdraw                   = Hsd*Hsd';
        Ht((i-1)*M+1:i*M,:)     = Hdraw;    % H(t)
        Htsd((i-1)*M+1:i*M,:)   = Hsd;      % Cholesky of H(t)
    end
    
    %% estimate kappas block
    % ---------------------------------------------------------------------
    if ~strcmp(options_kQ.prior,'fixed')
        [k_Q,kQ_acc_flag] = draw_k_Q(Qdraw,k_Q, Q_prmean,Q_prvar,tuning_para_kQ,options_kQ,M,p);
    end;
    if ~strcmp(options_kW.prior,'fixed')
         [k_W,kW_acc_flag] = draw_k_W(Wdraw,k_W, W_prmean,W_prvar,tuning_para_kW,options_kW,M,p);
    end;
    if M>1 % VAR or AR?
        if ~strcmp(options_kS.prior,'fixed')       
            [k_S,kS_acc_flag] = draw_k_S(Sblockdraw,k_S, S_prmean,S_prvar,tuning_para_kS,options_kS,M,p);
        end;
    end;
    
    %% automatic stabilization block
    % ---------------------------------------------------------------------
    if  irep/nadjust == round(irep/nadjust) && irep<=nburn && options_.acc_tuning==1 %% Tuning phase     
        if ~strcmp(options_kW.prior,'uniform')
            [tuning_para_kQ,tuning_para_kW,tuning_para_kS]=automatic_stabilization(tuning_para_kQ,tuning_para_kW,tuning_para_kS,jsux_kQ,jsux_kW,jsux_kS,jrep,options_,options_kQ,options_kW,options_kS,M);
        end
        jsux_kQ = 0;  jsux_kW = 0; jsux_kS = 0;jrep = 0;
    end
    
    jsux_kQ=jsux_kQ+kQ_acc_flag;
    jsux_kW=jsux_kW+kW_acc_flag;
    jsux_kS=jsux_kS+kS_acc_flag;
    
    isux_kQ=isux_kQ+kQ_acc_flag;
    isux_kW=isux_kW+kW_acc_flag;
    isux_kS=isux_kS+kS_acc_flag;    
    
    jrep = jrep + 1;
    
    %% Updating priors for Q,W,S with new scaling parameters
    % ---------------------------------------------------------------------
    % Q ~ IW(k2_Q*size(tau)*Var(B_OLS),size(stau))
    % 0: Primiceri setup 1: Minimal setup
    % ---------------------------------------------------------------------
    if Qprior_spec==0
        Q_prmean        = k_Q*tau*VB_OLS*k_Q;         
        Q_prvar         = tau;
    else
        Q_prmean        = k_Q*(sizeB+1)*VB_OLS*k_Q;        
        Q_prvar         = sizeB+1;                         
    end;        
    % W ~ IW(k2_W*(1+dim(W))*I_n,(1+dim(W)))
    % ---------------------------------------------------------------------
    W_prmean        = (k_W^2)*(1 + sizeW)*eye(M); % W_prmean = 0.00001*eye(M);
    W_prvar         = 1 + sizeW;                    % W_prvar = tau;
    if M>1 % AR or VAR?
        % S ~ IW(k2_S*(1+dimension(S)*Var(A_OLS),(1+dimension(S)))
        % -----------------------------------------------------------------
        S_prmean        = cell(M-1,1);
        S_prvar         = zeros(M-1,1);
        ind = 1;
        for ii = 2:M % S is block diagonal as in Primiceri (2005)
            S_prmean{ii-1}  = (k_S^2)*(1 + sizeS(ii-1))*VA_OLS(((ii-1)+(ii-3)*(ii-2)/2):ind,((ii-1)+(ii-3)*(ii-2)/2):ind);
            S_prvar(ii-1)   = 1 + sizeS(ii-1);
            ind             = ind + ii;
        end
    end;
    
    %% save all kQ, kW and kS
    % ---------------------------------------------------------------------
    if options_kQ.mvs_sampling==1;
        k_Q_post(irep,:)=[k_Q(1,1);k_Q(M+1,M+1)];
    else
        k_Q_post(irep,:)=k_Q;
    end;
    k_W_post(irep,1)=k_W;
    
    traceQ_post(irep,1)=trace(Qdraw);
    if M>1
        k_S_post(irep,1)=k_S;
        traceS_post(irep,1)=trace(Sdraw);
    end;
    traceW_post(irep,1)=trace(Wdraw);


    %% -------------------- SAVE AFTER-BURN-IN DRAWS ----------------------
     if ismember(irep,draw_index)
        % Get time-varying correlations and variances
        stemp6 = zeros(M,1);
        stemp5 = [];
        stemp7 = [];
        for i = 1:t
            stemp8 = corrvc(Ht((i-1)*M+1:i*M,:));
            stemp7a = [];
            ic = 1;
            for j = 1:M
                if j>1;
                    stemp7a = [stemp7a ; stemp8(j,1:ic)']; 
                    ic = ic+1;
                end
                stemp6(j,1) = sqrt(Ht((i-1)*M+j,j));
            end
            stemp5 = [stemp5 ; stemp6'];
            stemp7 = [stemp7 ; stemp7a'];
        end
        sig_post(:,:,iirep) = stemp5;           % diagonal of the VAR covariance matrix
        if M>1
            cor_post(:,:,iirep) = stemp7;       % off-diagonal elements of the VAR cov matrix
            At_post(:,:,iirep)      = Atdraw;   % lower triangular matrix A(t)
            ikc = 1;
            for kk = 2:M
                Sdraw(((kk-1)+(kk-3)*(kk-2)/2):ikc,((kk-1)+(kk-3)*(kk-2)/2):ikc)=Sblockdraw{kk-1};
                ikc = ikc + kk;
            end
            S_post(:,:,iirep) = Sdraw;          % covariance matrix S of A(t)
        end;
        
        Bt_post(:,:,iirep)      = Btdraw;       % regression coefficients B(t)       
        Sigt_post(:,:,iirep)    = Sigtdraw;     % diagonal std matrix SIGMA(t)
        Q_post(:,:,iirep)       = Qdraw;        % covariance matrix Q of B(t)
        W_post(:,:,iirep)       = Wdraw;        % covariance matrix W of SIGMA(t)
        
        %% ---------------    Forecast AFTER-BURN-IN DRAWS    -------------
        if options_.forehor>=1
            biga = zeros(M*p,M*p);
            for j = 1:p-1
                biga(j*M+1:M*(j+1),M*(j-1)+1:j*M) = eye(M);
            end
            if options_.const==1;               % 1: with constants
                bbtemp = Btdraw(M+1:end,end);   % get the draw of B(t) at time i=1,...,T  (exclude intercept)
            else
                bbtemp = Btdraw(1:end,end);
            end;
            splace = 0;
            for ii = 1:p
                for iii = 1:M
                    biga(iii,(ii-1)*M+1:ii*M) = bbtemp(splace+1:splace+M,1)';
                    splace = splace + M;
                end
            end
            if options_.const==1;        % 1: with constants
                intercepts=Btdraw(1:M,end);
            else
                intercepts=zeros(M,1);
            end;
            ctemp1 = biga;
            tmp=yFF; 
            for k=1:options_.forehor
                tmp= [intercepts;zeros(M*(p-1),1)]+ctemp1*tmp;
                forecast(:,k,iirep)=tmp(1:M,1); % forecast is a M x h x nsave array
            end;
        end
        iirep = iirep + 1;
    end % END saving after burn-in results
end %END main Gibbs loop (for irep = 1:nrep+nburn)

options_kQ.tuning_para_final = tuning_para_kQ; % saving final tuning parameters for the automatic stabilization algorithm
options_kW.tuning_para_final = tuning_para_kW;
options_kS.tuning_para_final = tuning_para_kS;

% Average accepatance ratio after burn in
% -------------------------------------------------------------------------
average_acc.kQ = isux_kQ/irep;
average_acc.kW = isux_kW/irep;
average_acc.kS = isux_kS/irep;
% disp(['Average Acceptance Ratio kW - Final: ' num2str(average_acc.kW)])
% disp(['Average Acceptance Ratio kQ - Final: ' num2str(average_acc.kQ)])
% if M>1
%     disp(['Average Acceptance Ratio kS - Final: ' num2str(average_acc.kS)])
% end;

% Save last set of last posterior draws (Potentially usefull fir future run of chains)
if M>1
    save([save_path save_str '_old_init'], 'Ht','Htsd','Wdraw', 'Sdraw', 'Qdraw', 'Sigtdraw', 'Atdraw', 'Btdraw', 'sigt', 'statedraw', 'Sblockdraw',...
        'k_Q','k_W','k_S','options_','options_kQ','options_kW','options_kS');
    save([save_path save_str 'kappas_posterior'],'k_Q_post','k_W_post','k_S_post','options_','options_kQ','options_kW','options_kS','average_acc',...
        'traceQ_post','traceW_post','traceS_post');
else
    save([save_path save_str '_old_init'], 'Ht','Htsd','Wdraw', 'Qdraw', 'Sigtdraw', 'Btdraw', 'sigt', 'statedraw',...
        'k_Q','k_W','k_S','options_','options_kQ','options_kW','options_kS');
    save([save_path save_str 'kappas_posterior'],'k_Q_post','k_W_post','options_','options_kQ','options_kW','average_acc',...
        'traceQ_post','traceW_post');
end;
    
% Save forecast means/medians
if options_.forehor>=1
    save([save_path save_str 'forecasts'],'forecast');
end

if options_.do_save
    if M>1
    save([save_path save_str], 'Wdraw', 'Sdraw', 'Qdraw', 'Sigtdraw', 'Atdraw', 'Btdraw', 'sigt', 'statedraw', 'Sblockdraw',...
        'k_Q_post','k_W_post','k_S_post','options_','options_kQ','options_kW','options_kS','average_acc','Bt_post','At_post',...
        'Sigt_post','Q_post','S_post','W_post','traceQ_post','traceW_post','traceS_post','sig_post','cor_post');
    else
            save([save_path save_str], 'Wdraw', 'Qdraw', 'Sigtdraw', 'Btdraw', 'sigt', 'statedraw',...
        'k_Q_post','options_','options_kQ','options_kW','options_kS','average_acc','Bt_post',...
        'Sigt_post','Q_post','W_post','traceQ_post','traceW_post','sig_post');
    end
end

% fprintf('\n')
% fprintf('====================================================================================\n')
% fprintf('Running full code took:\t%8.2f\t hours!\n',toc/3600)
% fprintf('====================================================================================\n')
% fprintf('\n')

