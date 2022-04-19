function [models, LLK, dhat, RunTime, allmodels] = catmip (PriorSampFunc, LikeFunc, mbnds, varargin)
% [models, LLK] = catmip (PriorSampFunc, LikeFunc, varargin)
%
% Uses the CATMIP algorithm to derive the posterior distribution for a
% non-linear inverse problems (tempering, resampling and metropolis)
% 
% From Minson et al., 2013: Bayesian inversion for finite fault earthquake
% source models I-theory and algorithm
%
% INPUTS
% PriorSampFunc SAMPLES the prior distribution (output is Nvar x Niter)
% LikeFunc      Calculates the log-likelihood of one model (output is scalar)
% mbnds         model bnds [Nvar x 2]
% varargin      name-value pairs given in default_opts
% 
% OUTPUTS
% models    posterior distribution of models
% LLK       log-likelihood of models
% dhat      predicted data
% RunTime   here i've only figured out how to store the first 1000 runtimes
% allmodels large array containing all models from each tempering step.
%
% YQW Dec 27, 2019. Based on Sarah Minson's mcmc demo code.
% 

% parse options
opts  = default_opts(varargin{:});
Cmsqr = opts.Cmsqr;
dbeta = opts.dbeta;

if (opts.Parallel)
    ParpoolObj = parpool(min([opts.Ncores,opts.Niter])); 
    NumWorkers = ParpoolObj.NumWorkers;
else
    NumWorkers = 0;
end

% preassign some variables to avoid broadcast variables in parfor loop
Nsteps     = opts.Nsteps;
Ndatasets  = opts.Ndatasets;
lowerbound = mbnds(:,1);
upperbound = mbnds(:,2);

% set initial tempering values
beta = 0; c = 0; m = 0;

% start by sampling from the prior distribution
models  = PriorSampFunc(opts.Niter);

% collect the length of models and data vectors
Nparam  = size(models,2);
[~, dhattmp] = LikeFunc(models(1,:));
Ndata = length(dhattmp);

% initialize variables
LLK     = zeros(opts.Niter, 1);      % log likelihood
dhat    = zeros(opts.Niter, Ndata);  % predicted data
RunTime = zeros(opts.Niter, 1);      % collect runtime


% run forward models
fprintf('Running initial %d models from prior...\n', opts.Niter);
parfor (mi = 1:opts.Niter, NumWorkers)
    tm = tic;
    [LLK(mi), dhat(mi,:)] = LikeFunc(models(mi,:));
    RunTime(mi) = toc(tm);
end
fprintf('Finished initial %d models from prior. Time taken = %.4f hours.\n\n',...
    opts.Niter, sum(RunTime)/3600);

fprintf('m\tCm^2\tCOV\tbeta\t\tNaccept\t\tNreject\n');
fprintf('----------------------------------');
fprintf('----------------------------------\n');
fprintf('%d\t%.4f\t%.4f\t%.2e\t%.4e\t%.4e\n',m,Cmsqr,c,beta,0,0);

if (opts.SaveFile)
    save(opts.FileName, 'models', 'LLK', 'dhat', 'beta', 'RunTime');
end

% track an array of all models to see transition of posterior during tempering
allmodels = zeros(opts.Niter, Nparam, 1000);
allmodels(:,:,1) = models;

% start loop for tempering
while beta<1
    m = m+1; % tempering step
    
    % update temperature
    [w,beta,dbeta,c] = calc_beta(LLK,beta,dbeta);
    
    fprintf('%d\t%.4f\t%.4f\t%.2e\t',m,Cmsqr,c,beta);
    
    % Resample to match new PriorFunc PDF
    count       = histcounts(rand([opts.Niter 1]),[0; cumsum(w)]);
    count(end-1)=sum(count(end-1:end));
    count       = count(1:end-1);
    
    ind = []; 
    for i=1:length(count)
        ind = [ind; repmat(i,count(i),1)]; 
    end
    
    models  = models(ind,:);
    LLK     = LLK(ind);
    dhat    = dhat(ind,:);
    
    % now update samples using a Metropolis chain to explore model space.
    % PROCEDURE:
    % 1. Calculate p=w/sum(w)
    % 2. Calculate the expected value: E = sum(p_i*models_i)
    % 3. Calcluate Sm = sum{p_i*models_i*models_i^T} - E*E^T
    % 4. Return Cm^2 * Sm
    p  = w/sum(w);
    E  = sum(repmat(p,1,Nparam).*models, 1);
    Sm = zeros(Nparam);
    for i=1:opts.Niter
        Sm = Sm+p(i)*models(i,:)'*models(i,:); 
    end
    Sm = Sm-E'*E;
    Sm = Cmsqr*Sm;
    Sm = 0.5*(Sm + Sm'); % Make sure that Sm is symmetric
    
    % Run opts.Nsteps of Metropolis on each sample
    IOacc = zeros(opts.Niter, opts.Nsteps-1);
    
    % Loop over samples: each sample is the seed for a Markov chain
    parfor (ii = 1:opts.Niter, NumWorkers)
        
        X       = zeros(Nsteps  , Nparam);
        Xllk    = zeros(Nsteps  , 1     );
        Xdhat   = zeros(Nsteps  , Ndata );
        XIOacc  = zeros(Nsteps-1, 1     );
        
        % Our current sample is the seed for the chain
        X(1,:)      = models(ii,:);
        Xllk(1)     = LLK(ii);
        Xdhat(1,:)  = dhat(ii,:);
        z           = mvnrnd(zeros(1,Nparam), Sm, Nsteps);
        
        % Run Metropolis
        for k=1:Nsteps-1 
            % current step
            x   = X(k,:);
            px  = Xllk(k);
            
            % proposed step
            y   = X(k,:) + z(k,:);
            
            if all(y(:)>=lowerbound) && all(y(:)<=upperbound)
                [py,dy] = LikeFunc(y);
            else
                py = -1e10;
                dy = cell(1,Ndatasets);
            end
            
            % compare posterior probabilities
            r   = beta*(py-px);
            u   = log(rand);
            
            if u<=r
                X(k+1,:)     = y; 
                Xllk(k+1,:)  = py; 
                Xdhat(k+1,:) = dy;
                XIOacc(k)    = 1;
            else
                X(k+1,:)     = x; 
                Xllk(k+1,:)  = px; %Nreject=Nreject+1;
                Xdhat(k+1,:) = Xdhat(k,:);
            end
        end
        
        % Save only the last sample from the Markov chain
        % Now the original sample has been updated by MCMC
        models(ii,:) = X(end,:); 
        LLK(ii)      = Xllk(end); 
        dhat(ii,:)   = Xdhat(end,:);
        IOacc(ii,:)  = XIOacc';
    end
    
    
    Naccept = sum(IOacc(:));
    Nreject = length(IOacc(:))-Naccept;
    allmodels(:,:,m+1) = models;

    fprintf('%.4e\t%.4e\n',Naccept,Nreject);
    
    % Rescale step size by acceptance rate per Matt Muto
    accRatio = Naccept/(Naccept + Nreject);
    kc = (8*accRatio + 1)/9;
    % //kc = max(kc,0.2);   kc = min(kc,1);
    if (kc < 0.001); kc = 0.001; end
    if (kc > 1.0); kc = 1.0; end
    Cmsqr = kc * kc;
    
    if (opts.SaveFile)
        allmodelstmp = allmodels(:,:,1:m+1);
        save(opts.FileName, 'models', 'allmodelstmp', 'LLK', 'dhat', 'beta', 'RunTime');
    end
    
    if (1-beta < 0.005); fprintf('mstop=%d\n',m); break; end
    
end

allmodels(:,:,m+2:end) = [];

end

function opts = default_opts (varargin)

p = inputParser;

p.addParameter('Niter',     1000,   @isnumeric);
p.addParameter('Nsteps',    5,      @isnumeric);
p.addParameter('Cmsqr',     0.1^2,  @isnumeric);
p.addParameter('dbeta',     1e-5,   @isnumeric);
p.addParameter('Parallel',  false,  @islogical);
p.addParameter('Ncores',    6,      @isnumeric);
p.addParameter('SaveFile',  false,  @islogical);
p.addParameter('FileName',  '',     @ischar);
p.addParameter('Ndatasets', 1,      @isnumeric);

p.parse(varargin{:});
opts = p.Results;
if ~isempty(opts.FileName), opts.SaveFile = true; end 

end

function [w,beta,dbeta,c] = calc_beta (LLK,beta,dbeta)
% calculates the next optimum cooling temperature. Most efficient sampling 
% achieved when new beta is chosen so that cov(w)=1 (Beck and Zuev, 2013)

try
    dbeta2 = fzero(@(db) (calc_c(db, LLK)-1), [0,1-beta]);
catch
    dbeta2 = 1-beta;
end

% dbeta2  = min([dbeta2, 1-beta]);
beta    = beta + dbeta2;

w = calc_w(dbeta2, LLK);
c = std(w,'omitnan')/mean(w,'omitnan');
end

function c = calc_c (dbeta, LLK)
w = calc_w(dbeta, LLK);
c = std(w,'omitnan')/mean(w,'omitnan');
end

function w = calc_w (dbeta, LLK)
logw = dbeta*LLK;
w    = exp(logw - logsumexp(logw)); % need to normalize w
end

function s = logsumexp(x)
xmax = max(x,[],'omitnan');
s = xmax + log(sum(exp(x-xmax)));
end
