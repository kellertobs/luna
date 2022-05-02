function [mdls, pp, dhat, runtime, allmodels] = catmip (priorfunc, priorsampfunc, likefunc, varargin)
% [models, LLK] = catmip (priorfunc, priorsampfunc, likefunc, varargin)
%
% Uses the CATMIP algorithm to derive the posterior distribution for a
% non-linear inverse problems (tempering, resampling and metropolis)
% 
% From Minson et al., 2013: Bayesian inversion for finite fault earthquake
% source models I-theory and algorithm
%
% INPUTS
% priorfunc     calculates the log-prior probability of one model (output is scalar)     
% priorsampfunc SAMPLES the prior distribution (output is Nvar x Niter)
% likefunc      calculates the log-likelihood of one model (output is scalar)
% mbnds         model bnds [Nvar x 2]
% varargin      name-value pairs given in default_opts
% 
% OUTPUTS
% mdls      posterior distribution of models
% pp        posterior probability of models
% dhat      predicted data
% runtime   here i've only figured out how to store the first 1000 runtimes
% allmodels large array containing all models from each tempering step.
%
% YQW Dec 27, 2019. Based on Sarah Minson's mcmc demo code.
% Major edit 21 April 2022 to allow for different prior distributions
% 

% parse options
opts  = default_opts(varargin{:});
Cmsqr = opts.Cmsqr;
dbeta = opts.dbeta;
Nstps = opts.Nsteps;

totaltime = tic;

if (opts.Parallel)
    delete(gcp('nocreate'));
    ParpoolObj = parpool(min([opts.Ncores,opts.Niter])); 
    NumWorkers = ParpoolObj.NumWorkers;
else
    NumWorkers = 0;
end

% set initial tempering values
beta = 0; c = 0; m = 0;

% start by sampling from the prior distribution
mdls  = priorsampfunc(opts.Niter);

% collect the length of models and data vectors
Nparam    = size(mdls,2);
[~, dtmp] = likefunc(mdls(1,:));
Ndata     = length(dtmp);

% initialize variables
pp      = zeros(opts.Niter, 1    );  % posterior probability
dhat    = zeros(opts.Niter, Ndata);  % predicted data
runtime = zeros(opts.Niter, 1    );  % collect runtime


% run forward models
fprintf('Running initial %d models from prior...\n', opts.Niter);
parfor (mi = 1:opts.Niter, NumWorkers)
    tm = tic;
    pr = priorfunc(mdls(mi,:));
    
    if pr == -inf   
        % bad params, don't even calculate forward model
        pp(mi) = pr;
    else
        [lk, dhat(mi,:)] = likefunc(mdls(mi,:));
        pp(mi)           = lk + pr;
    end
    
    runtime(mi)      = toc(tm);
end
fprintf('Finished initial %d models from prior. Time taken = %.4f hours.\n\n',...
    opts.Niter, sum(runtime)/3600);

fprintf('m\tCm^2\tCOV\tbeta\t\tNaccept\t\tNreject\t\tAcceptRatio(pct)\n');
fprintf('--------------------------------------------');
fprintf('--------------------------------------------\n');
fprintf('%d\t%.4f\t%.4f\t%.2e\t%.4e\t%.4e\n',m,Cmsqr,c,beta,0,0);

if (opts.SaveFile), save(opts.FileName, 'mdls','pp','dhat','beta','RunTime');end

% track an array of all models to see transition of posterior during tempering
allmodels = zeros(opts.Niter, Nparam, 20);
allmodels(:,:,1) = mdls;

% start loop for tempering
while beta<1
    m = m+1; % tempering step
    
    % update temperature
    [w,beta,dbeta,c] = calc_beta(pp,beta,dbeta);
    
    fprintf('%d\t%.4f\t%.4f\t%.2e\t',m,Cmsqr,c,beta);
    
    % Resample to match new PriorFunc PDF
    count       = histcounts(rand([opts.Niter 1]),[0; cumsum(w)]);
    count(end-1)=sum(count(end-1:end));
    count       = count(1:end-1);
    
    ind = []; 
    for i=1:length(count)
        ind = [ind; repmat(i,count(i),1)]; 
    end
    
    % reassign models, pp and dhat according to resampling
    mdls  = mdls(ind,:);
    pp    = pp(ind);
    dhat  = dhat(ind,:);
    
    
    
    % now update samples using a Metropolis chain to explore model space
    % using a Gaussian proposal density with covariance Sm (to make z below)
    % PROCEDURE:
    % 1. Calculate p=w/sum(w)
    % 2. Calculate the expected value: E = sum(p_i*models_i)
    % 3. Calcluate Sm = sum{p_i*models_i*models_i^T} - E*E^T
    % 4. Return Cm^2 * Sm
    p  = w/sum(w);
    E  = sum(repmat(p,1,Nparam).*mdls, 1);
    Sm = zeros(Nparam);
    for i = 1:opts.Niter
        Sm = Sm + p(i)*mdls(i,:)'*mdls(i,:); 
    end
    Sm = Sm - E'*E;
    Sm = Cmsqr*Sm;
    Sm = 0.5*(Sm + Sm'); % Make sure that Sm is symmetric
    
    % Run opts.Nsteps of Metropolis on each sample
    IOacc = zeros(opts.Niter, opts.Nsteps-1);
    
    % Loop over samples: each sample is the seed for a Markov chain
    parfor (ii = 1:opts.Niter, NumWorkers)
        
        Xmd = zeros(Nstps  , Nparam);      % matrix of models
        Xpp = zeros(Nstps  , 1     );      % posterior probabilities
        Xdh = zeros(Nstps  , Ndata );      % store data
        Xac = zeros(Nstps-1, 1     );      % accepted model? 1 or 0
        
        % Our current sample is the seed for the chain
        Xmd(1,:) = mdls(ii,:);
        Xpp(1  ) =   pp(ii);
        Xdh(1,:) = dhat(ii,:);
        z        = mvnrnd(zeros(1,Nparam), Sm, Nstps);
        
        % preset dy for the chain just for initialization
        dy = Xdh(1,:);
        
        % Run Metropolis
        for k=1:Nstps-1 
            % current step
            x  = Xmd(k,:);
            px = Xpp(k  );
            dx = Xdh(k,:);
            
            % proposed step
            y   = Xmd(k,:) + z(k,:);
            
            % calculate prior probability
            pr = priorfunc(y);
            
            if pr==-inf
                % bad params, don't even calculate forward model
                py = pr;
            else
                % calc likelihood and prior == posterior probability of y
                [lk,dy] = likefunc(y);
                py      = lk + pr;
            end
            
            % compare posterior probabilities
            r   = beta*(py-px);
            u   = log(rand);
            
            if u<=r
                Xmd(k+1,:) = y; 
                Xpp(k+1,:) = py; 
                Xdh(k+1,:) = dy;
                Xac(k    ) = 1;
            else
                Xmd(k+1,:) = x; 
                Xpp(k+1,:) = px; 
                Xdh(k+1,:) = dx;
            end
        end
        
        % Save only the last sample from the Markov chain
        % Now the original sample has been updated by MCMC
        mdls(ii,:)   = Xmd(end,:); 
        pp(ii)       = Xpp(end  ); 
        dhat(ii,:)   = Xdh(end,:);
        IOacc(ii,:)  = Xac';
    end
    
    allmodels(:,:,m+1) = mdls;
    
    Nacc     = sum(IOacc(:));
    Nrej     = length(IOacc(:)) - Nacc;
    accRatio = Nacc/(Nacc + Nrej);
    
    fprintf('%.4e\t%.4e\t%.2f\n', Nacc, Nrej, 100*accRatio);
    
    % Rescale step size by acceptance rate per Matt Muto
    kc = (8*accRatio + 1)/9;
    % //kc = max(kc,0.2);   kc = min(kc,1);
    if (kc < 0.001); kc = 0.001; end
    if (kc > 1.000); kc = 1.000; end
    Cmsqr = kc * kc;
    
    if (opts.SaveFile)
        allmodelstmp = allmodels(:,:,1:m+1);
        save(opts.FileName, 'mdls', 'allmodelstmp', 'pp', 'dhat', 'beta', 'RunTime');
    end
    
    if (1-beta < 0.005); fprintf('mstop=%d\n\n',m); break; end
    
end

allmodels(:,:,m+2:end) = [];

% shut down parallel pool
if (opts.Parallel), delete(ParpoolObj); end

totalrt = toc(totaltime);
fprintf('Inversion completed. Total time taken = %.4f hours.\n\n',sum(totalrt)/3600);


end

function opts = default_opts (varargin)

p = inputParser;

p.addParameter('Niter'   ,  1000    , @isnumeric);
p.addParameter('Nsteps'  ,  5       , @isnumeric);
p.addParameter('Cmsqr'   ,  0.1^2   , @isnumeric);
p.addParameter('dbeta'   ,  1e-5    , @isnumeric);
p.addParameter('Parallel',  false   , @islogical);
p.addParameter('Ncores'  ,  6       , @isnumeric);
p.addParameter('SaveFile',  false   , @islogical);
p.addParameter('FileName',  ''      , @ischar);
p.addParameter('Ndatasets', 1       , @isnumeric);

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
beta = beta + dbeta2;
w    = calc_w(dbeta2, LLK);
c    = std(w,'omitnan')/mean(w,'omitnan');
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
s    = xmax + log(sum(exp(x-xmax)));
end
