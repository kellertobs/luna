function [ms] = priorsamp (mbnds, Niter)
% 
% [ms] = priorsamp (mbnds, Niter)
% samples a set of models from a uniform prior distribution
% 
% INPUTS
% mbnds     matrix of lower, upper bounds on prior [Nvar x 2]
% Niter     number of samples
% 
% OUTPUTS
% ms        matrix of model parameters [Niter x Nvar]

% check size of mbnds
if size(mbnds,2)~=2, mbnds = mbnds'; end

lower = repmat(mbnds(:,1)', Niter, 1);
upper = repmat(mbnds(:,2)', Niter, 1);
ms    = unifrnd(lower, upper);
end
