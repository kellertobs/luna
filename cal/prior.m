function [prior] = prior (model, bnds)

% return log prior probability of the model (scalar), assumes uniform prior
% model     vector of Nparams x 1 
% bnds      vector of Nparams x 2 [lower, upper]

pvar  = double(model>=bnds (:,1) & model<=bnds(:,2));
prior = sum(log(pvar)); 


end