function [L] = likelihood (mdl, exp, sig)

% returns log likelihood of the model (scalar),
% assumes normally distributed data errors

if mdl.flag == 1
    % solution found
    data  = [exp.fphs(:); exp.oxds(:)];
    dhat  = [mdl.fphs(:); mdl.oxds(:)];
    sigma = [sig.fphs(:); sig.oxds(:)];
    
    L = sum(-0.5*((data-dhat)./sigma).^2);
    
else
    % newton solver did not converge, assign a very small likelihood
    L = -1e8;
end
end