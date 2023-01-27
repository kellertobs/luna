function [L] = likelihood (mdl, exp, sig, Tsol, Tliq, Tsl_sigma, wgt)

% returns log likelihood of the model (scalar),
% assumes normally distributed data errors

if mdl.flag == 1
    % solution found

    % experimental data from Schmidt & Kraettli (2020)
    expt_data  = [exp.phs(:); exp.oxd(:)];
    expt_dhat  = [mdl.phs(:); mdl.oxd(:)];
    expt_sigma = [sig.phs(:); sig.oxd(:)];

    % solidus and liquidus info
    Tsl_data   = [    Tsol(:);     Tliq(:)];
    Tsl_dhat   = [mdl.Tsol(:); mdl.Tliq(:)];
    
    L = wgt(1)*sum( -0.5*((expt_data - expt_dhat)./expt_sigma).^2 ) + ...
        wgt(2)*sum( -0.5*(( Tsl_data -  Tsl_dhat)./ Tsl_sigma).^2 );
else
    % newton solver did not converge, assign a very small likelihood
    L = -1e8;
end
end