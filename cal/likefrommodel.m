function [L, tmp] = likefrommodel (model, cal0, oxds, Temp, Pres, stages, hasolv, haspxn, hasplg, exp, sig)
% 
% [L] = likefrommodel (model, cal0, oxds, Temp, Pres, stages, hasolv, haspxn, hasplg, exp, sig)
% 
% calculates the likelihood of the predicted data given input model, needed
% to run catmip
% really just combines the runmodel and likelihood functions
% 
% 

mdl = runmodel(model, cal0, oxds, Temp, Pres, stages, hasolv, haspxn, hasplg);
L   = likelihood(mdl, exp, sig);
tmp = 1;

end