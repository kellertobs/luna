function [L, tmp] = likefrommodel (model, cal0, c0, Tmp, Prs, stages, hasolv, haspxn, hasplg, hasspn, hasqtz, exp, sig, Tsol, Tliq, Psl, Tsl_sigma, wgt)
% 
% [L] = likefrommodel (model, cal0, oxds, Temp, Pres, stages, hasolv, haspxn, hasplg, hasspn, hasqtz, exp, sig)
% 
% calculates the likelihood of the predicted data given input model, needed
% to run catmip
% really just combines the runmodel and likelihood functions
% 
% 

mdl = runmodel(model, cal0, c0, Tmp, Prs, stages, hasolv, haspxn, hasplg, hasspn, hasqtz, Psl);
L   = likelihood(mdl, exp, sig, Tsol, Tliq, Tsl_sigma, wgt);
tmp = 1;

end