function [cal] = model2cal (cal, model)

im  = 0;

% set pure component melting points T_m^i at P=0
for ic = 1:cal.ncmp
    im = im+1;
    cal.T0(ic)  =  model(im);
end

% set first coeff. for P-dependence of T_m^i [GPa]
for ic = 1:cal.ncmp
    im = im+1;
    cal.A(ic)  =  model(im);
end

% set second coeff. for P-dependence of T_m^i [1]
for ic = 1:cal.ncmp
    im = im+1;
    cal.B(ic)  =  model(im);
end

% set coeff. for T-dependence of partition coefficients K^i [1/K]
for ic = 1:cal.ncmp
    im = im+1;
    cal.r(ic)  =  model(im);
end

end