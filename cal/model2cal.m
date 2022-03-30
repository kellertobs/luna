function [cal] = model2cal (cal0, model)

cal = cal0;

% set pure component melting points T_m^i at P=0
cal.T0(cal.fo )  =  model(1);
cal.T0(cal.fay)  =  model(2);
cal.T0(cal.opx)  =  model(3);
cal.T0(cal.cpx)  =  model(4);
cal.T0(cal.an )  =  model(5);
cal.T0(cal.ab )  =  model(6);

% set first coeff. for P-dependence of T_m^i [GPa]
cal.A(cal.fo )   = model(7);
cal.A(cal.fay)   = model(8);
cal.A(cal.opx)   = model(9);
cal.A(cal.cpx)   = model(10);
cal.A(cal.an )   = model(11);
cal.A(cal.ab )   = model(12);

% set second coeff. for P-dependence of T_m^i [1]
cal.B(cal.fo )   =  model(13);
cal.B(cal.fay)   =  model(14);
cal.B(cal.opx)   =  model(15);
cal.B(cal.cpx)   =  model(16);
cal.B(cal.an )   =  model(17);
cal.B(cal.ab )   =  model(18);

% set coeff. for T-dependence of partition coefficients K^i [1/K]
cal.r(cal.fo )   =  model(19);
cal.r(cal.fay)   =  model(20);
cal.r(cal.opx)   =  model(21);
cal.r(cal.cpx)   =  model(22);
cal.r(cal.an )   =  model(23);
cal.r(cal.ab )   =  model(24);

end