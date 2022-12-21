% set number of components
cal.nc  =  6;

% set string with component names
cal.CompStr  =  {'fo','fay','opx','cpx','an','ab'};
for i = 1:cal.nc; cal.(cal.CompStr{i}) = i; end

% set melting model end-member oxide compositions
cal.oxdStr =  {'SiO$_2$','Al$_2$O$_3$','FeO','MgO','CaO','Na$_2$O'};
cal.oxds   = [  42.7    0.0    0.0   57.3    0.0      0
                30.0    0.0   70.0    0.0    1.0      0
                55.0    7.5    3.0   34.0    0.4    0.1
                49.0    2.0   27.7    5.1   16.0    0.2
                44.2   35.6    0.4    0.2   19.1    0.5
                71.6   17.8    1.1    0.1    2.4    7.0 ];
%                 67.4   20.4    1.5    0.1    0.6   10.0  ];
cal.oxds = cal.oxds./sum(cal.oxds,2)*100;

% set bulk composition
cal.c0 = bulkmoon;

% specify calibration parameters

% set pure component melting points T_m^i at P=0
cal.T0(cal.fo )  =  1890;
cal.T0(cal.fay)  =  1205;
cal.T0(cal.opx)  =  1550;
cal.T0(cal.cpx)  =  1050;
cal.T0(cal.an )  =  1400;
cal.T0(cal.ab )  =   950;

% set first coeff. for P-dependence of T_m^i [GPa]
cal.A(cal.fo )   =  10.8;
cal.A(cal.fay)   =  15.8;
cal.A(cal.opx)   =  10.0;
cal.A(cal.cpx)   =  12.7;
cal.A(cal.an )   =  12.0;
cal.A(cal.ab )   =  14.2;

% set second coeff. for P-dependence of T_m^i [1]
cal.B(cal.fo )   =  3.70;
cal.B(cal.fay)   =  1.59;
cal.B(cal.opx)   =  2.80;
cal.B(cal.cpx)   =  1.78;
cal.B(cal.an )   =  2.00;
cal.B(cal.ab )   =  1.66;

% set entropy gain of fusion DeltaS [J/K]
cal.dS           =  350;

% set coeff. for T-dependence of partition coefficients K^i [1/K]
cal.r(cal.fo )   =  60;
cal.r(cal.fay)   =  45;
cal.r(cal.opx)   =  55;
cal.r(cal.cpx)   =  15;
cal.r(cal.an )   =  45;
cal.r(cal.ab )   =  30;


% set up calibration parameter bounds for mcmc
bnds.T0(cal.fo ,:) = cal.T0(cal.fo ) + [-1,+1]* 10;
bnds.T0(cal.fay,:) = cal.T0(cal.fay) + [-1,+1]* 10;
bnds.T0(cal.opx,:) = cal.T0(cal.opx) + [-1,+1]* 20;
bnds.T0(cal.cpx,:) = cal.T0(cal.cpx) + [-1,+1]*100;
bnds.T0(cal.an ,:) = cal.T0(cal.an ) + [-1,+1]* 30;
bnds.T0(cal.ab ,:) = cal.T0(cal.ab ) + [-1,+1]*100;

bnds.A(cal.fo ,:)  = cal.A(cal.fo )  + [-1,+1]*10;
bnds.A(cal.fay,:)  = cal.A(cal.fay)  + [-1,+1]*10;
bnds.A(cal.opx,:)  = cal.A(cal.opx)  + [-1,+1]*10;
bnds.A(cal.cpx,:)  = cal.A(cal.cpx)  + [-1,+1]*10;
bnds.A(cal.an ,:)  = cal.A(cal.an )  + [-1,+1]*10;
bnds.A(cal.ab ,:)  = cal.A(cal.ab )  + [-1,+1]*10;

bnds.B(cal.fo ,:)  = cal.B(cal.fo )  + [-1,+1]*0.2;
bnds.B(cal.fay,:)  = cal.B(cal.fay)  + [-1,+1]*0.2;
bnds.B(cal.opx,:)  = cal.B(cal.opx)  + [-1,+1]*0.2;
bnds.B(cal.cpx,:)  = cal.B(cal.cpx)  + [-1,+1]*0.2;
bnds.B(cal.an ,:)  = cal.B(cal.an )  + [-1,+1]*0.2;
bnds.B(cal.ab ,:)  = cal.B(cal.ab )  + [-1,+1]*0.2;

bnds.r(cal.fo ,:)  = cal.r(cal.fo )  + [-1,+1]*5;
bnds.r(cal.fay,:)  = cal.r(cal.fay)  + [-1,+1]*5;
bnds.r(cal.opx,:)  = cal.r(cal.opx)  + [-1,+1]*5;
bnds.r(cal.cpx,:)  = cal.r(cal.cpx)  + [-1,+1]*5;
bnds.r(cal.an ,:)  = cal.r(cal.an )  + [-1,+1]*5;
bnds.r(cal.ab ,:)  = cal.r(cal.ab )  + [-1,+1]*5;

bnds.mat           = [bnds.T0; bnds.A; bnds.B; bnds.r];

% variable names
vname = strcat(repelem({'T0';'A';'B';'r'},cal.nc,1), ',', repmat(cal.CompStr',4,1));
