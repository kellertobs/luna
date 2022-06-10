% set number of components
cal.nc  =  6;

% set string with component names
cal.CompStr  =  {'fo','fay','opx','cpx','an','ab'};
for i = 1:cal.nc; cal.(cal.CompStr{i}) = i; end

% set melting model end-member oxide compositions
cal.oxdStr =  {'SiO2','Al2O3','FeO','MgO','CaO','Na2O'};
cal.oxds   = [  42.7    0.0    0.0   57.3    0.0      0
                30.0    0.0   70.0    0.0    1.0      0
                55.0    7.5    3.0   34.0    0.4    0.1
                48.0    2.0   28.7    5.1   16.0    0.2
                44.2   35.6    0.4    0.2   19.1    0.5
                67.4   20.4    1.5    0.1    0.6   10.0  ];
cal.oxds = cal.oxds./sum(cal.oxds,2)*100;

% specify calibration parameters

% set pure component melting points T_m^i at P=0
cal.T0(cal.fo )  =  1700;
cal.T0(cal.fay)  =  1050;
cal.T0(cal.opx)  =  1300;
cal.T0(cal.cpx)  =  1000;
cal.T0(cal.an )  =  1150;
cal.T0(cal.ab )  =   900;

% set first coeff. for P-dependence of T_m^i [GPa]
cal.A(cal.fo )   =  11;
cal.A(cal.fay)   =   6;
cal.A(cal.opx)   =   4;
cal.A(cal.cpx)   =   3;
cal.A(cal.an )   =   8;
cal.A(cal.ab )   =   4;

% set second coeff. for P-dependence of T_m^i [1]
cal.B(cal.fo )   =  1.9;
cal.B(cal.fay)   =  1.7;
cal.B(cal.opx)   =  1.9;
cal.B(cal.cpx)   =  1.9;
cal.B(cal.an )   =  2.1;
cal.B(cal.ab )   =  2.0;

% set entropy gain of fusion DeltaS [J/K]
cal.dS           =  300;

% set coeff. for T-dependence of partition coefficients K^i [1/K]
cal.r(cal.fo )   =  64;
cal.r(cal.fay)   =  51;
cal.r(cal.opx)   =  18;
cal.r(cal.cpx)   =  13;
cal.r(cal.an )   =   5;
cal.r(cal.ab )   =  19;
