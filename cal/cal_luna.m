% set number of components
cal.nc  =  6;

% set string with component names
cal.CompStr  =  {'fo','fay','opx','cpx','an','ab'};
for i = 1:cal.nc; cal.(cal.CompStr{i}) = i; end

% set melting model end-member oxide compositions
cal.oxdStr =  {'SiO$_2$','Al$_2$O$_3$','FeO','MgO','CaO','Na$_2$O'};
%             SiO2  Al2O3    FeO    MgO    CaO   Na2O
cal.oxd = [    42.7    0.0    0.0   57.3    0.0    0.0
               30.0    0.0   69.0    0.0    1.0    0.0
               55.0    7.5    3.0   34.0    0.4    0.1
               49.0    2.0   27.7    5.1   16.0    0.2
               44.2   35.6    0.4    0.2   19.1    0.5
               67.4   20.4    1.5    0.1    0.6   10.0  ];
cal.oxd = cal.oxd./sum(cal.oxd,2)*100;

% specify calibration parameters

% set pure component melting points T_m^i at P=0
cal.T0(cal.fo )  =  1749;
cal.T0(cal.fay)  =  1049;
cal.T0(cal.opx)  =  1620;
cal.T0(cal.cpx)  =  960;
cal.T0(cal.an )  =  1295;
cal.T0(cal.ab )  =  1195;

% set first coeff. for P-dependence of T_m^i [GPa]
cal.A(cal.fo )   =  9.56;
cal.A(cal.fay)   =  5.88;
cal.A(cal.opx)   =  7.94;
cal.A(cal.cpx)   =  6.15;
cal.A(cal.an )   =  6.5;
cal.A(cal.ab )   =  6.1;

% set second coeff. for P-dependence of T_m^i [1]
cal.B(cal.fo )   =  1.73;
cal.B(cal.fay)   =  1.75;
cal.B(cal.opx)   =  1.58;
cal.B(cal.cpx)   =  1.91;
cal.B(cal.an )   =  1.79;
cal.B(cal.ab )   =  1.62;

% set entropy gain of fusion DeltaS [J/K]
cal.dS           =  300;

% set coeff. for T-dependence of partition coefficients K^i [1/K]
cal.r(cal.fo )   =  103;
cal.r(cal.fay)   =  87;
cal.r(cal.opx)   =  97;
cal.r(cal.cpx)   =  21;
cal.r(cal.an )   =  96;
cal.r(cal.ab )   =  31;
