% set number of components
cal.nc  =  6;

% set string with component names
cal.CompStr  =  {'fo','fay','opx','cpx','an','ab'};
for i = 1:cal.nc; cal.(cal.CompStr{i}) = i; end

% set melting model end-member oxide compositions
cal.oxdStr =  {'SiO$_2$','Al$_2$O$_3$','FeO','MgO','CaO','Na$_2$O'};
cal.oxd    = [  42.7    0.0    0.0   57.3    0.0    0.0
                30.0    0.0   70.0    0.0    1.0    0.0
                55.0    7.5    3.0   34.0    0.4    0.1
                49.0    2.0   27.7    5.1   16.0    0.2
                44.2   35.6    0.4    0.2   19.1    0.5
                70.6   18.4    1.4    0.1    0.5    9.0 ];  % 0.9 alb + 0.1 qtz
%                 67.4   20.4    1.5    0.1    0.6   10.0  ];
cal.oxd = cal.oxd./sum(cal.oxd,2)*100;

% specify calibration parameters

% set pure component melting points T_m^i at P=0
cal.T0(cal.fo )  =  1891;
cal.T0(cal.fay)  =  1197;
cal.T0(cal.opx)  =  1538;
cal.T0(cal.cpx)  =  1012;
cal.T0(cal.an )  =  1374;
cal.T0(cal.ab )  =   992;

% set first coeff. for P-dependence of T_m^i [GPa]
cal.A(cal.fo )   =  11.6465;
cal.A(cal.fay)   =  15.7713;
cal.A(cal.opx)   =   9.9248;
cal.A(cal.cpx)   =  13.6101;
cal.A(cal.an )   =  12.0781;
cal.A(cal.ab )   =  14.7734;

% set second coeff. for P-dependence of T_m^i [1]
cal.B(cal.fo )   =  3.7265;
cal.B(cal.fay)   =  1.5034;
cal.B(cal.opx)   =  2.9294;
cal.B(cal.cpx)   =  1.9529;
cal.B(cal.an )   =  2.0868;
cal.B(cal.ab )   =  1.7408;

% set entropy gain of fusion DeltaS [J/K]
cal.dS           =  350;

% set coeff. for T-dependence of partition coefficients K^i [1/K]
cal.r(cal.fo )   =  64.7191;
cal.r(cal.fay)   =  48.6601;
cal.r(cal.opx)   =  57.4269;
cal.r(cal.cpx)   =  10.0294;
cal.r(cal.an )   =  49.4751;
cal.r(cal.ab )   =  26.5824;

