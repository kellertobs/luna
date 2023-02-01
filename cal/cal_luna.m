% set number of components
cal.nc  =  6;

% set string with component names
cal.cmpStr  =  {'for','fay','opx','cpx','ant','eut'};
for i = 1:cal.nc; cal.(cal.cmpStr{i}) = i; end

% set mineral end-member oxide compositions
cal.oxdStr =  {'SiO$_2$','TiO$_2$','Al$_2$O$_3$','FeO','MgO','CaO','Na$_2$O'};
cal.mnr_oxd = [  42.7    0.05   0.0    0.0   57.25   0.0      0     % forsterite
                 30.0    0.20   0.0   69.0    0.0    0.8      0     % fayalite
                 55.6    0.10   5.0    3.9   34.3    1.0    0.1     % Mg-pyroxene
                 50.5    0.35   8.2   12.2   25.45   3.2    0.1     % Mg-Fe-Ca-pyroxene
                 48.3    1.20   2.5   26.0    5.7   16.1    0.2     % Fe-Ca-pyroxene
                 44.2    0.08  35.52   0.4    0.2   19.1    0.5     % anorthite
                 67.4    0.20  20.4    1.3    0.1    0.6   10.0     % albite
                  0.5   20.30   4.0   73.3    1.5    0.4    0.0     % spinel
                99.95    0.05   0.0    0.0    0.0    0.0    0.0 ];  % quartz
cal.mnr_oxd = cal.mnr_oxd./sum(cal.mnr_oxd,2)*100;

% set mineral proportions in eutectic component
cal.mnrStr  = {'mnr_for','mnr_fay','mnr_px1','mnr_px2','mnr_px3','mnr_ant','mnr_alb','mnr_spn','mnr_qtz'};
for i = 1:length(cal.mnrStr); cal.(cal.mnrStr{i}) = i; end

cal.eut_mnr = [ 0.0   0.0   0.0   0.0   0.45  0.17  0.11  0.20  0.07];
cal.eut_mnr = cal.eut_mnr./sum(cal.eut_mnr,2);

% get oxide composition of eutectic component
cal.oxd = [cal.mnr_oxd([1 2 3 4 6],:); cal.eut_mnr*cal.mnr_oxd];

% specify calibration parameters

% set pure component melting points T_m^i at P=0
cal.T0(cal.for)  =  1890;
cal.T0(cal.fay)  =  1200;
cal.T0(cal.opx)  =  1425;
cal.T0(cal.cpx)  =  1166;
cal.T0(cal.ant)  =  1300;
cal.T0(cal.eut)  =   992;

% set first coeff. for P-dependence of T_m^i [GPa]
cal.A(cal.for)   =  22;
cal.A(cal.fay)   =   9;
cal.A(cal.opx)   =  16;
cal.A(cal.cpx)   =   5;
cal.A(cal.ant)   =  10;
cal.A(cal.eut)   =   4;

% set second coeff. for P-dependence of T_m^i [1]
cal.B(cal.for)   =  2.9;
cal.B(cal.fay)   =  2.2;
cal.B(cal.opx)   =  2.4;
cal.B(cal.cpx)   =  1.7;
cal.B(cal.ant)   =  2.4;
cal.B(cal.eut)   =  1.9;

% set entropy gain of fusion DeltaS [J/K]
cal.dS           =  300;

% set coeff. for T-dependence of partition coefficients K^i [1/K]
cal.r(cal.for)   =  30.5;
cal.r(cal.fay)   =  26.0;
cal.r(cal.opx)   =  69.9;
cal.r(cal.cpx)   =  15.0;
cal.r(cal.ant)   =  79.7;
cal.r(cal.eut)   =  10.8;

% specify geochemical model parameters
cal.nte          =  4;  % number of trace elements
cal.nir          =  2;  % number of isotope ratios
cal.Kte_cmp      =  [0.01  0.01  0.01  0.01  0.01  0.01  % component-wise trace element partition coefficients
                     0.10  0.10  0.10  0.10  0.10  0.10
                     3.00  3.00  3.00  3.00  3.00  3.00
                     10.0  10.0  10.0  10.0  10.0  10.0];