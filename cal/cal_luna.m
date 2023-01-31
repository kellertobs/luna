% set number of components
cal.nc  =  6;

% set string with component names
cal.cmpStr  =  {'for','fay','opx','cpx','ant','eut'};
for i = 1:cal.nc; cal.(cal.cmpStr{i}) = i; end

% set mineral end-member oxide compositions
cal.oxdStr =  {'SiO$_2$','Al$_2$O$_3$','FeO','MgO','CaO','Na$_2$O'};
cal.mnr_oxd = [  42.7    0.0    0.0   57.3    0.0      0     % forsterite
                 30.0    0.0   69.2    0.0    0.8      0     % fayalite
                 55.6    5.0    4.0   34.3    1.0    0.1     % Mg-pyroxene
                 50.5    8.2   12.2   24.5    4.5    0.1     % Mg-Fe-Ca-pyroxene
                 48.3    3.0   26.8    5.6   16.1    0.2     % Fe-Ca-pyroxene
                 44.2   35.6    0.4    0.2   19.1    0.5     % anorthite
                 67.4   20.4    1.5    0.1    0.6   10.0     % albite
                  0.5    4.0   93.6    1.5    0.4    0.0     % spinel
                100.0    0.0    0.0    0.0    0.0    0.0 ];  % quartz
cal.mnr_oxd = cal.mnr_oxd./sum(cal.mnr_oxd,2)*100;

% set mineral proportions in eutectic component
cal.mnrStr  = {'mnr_for','mnr_fay','mnr_px1','mnr_px2','mnr_px3','mnr_ant','mnr_alb','mnr_spn','mnr_qtz'};
for i = 1:length(cal.mnrStr); cal.(cal.mnrStr{i}) = i; end

cal.eut_mnr = [ 0.0   0.0   0.0   0.0   0.55  0.15  0.10  0.14  0.06];
cal.eut_mnr = cal.eut_mnr./sum(cal.eut_mnr,2);

% get oxide composition of eutectic component
cal.oxd = [cal.mnr_oxd([1 2 3 4 6],:); cal.eut_mnr*cal.mnr_oxd];

% set bulk composition
cal.c0 = [0.425    0.09    0.125    0.16    0.10    0.1];

% specify calibration parameters

% set pure component melting points T_m^i at P=0
cal.T0(cal.for)  =  1.9499e+03;
cal.T0(cal.fay)  =  1.1506e+03;
cal.T0(cal.opx)  =  1.5144e+03;
cal.T0(cal.cpx)  =  1.1026e+03;
cal.T0(cal.ant)  =  1.3009e+03;
cal.T0(cal.eut)  =  977.9228;

% set first coeff. for P-dependence of T_m^i [GPa]
cal.A(cal.for)   =  12.5794;
cal.A(cal.fay)   =   6.9059;
cal.A(cal.opx)   =  23.7030;
cal.A(cal.cpx)   =   6.6469;
cal.A(cal.ant)   =  24.8383;
cal.A(cal.eut)   =   7.9110;

% set second coeff. for P-dependence of T_m^i [1]
cal.B(cal.for)   =  3.1760;
cal.B(cal.fay)   =  1.0484;
cal.B(cal.opx)   =  2.9818;
cal.B(cal.cpx)   =  1.5006;
cal.B(cal.ant)   =  2.6789;
cal.B(cal.eut)   =  1.0041;

% set entropy gain of fusion DeltaS [J/K]
cal.dS           =  350;

% set coeff. for T-dependence of partition coefficients K^i [1/K]
cal.r(cal.for)   =  53.4197;
cal.r(cal.fay)   =  64.5492;
cal.r(cal.opx)   =  45.4267;
cal.r(cal.cpx)   =  15.1006;
cal.r(cal.ant)   =  49.3653;
cal.r(cal.eut)   =  24.6449;

% specify geochemical model parameters
cal.nte          =  4;  % number of trace elements
cal.nir          =  2;  % number of isotope ratios
cal.Kte_cmp      =  [0.01  0.01  0.01  0.01  0.01  0.01  % component-wise trace element partition coefficients
                     0.10  0.10  0.10  0.10  0.10  0.10
                     3.00  3.00  3.00  3.00  3.00  3.00
                     10.0  10.0  10.0  10.0  10.0  10.0];