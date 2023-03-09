% set number of components, mineral end-members, oxides
cal.ncmp  =  4;
cal.nmnr  =  10;
cal.noxd  =  7;

% set strings with component, mineral, oxide names
cal.cmpStr  =  {'dun','pxn','bas','eut'};
for i = 1:cal.ncmp; cal.(cal.cmpStr{i}) = i; end

cal.oxdStr =  {'SiO$_2$','TiO$_2$','Al$_2$O$_3$','FeO','MgO','CaO','Na$_2$O'};
for i = 1:cal.noxd; cal.(cal.oxdStr{i}(1:2)) = i; end

cal.mnrStr  = {'for','fay','ens','hyp','aug','pig','ant','alb','ilm','qtz'};
for i = 1:length(cal.mnrStr); cal.(cal.mnrStr{i}) = i; end

% set mineral end-member oxide compositions
cal.mnr_oxd = [  42.7    0.05   0.0    0.0   57.25   0.0    0.0     % forsterite
                 30.0    0.20   0.0   69.0    0.0    0.8    0.0     % fayalite
                 56.0    0.08   4.49   4.33	 33.92   1.11	0.07    % enstatite
                 50.0    0.15   8.85  12.17  26.08   2.70   0.05    % hypersthene
                 51.0    0.08   7.16   2.88  22.25  16.43   0.20    % Mg-augite
                 47.6    1.68   2.05  33.27   5.19  10.08   0.13    % Fe-pigeonite
                 44.2    0.08  35.52   0.4    0.2   19.1    0.5     % anorthite
                 67.4    0.20  20.4    1.1    0.3    2.1    8.5     % albite
                 0.32   20.56   3.07  74.63   1.14   0.28   0.0     % ilmenite
                99.95    0.05   0.0    0.0    0.0    0.0    0.0 ];  % quartz
cal.mnr_oxd = cal.mnr_oxd./sum(cal.mnr_oxd,2)*100;

% set mineral proportions in eutectic component
cal.cmp_mnr = [ 1.0000	0.0000	0.0000	0.0000	0.0000	0.0000	0.0000	0.0000	0.0000	0.0000      % dunite
                0.1600	0.1800	0.5400	0.1200	0.0000	0.0000	0.0000	0.0000	0.0000	0.0000      % pyroxenite
                0.0200	0.0700	0.0000	0.0500	0.2700	0.1000	0.3300	0.1000	0.0600	0.0000      % basalt
                0.0000	0.0000	0.0000	0.0000	0.0000	0.4500	0.1300	0.1400	0.1800	0.1000 ];   % eutectic
cal.cmp_mnr = cal.cmp_mnr./sum(cal.cmp_mnr,2);

% get oxide composition of eutectic component
cal.cmp_oxd = cal.cmp_mnr*cal.mnr_oxd;


% set pure component melting points T_m^i at P=0
cal.T0(cal.for)  =  1890;
cal.T0(cal.pxn)  =  1380;
cal.T0(cal.bas)  =  1090;
cal.T0(cal.eut)  =  1021;

% set first coeff. for P-dependence of T_m^i [GPa]
cal.A(cal.for)   =   6.1;
cal.A(cal.pxn)   =   4.7;
cal.A(cal.bas)   =   2.85;
cal.A(cal.eut)   =   2.7;

% set second coeff. for P-dependence of T_m^i [1]
cal.B(cal.for)   =  8.9;
cal.B(cal.pxn)   =  3.3;
cal.B(cal.bas)   =  2.5;
cal.B(cal.eut)   =  2.5;

% set entropy gain of fusion DeltaS [J/K]
cal.dS           =  300;

% set coeff. for T-dependence of partition coefficients K^i [1/K]
cal.r(cal.for)   =  22.0;
cal.r(cal.pxn)   =  18.0;
cal.r(cal.bas)   =   8.0;
cal.r(cal.eut)   =   6.0;

% specify geochemical model parameters
cal.nte          =  4;  % number of trace elements
cal.nir          =  2;  % number of isotope ratios
cal.Kte_cmp      =  [0.01,0.10,3.0,10.0].'.*ones(cal.nte,cal.nmnr); % mineral trace element partition coefficients
