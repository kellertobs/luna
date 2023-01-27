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
cal.T0(cal.for)  =  1900;
cal.T0(cal.fay)  =  1200;
cal.T0(cal.opx)  =  1550;
cal.T0(cal.cpx)  =  1150;
cal.T0(cal.ant)  =  1350;
cal.T0(cal.eut)  =   950;

% set first coeff. for P-dependence of T_m^i [GPa]
cal.A(cal.for)   =  10.0;
cal.A(cal.fay)   =  10.0;
cal.A(cal.opx)   =  15.0;
cal.A(cal.cpx)   =  10.0;
cal.A(cal.ant)   =  15.0;
cal.A(cal.eut)   =  10.0;

% set second coeff. for P-dependence of T_m^i [1]
cal.B(cal.for)   =  3.50;
cal.B(cal.fay)   =  1.50;
cal.B(cal.opx)   =  2.50;
cal.B(cal.cpx)   =  2.00;
cal.B(cal.ant)   =  2.20;
cal.B(cal.eut)   =  1.50;

% set entropy gain of fusion DeltaS [J/K]
cal.dS           =  350;

% set coeff. for T-dependence of partition coefficients K^i [1/K]
cal.r(cal.for)   =  55;
cal.r(cal.fay)   =  55;
cal.r(cal.opx)   =  55;
cal.r(cal.cpx)   =  25;
cal.r(cal.ant)   =  55;
cal.r(cal.eut)   =  25;


% set up calibration parameter bounds for mcmc
bnds.T0(cal.for,:) = cal.T0(cal.for) + [-1,+1]*50;
bnds.T0(cal.fay,:) = cal.T0(cal.fay) + [-1,+1]*50;
bnds.T0(cal.opx,:) = cal.T0(cal.opx) + [-1,+1]*50;
bnds.T0(cal.cpx,:) = cal.T0(cal.cpx) + [-1,+1]*50;
bnds.T0(cal.ant,:) = cal.T0(cal.ant) + [-1,+1]*50;
bnds.T0(cal.eut,:) = cal.T0(cal.eut) + [-1,+1]*50;

bnds.A(cal.for,:)  = cal.A(cal.for)  + [-1,+1]*10;
bnds.A(cal.fay,:)  = cal.A(cal.fay)  + [-1,+1]*10;
bnds.A(cal.opx,:)  = cal.A(cal.opx)  + [-1,+1]*10;
bnds.A(cal.cpx,:)  = cal.A(cal.cpx)  + [-1,+1]*10;
bnds.A(cal.ant,:)  = cal.A(cal.ant)  + [-1,+1]*10;
bnds.A(cal.eut,:)  = cal.A(cal.eut)  + [-1,+1]*10;

bnds.B(cal.for,:)  = cal.B(cal.for)  + [-1,+1]*0.5;
bnds.B(cal.fay,:)  = cal.B(cal.fay)  + [-1,+1]*0.5;
bnds.B(cal.opx,:)  = cal.B(cal.opx)  + [-1,+1]*0.5;
bnds.B(cal.cpx,:)  = cal.B(cal.cpx)  + [-1,+1]*0.5;
bnds.B(cal.ant,:)  = cal.B(cal.ant)  + [-1,+1]*0.5;
bnds.B(cal.eut,:)  = cal.B(cal.eut)  + [-1,+1]*0.5;

bnds.r(cal.for,:)  = cal.r(cal.for)  + [-1,+1]*10;
bnds.r(cal.fay,:)  = cal.r(cal.fay)  + [-1,+1]*10;
bnds.r(cal.opx,:)  = cal.r(cal.opx)  + [-1,+1]*10;
bnds.r(cal.cpx,:)  = cal.r(cal.cpx)  + [-1,+1]*10;
bnds.r(cal.ant,:)  = cal.r(cal.ant)  + [-1,+1]*10;
bnds.r(cal.eut,:)  = cal.r(cal.eut)  + [-1,+1]*10;

bnds.mat           = [bnds.T0; bnds.A; bnds.B; bnds.r];

% variable names
vname = strcat(repelem({'T0';'A';'B';'r'},cal.nc,1), ',', repmat(cal.cmpStr',4,1));
