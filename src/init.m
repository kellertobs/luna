% create output directory
if ~isfolder([opdir,'/',runID])
    mkdir([opdir,'/',runID]);
end

% save input parameters and runtime options (unless restarting)
if restart == 0 
    parfile = [opdir,'/',runID,'/',runID,'_par'];
    save(parfile);
end

fprintf('\n\n')
fprintf('*************************************************************\n');
fprintf('*****  RUN LUNA MODEL | %s  ***************\n',datetime('now'));
fprintf('*************************************************************\n');
fprintf('\n   run ID: %s \n\n',runID);

load ocean;                  % load custom colormap
run(['../cal/cal_',calID]);  % load melt model calibration
calibrt =  0;                % not in calibrate mode
TINY    =  1e-16;            % minimum cutoff phase, component fractions
BCA     =  {'','periodic'};  % boundary condition on advection (top/bot, sides)
BCD     =  {'','periodic'};  % boundary condition on advection (top/bot, sides)
bnchm   =  0;                % not a benchmark run

Dsx     = -cal.dS;           % use entropy change of crystallisation set in cal

% calculate dimensionless numbers characterising the system dynamics
wtm    = [c0*cal.cmp_oxd,0,0]; % 8 major elements + H2O
rhom0  = DensityX(wtm,T0,Ptop/1e8);
var0.c = c0;          % in wt
var0.T = T0.*exp(aT./rhom0./cP.*Ptop);  % in C
var0.P = Ptop/1e9;    % convert to GPa
var0.f = 1;           % in wt

[phs0,cal0]  =  meltmodel(var0,cal,'E');

m0  = phs0.f; x0 = 1-m0;
cm0 = phs0.cl;
cx0 = phs0.cs;

% update phase oxide compositions
cm0_oxd = cm0*cal.cmp_oxd;
cx0_oxd = cx0*cal.cmp_oxd;
cm0_mnr = cm0*cal.cmp_mnr;
cx0_mnr = cx0*cal.cmp_mnr;

cm1_oxd = (0.99.*cm0_oxd + 0.01.*cal.cmp_oxd(1,:));
cm2_oxd = (0.99.*cm0_oxd + 0.01.*cal.cmp_oxd(4,:));

wtm    = [cm0_oxd,0,0]; % 8 major elements + H2O
etam0  = Giordano08(wtm,T0);
rhom0  = DensityX(wtm,T0,Ptop/1e8);
Drhox0 = abs(sum(cx0_mnr./rhox0).^-1 - rhom0);

wtm    = [cm1_oxd,0,0];
rhom1  = DensityX(wtm,T0,Ptop/1e8);

wtm    = [cm2_oxd,0,0];
rhom2  = DensityX(wtm,T0,Ptop/1e8);

DrhoT  = rhom0.*aT*max([abs(Ttop-Tbot)/100,abs(T0-T1),T0/100]);
Drhoc  = abs(rhom1-rhom2);
Drhox  = 0.01*Drhox0;
Drho0  = DrhoT + Drhoc + Drhox;

uT    = DrhoT*g0*(D/10)^2/etam0/etareg;
uc    = Drhoc*g0*(D/10)^2/etam0/etareg;
ux    = Drhox*g0*(D/10)^2/etam0/etareg;
u0    = Drho0*g0*(D/10)^2/etam0/etareg;

wx0   = Drhox0*g0*d0^2/etam0;

kW0   = u0/10*h/10;
ud0   = (kT0+rhom0.*cP.*(kW0 + mink))/rhom0/cP/(D/10);

uwT   = bnd_w/tau_T; 

RaT   = uT/ud0;
Rac   = uc/ud0;
Rax   = ux/ud0;
Ra    = u0/ud0;

Rux   = wx0/u0;

RwT   = uwT/u0;

Re    = u0*rhom0*(D/10)/etam0/etareg;
Rex   = wx0*rhom0*d0/etam0;

fprintf('    crystal Re: %1.3e \n'  ,Rex);
fprintf('     system Re: %1.3e \n\n',Re );

fprintf('    thermal Ra: %1.3e \n'  ,RaT);
fprintf('   chemical Ra: %1.3e \n'  ,Rac);
fprintf('    crystal Ra: %1.3e \n'  ,Rax);
fprintf('   combined Ra: %1.3e \n\n',Ra );

fprintf('    crystal Ru: %1.3e \n\n',Rux);

fprintf('    thermal Rw: %1.3e \n\n',RwT);


% get coordinate arrays
Xc        = -h/2:h:L+h/2;
Zc        = -h/2:h:D+h/2;
Xu        = (Xc(1:end-1)+Xc(2:end))./2;
Zw        = (Zc(1:end-1)+Zc(2:end))./2;
[XXu,ZZu] = meshgrid(Xu,Zc);
[XXw,ZZw] = meshgrid(Xc,Zw);
Xc        = Xc(2:end-1);
Zc        = Zc(2:end-1);
[XX,ZZ]   = meshgrid(Xc,Zc);

Nx = length(Xc);
Nz = length(Zc);

inz = 2:Nz-1;
inx = 2:Nx-1;

% get smoothed initialisation field
rng(seed);
smth = smth*Nx*Nz*1e-4;
rp   = randn(Nz,Nx);
for i = 1:round(smth)
    rp = rp + diffus(rp,1/8*ones(size(rp)),1,[1,2],BCD);
    rp = rp - mean(mean(rp));
end
rp = rp./max(abs(rp(:)));


% get mapping arrays
NP = (Nz+2) * (Nx+2);
NW = (Nz+1) * (Nx+2);
NU = (Nz+2) * (Nx+1);
MapP = reshape(1:NP,Nz+2,Nx+2);
MapW = reshape(1:NW,Nz+1,Nx+2);
MapU = reshape(1:NU,Nz+2,Nx+1) + NW;

% set up shape functions for transient boundary layers
topshape = zeros(size(ZZ));
botshape = zeros(size(ZZ));
switch bndmode
    case 0  % none
    case 1  % top only
        topshape = exp( ( -ZZ+h/2)/bnd_w);
    case 2  % bot only
        botshape = exp(-(D-ZZ-h/2)/bnd_w);
    case 3  % top/bot only
        topshape = exp( ( -ZZ+h/2)/bnd_w);
        botshape = exp(-(D-ZZ-h/2)/bnd_w);
end

% set velocity boundaries
sds =  0;  % hindered slip
top = -1;  % free slip
bot = -1;  % free slip

% load calibration
run(['../cal/cal_',calID,'.m']);

% initialise solution fields
Tp  =  T0 + (T1-T0) .* (1+erf((ZZ/D-zlay)/wlay_T))/2 + dT.*rp;  if bndinit && ~isnan(Twall); Tp = Tp + (Twall-Tp).*bndshape; end % temperature [C]
c   =  zeros(Nz,Nx,cal.ncmp); cxq = c; cmq = c;
c0 = c0./sum(c0);  c1 = c1./sum(c1);
for i=1:cal.ncmp
    c(:,:,i)   =  c0(i) + (c1(i)-c0(i)) .* (1+erf((ZZ/D-zlay)/wlay_c))/2 + dc(i).*rp;  
    if bndinit && ~isnan(cwall); c(inz,inx,i) = c(inz,inx,i) + (cwall-c(inz,inx,i)).*bndshape; end % major component
end

te = zeros(Nz,Nx,cal.nte);
for i = 1:length(te0)
    te(:,:,i)  =  te0(i) + (te1(i)-te0(i)) .* (1+erf((ZZ/D-zlay)/wlay_c))/2 + dte(i).*rp;  % trace elements
    if any(bndinit(:)) && ~isnan(tewall(i)); te(inz,inx,i)  = te(inz,inx,i) + (tewall(i)-te(inz,inx,i)).*bndshape; end 
end
ir = zeros(Nz,Nx,cal.nir);
for i = 1:length(ir0)
    ir(:,:,i)  =  ir0(i) + (ir1(i)-ir0(i)) .* (1+erf((ZZ/D-zlay)/wlay_c))/2 + dir(i).*rp;  % isotope ratios
    if any(bndinit(:)) && ~isnan(irwall(i)); ir(inz,inx,i)  = ir(inz,inx,i) + (irwall(i)-ir(inz,inx,i)).*bndshape; end 
end

U   =  zeros(Nz+2,Nx+1);  UBG = U; Ui = U;
W   =  zeros(Nz+1,Nx+2);  WBG = W; Wi = W; wx = 0.*W; wm = 0.*W;
P   =  zeros(Nz+2,Nx+2);  Vel = 0.*Tp;
SOL = [W(:);U(:);P(:)];

% initialise auxiliary fields
Wx  = W;  Ux  = U;
Wm  = W;  Um  = U;

eIIref =  1e-6;  
Div_V  =  0.*Tp;  advn_rho = 0.*Tp;  drhodt = 0.*Tp;  drhodto = drhodt;
exx    =  0.*Tp;  ezz = 0.*Tp;  exz = zeros(Nz-1,Nx-1);  eII = 0.*Tp;  
txx    =  0.*Tp;  tzz = 0.*Tp;  txz = zeros(Nz-1,Nx-1);  tII = 0.*Tp; 
VolSrc =  0.*Tp;  MassErr = 0;

rho    =  rhom0(1).*ones(size(Tp));
rhoref =  mean(rho,'all');
Pt     = Ptop + rhoref.*g0.*ZZ;
mq  =  ones(size(Tp));  m = mq;
xq  = zeros(size(Tp));  x = xq;

% get volume fractions and bulk density
step   = 0;
res    = 1;  tol = 1e-12;
cal.Tliq = reshape(Tp,[],1);
while res > tol
    Pti = Pt;
    
    rhoref = mean(rho,'all');
    Adbt   = aT./rhoref;
    if Nz==1; Pt = Ptop.*ones(size(Tp)); else
        rhofz  = (rho(1:end-1,:)+rho(2:end,:))/2;
        Pt(1,:)     = repmat(mean(rhofz(1,:),2).*g0.*h/2,1,Nx) + Ptop;
        Pt(2:end,:) = Pt(1,:) + repmat(cumsum(mean(rhofz,2).*g0.*h),1,Nx);
    end

    T  = (Tp+273.15).*exp(Adbt./cP.*Pt);

%     T   = max(T,reshape(cal.Tliq+273.15,Nz,Nx));

    var.c = reshape(c,Nx*Nz,cal.ncmp);  % in wt
    var.T = reshape(T,Nx*Nz,1)-273.15;  % convert to C
    var.P = reshape(Pt,Nx*Nz,1)/1e9;    % convert to GPa
    var.f = 1-reshape(xq,Nx*Nz,1);      % in wt
    
    [phs,cal] =  meltmodel(var,cal,'E');
    
    mq   = max(0,min(1,reshape(phs.f,Nz,Nx)));

    xq = 1-mq;
    x  = xq;  m = mq;

    cxq = reshape(phs.cs,Nz,Nx,cal.ncmp);
    cmq = reshape(phs.cl,Nz,Nx,cal.ncmp);

    cmq = cmq./sum(cmq,3);
    cxq = cxq./sum(cxq,3);

    cx = cxq;  cm = cmq;

    update;

    res  = norm(Pt(:)-Pti(:),2)./norm(Pt(:),2);
end
rhoo = rho; 
dto   = dt;
Pto   = Pt;
  
% get bulk enthalpy, silica, volatile content densities
X  = rho.*x;  Xo = X;  M = rho-X;  Mo = M;  RHO = M+X;
S  = rho.*(cP.*log(T/(T0+273.15)) + x.*Dsx - Adbt.*(Pt-Ptop));  So = S;
S0 = rho.*(cP.*log((T0+273.15)) - Adbt.*(Ptop)); 
s  = S./rho;
C  = 0.*c;  Co = C;
for i = 1:cal.ncmp; C(:,:,i) = rho.*(m.*cm(:,:,i) + x.*cx(:,:,i)); end

% get phase entropies
sm = S./rho - x.*Dsx;
sx = sm + Dsx;

% get trace element phase compositions
for k = 1:cal.nte
    tem(:,:,k)  = te(:,:,k)./(m + x.*Kte(k) );
    tex(:,:,k)  = te(:,:,k)./(m./Kte(k)  + x);
end

% get geochemical component densities
for k = 1:cal.nte
    TE(:,:,k)  = rho.*(m.*tem(:,:,k) + x.*tex(:,:,k));
end
TEo = TE;

for k = 1:cal.nir
    IR(:,:,k)  = rho.*ir(:,:,k);
end
IRo = IR;

% initialise reaction/decay rates
Gx = 0.*x;  Gm = 0.*m;

% initialise auxiliary variables 
dSdt   = 0.*T; dSdto = dSdt; bnd_S = 0.*T; diss_h = 0.*T; 
dCdt   = 0.*c; dCdto = dCdt;
dXdt   = 0.*x; dXdto = dXdt;
dMdt   = 0.*m; dMdto = dMdt;
dTEdt  = 0.*te; dTEdto = dTEdt;
dIRdt  = 0.*ir; dIRdto = dIRdt;

% initialise timing and iterative parameters
time    =  0;
iter    =  0;
hist    = [];
dsumMdt = 0; dsumMdto = 0;
dsumSdt = 0; dsumSdto = 0;
dsumCdt = zeros(1,cal.ncmp); dsumCdto = zeros(1,cal.ncmp);

% overwrite fields from file if restarting run
if restart
    if     restart < 0  % restart from last continuation frame
        name = [opdir,'/',runID,'/',runID,'_cont.mat'];
    elseif restart > 0  % restart from specified continuation frame
        name = [opdir,'/',runID,'/',runID,'_',num2str(restart),'.mat'];
    end
    if exist(name,'file')
        fprintf('\n   restart from %s \n\n',name);
        load(name,'U','W','P','Pt','x','m','chi','mu','X','S','C','T','c','cm','cx','TE','IR','te','ir','dSdt','dCdt','dXdt','dTEdt','dIRdt','Gx','rho','eta','eII','tII','dt','time','step','VolSrc','wx','wm');
        name = [opdir,'/',runID,'/',runID,'_hist'];
        load(name,'hist');

        xq  = x;
        SOL = [W(:);U(:);P(:)];

        So = S;
        Co = C;
        Xo = X;
        Mo = M;
        rhoo = rho;
        TEo = TE;
        IRo = IR;
        dSdto = dSdt;
        dCdto = dCdt;
        dXdto = dXdt;
        dMdto = dMdt;
        drhodto = drhodt;
        dTEdto = dTEdt;
        dIRdto = dIRdt;
        Div_Vo = Div_V;
        dto    = dt;

        update; output; restart = 0;

    else % continuation file does not exist, start from scratch
        fprintf('\n   !!! restart file does not exist !!! \n   => starting run from scratch %s \n\n',name);
        update;
        fluidmech;
        history;
        output;
    end
else
    % complete, plot, and save initial condition
    update;
    fluidmech;
    history;
    output;
end

time  = time+dt;
step  = step+1;