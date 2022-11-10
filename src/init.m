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
TINY    =  1e-12;            % minimum cutoff phase, component fractions
BCA     =  {'',''};          % boundary condition on advection (top/bot, sides)

Dsx     = -cal.dS;           % use entropy change of crystallisation set in cal

% calculate dimensionless numbers characterising the system dynamics
var0.c = c0;          % in wt
var0.T = T0;          % convert to C
var0.P = Ptop/1e9;    % convert to GPa
var0.f = 1;           % in wt

[phs0,cal0]  =  meltmodel(var0,cal,'E');

m0  = phs0.f; x0 = 1-m0;
cm0 = phs0.cl;

% update phase oxide compositions
oxdm0 = cm0*cal.oxd;

wtm([1 3 4 6 7 8 11 12]) = [oxdm0,0,0]; % SiO2
eta0 = grdmodel08(wtm,T0);
rho0 = sum(cm0./rhom0).^-1;

DrhoT = aT*rho0*abs(Ttop-Tbot);
Drhoc = 0.1*max(abs(rhom0-rho0)); 
Drhox = 0.1*mean(abs(rhox0-rhom0));
Drho0 = DrhoT + Drhoc + Drhox;

uT    = DrhoT*g0*(D/10)^2/eta0/etareg;
uc    = Drhoc*g0*(D/10)^2/eta0/etareg;
ux    = Drhox*g0*(D/10)^2/eta0/etareg;
u0    = Drho0*g0*(D/10)^2/eta0/etareg;

wx0   = mean(abs(rhox0-rhom0))*g0*dx^2/eta0;

ud0   = kT0*dffreg/rho0/cP/(D/10);

uwT   = dw/tau_T; 

RaT   = uT/ud0;
Rac   = uc/ud0;
Rax   = ux/ud0;
Ra    = u0/ud0;

Rux   = wx0/u0;

RwT   = uwT/u0;

Re    = u0*rho0*(D/10)/eta0/etareg;
Rex   = wx0*rho0*dx/eta0;

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
[XX,ZZ]   = meshgrid(Xc,Zc);
Xf        = (Xc(1:end-1)+Xc(2:end))./2;
Zf        = (Zc(1:end-1)+Zc(2:end))./2;
[XXu,ZZu] = meshgrid(Xf,Zc);
[XXw,ZZw] = meshgrid(Xc,Zf);

Nx = length(Xc);
Nz = length(Zc);

inz = 2:Nz-1;
inx = 2:Nx-1;

% get smoothed initialisation field
rng(seed);
rp = randn(Nz,Nx);
for i = 1:round(smth)
    rp(2:end-1,2:end-1) = rp(2:end-1,2:end-1) + diff(rp(:,2:end-1),2,1)./8 + diff(rp(2:end-1,:),2,2)./8;
    rp = rp - mean(mean(rp(2:end-1,2:end-1)));
    rp([1 end],:) = 0;
    rp(:,[1 end]) = 0;
end
rp = rp./max(abs(rp(:)));

% get mapping arrays
NP =  Nz   * Nx   ;
NW = (Nz-1)* Nx   ;
NU =  Nz   *(Nx-1);
MapP = reshape(1:NP,Nz  ,Nx  );
MapW = reshape(1:NW,Nz-1,Nx  );
MapU = reshape(1:NU,Nz  ,Nx-1) + NW;

switch bndmode
    case 0  % none
        topshape = zeros(size(ZZ(inz,inx)));
        botshape = zeros(size(ZZ(inz,inx)));
    case 1  % top only
        topshape = exp( ( -ZZ(inz,inx)+h/2)/dw);
        botshape = zeros(size(ZZ(inz,inx)));
    case 2  % bot only
        topshape = zeros(size(ZZ(inz,inx)));
        botshape = exp(-(D-ZZ(inz,inx)-h/2)/dw);
    case 3  % top/bot only
        topshape = exp( ( -ZZ(inz,inx)+h/2)/dw);
        botshape = exp(-(D-ZZ(inz,inx)-h/2)/dw);
end
topshape = max(0,min(1,topshape));
botshape = max(0,min(1,botshape));

% set all boundaries to free slip
sds = -1;
top = -1;
bot = -1;
%                              sds = -1;      % free slip sides for other types
% if bndmode==1 || bndmode>=3; top = +1;      % no slip top for 'top only(1)', 'top/bot(3)', 'all sides(4)'
% else;                        top = -1; end  % free slip for other types
% if bndmode>=2;               bot = +1;      % no slip bot for 'bot only(2)', 'top/bot(3)', 'all sides(4)'
% else;                        bot = -1; end  % free slip for other types

% load calibration
run(['../cal/cal_',calID,'.m']);

% initialise solution fields
Tp  =  T0 + (T1-T0) .* (1+erf((ZZ/D-zlay)/wlay_T))/2 + dT.*rp;  if bndinit && ~isnan(Twall); Tp = Tp + (Twall-Tp).*bndshape; end % temperature [C]
c   =  zeros(Nz,Nx,cal.nc); cxq = c; cmq = c;
c0 = c0./sum(c0);  c1 = c1./sum(c1);
for i=1:cal.nc
    c(:,:,i)   =  c0(i) + (c1(i)-c0(i)) .* (1+erf((ZZ/D-zlay)/wlay_c))/2 + dc(i).*rp;  
    if bndinit && ~isnan(cwall); c(inz,inx,i) = c(inz,inx,i) + (cwall-c(inz,inx,i)).*bndshape; end % major component
end

te = zeros(Nz,Nx,length(te0));
for i = 1:length(te0)
    te(:,:,i)  =  te0(i) + (te1(i)-te0(i)) .* (1+erf((ZZ/D-zlay)/wlay_c))/2 + dte(i).*rp;  % trace elements
    if any(bndinit(:)) && ~isnan(tewall(i)); te(inz,inx,i)  = te(inz,inx,i) + (tewall(i)-te(inz,inx,i)).*bndshape; end 
end
ir = zeros(Nz,Nx,length(ir0));
for i = 1:length(ir0)
    ir(:,:,i)  =  ir0(i) + (ir1(i)-ir0(i)) .* (1+erf((ZZ/D-zlay)/wlay_c))/2 + dir(i).*rp;  % isotope ratios
    if any(bndinit(:)) && ~isnan(irwall(i)); ir(inz,inx,i)  = ir(inz,inx,i) + (irwall(i)-ir(inz,inx,i)).*bndshape; end 
end

U   =  zeros(size((XX(:,1:end-1)+XX(:,2:end))));  Ui = U;  res_U = 0.*U;
W   =  zeros(size((XX(1:end-1,:)+XX(2:end,:))));  Wi = W;  res_W = 0.*W; wf = 0.*W; wc = 0.*W;
P   =  zeros(size(XX));  Pi = P;  res_P = 0.*P;  meanQ = 0;  
SOL = [W(:);U(:);P(:)];

% initialise auxiliary fields
Wx  = W;  Ux  = U;
Wm  = W;  Um  = U;

eIIref =  1e-6;  
Div_V  =  0.*P;  Div_rhoV = 0.*P(inz,inx);  Div_rhoVo = Div_rhoV;
exx    =  0.*P;  ezz = 0.*P;  exz = zeros(Nz-1,Nx-1);  eII = 0.*P;  
txx    =  0.*P;  tzz = 0.*P;  txz = zeros(Nz-1,Nx-1);  tII = 0.*P; 
VolSrc =  0.*P(inz,inx);  MassErr = 0;  drhodt = 0.*P;  drhodto = 0.*P;

rho    =  rhom0(1).*ones(size(Tp));
rhoref =  mean(rho(inz,inx),'all');
Pt     =  Ptop + rhoref.*g0.*ZZ;
mq  =  ones(size(Tp));  m = mq;
xq  = zeros(size(Tp));  x = xq;

% get volume fractions and bulk density
step   =  0;
theta  = 1/2;
res = 1;  tol = 1e-12;
while res > tol
    xi = x;
    
    rhoref =  mean(rho(inz,inx),'all');
    Pt     =  Ptop + rhoref.*g0.*ZZ;
    Adbt   =  aT./rhoref;
    if Nz<=10; Pt = Ptop.*ones(size(Tp)); end

    T    =  (Tp+273.15).*exp(Adbt./cP.*Pt);

    var.c = reshape(c(inz,inx,:),(Nx-2)*(Nz-2),cal.nc);  % in wt
    var.T = reshape(T(inz,inx),(Nx-2)*(Nz-2),1)-273.15;  % convert to C
    var.P = reshape(Pt(inz,inx),(Nx-2)*(Nz-2),1)/1e9;    % convert to GPa
    var.f = 1-reshape(xq(inz,inx),(Nx-2)*(Nz-2),1);      % in wt
    
    [phs,cal]     =  meltmodel(var,cal,'E');
    
    mq(inz,inx)   = max(0,min(1,reshape(phs.f,(Nz-2),(Nx-2))));
    mq([1 end],:) = mq([2 end-1],:);
    mq(:,[1 end]) = mq(:,[2 end-1]);

    xq = 1-mq;
    x  = xq;  m = mq;

    for i = 2:cal.nc
        cxq(inz,inx,i) = reshape(phs.cs(:,i),(Nz-2),(Nx-2));
        cmq(inz,inx,i) = reshape(phs.cl(:,i),(Nz-2),(Nx-2));
    end
    cmq(:,:,1) = 1 - sum(cmq(:,:,2:end),3);
    cxq(:,:,1) = 1 - sum(cxq(:,:,2:end),3);

    cxq([1 end],:,:) = cxq([2 end-1],:,:);
    cxq(:,[1 end],:) = cxq(:,[2 end-1],:);
    cmq([1 end],:,:) = cmq([2 end-1],:,:);
    cmq(:,[1 end],:) = cmq(:,[2 end-1],:);

    cx = cxq;  cm = cmq;

    update;

    res  = norm(x(:)-xi(:),2)./sqrt(length(x(:)));
end
dto   = dt;
Pto   = Pt;
  
% get bulk enthalpy, silica, volatile content densities
X  = rho.*x;
S  = rho.*(cP.*log(T/(T0+273.15)) + x.*Dsx - Adbt.*(Pt-Ptop));  
S0 = rho.*(cP.*log((T0+273.15)) - Adbt.*(Ptop)); 
s  = S./rho;
C  = 0.*c;
for i = 1:cal.nc; C(:,:,i) = rho.*(m.*cm(:,:,i) + x.*cx(:,:,i)); end

% get phase entropies
sm = S./rho - x.*Dsx;
sx = sm + Dsx;

% get trace element phase compositions
for k = 1:length(te0)
    tem(:,:,k)  = te(:,:,k)./(m + x.*Kte(k) );
    tex(:,:,k)  = te(:,:,k)./(m./Kte(k)  + x);
end

% get geochemical component densities
for k = 1:length(te0)
    TE(:,:,k)  = rho.*(m.*tem(:,:,k) + x.*tex(:,:,k));
end
for k = 1:length(ir0)
    IR(:,:,k)  = rho.*ir(:,:,k);
end

% initialise reaction/decay rates
Gx = 0.*x(inz,inx); 

% initialise auxiliary variables 
dSdt   = 0.*T(inz,inx); bnd_S = 0.*T(inz,inx); diss_h = 0.*T(inz,inx); 
dCdt   = 0.*c(inz,inx);
dXdt   = 0.*x(inz,inx);
dTEdt  = 0.*te(inz,inx);
dIRdt  = 0.*ir(inz,inx);

% initialise timing and iterative parameters
time    =  0;
iter    =  0;
hist    = [];
dsumMdt = 0;
dsumSdt = 0;
dsumCdt = zeros(1,cal.nc);

% overwrite fields from file if restarting run
if restart
    if     restart < 0  % restart from last continuation frame
        name = [opdir,'/',runID,'/',runID,'_cont.mat'];
    elseif restart > 0  % restart from specified continuation frame
        name = [opdir,'/',runID,'/',runID,'_',num2str(restart),'.mat'];
    end
    if exist(name,'file')
        fprintf('\n   restart from %s \n\n',name);
        load(name,'U','W','P','Pt','x','m','chi','mu','X','S','C','T','c','cm','cx','TE','IR','te','ir','dSdt','dCdt','dXdt','dTEdt','dIRdt','Gx','rho','eta','eII','tII','dt','time','step','hist','VolSrc','wx','wm');
        
        xq = x;
        SOL = [W(:);U(:);P(:)];
        rhoo = rho; Div_rhoVo = Div_rhoV;
        update; output;
    else % continuation file does not exist, start from scratch
        fprintf('\n   !!! restart file does not exist !!! \n   => starting run from scratch %s \n\n',name);
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