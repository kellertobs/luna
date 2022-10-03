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

% calculate dimensionless numbers characterising the system dynamics
var0.c = c0;          % in wt
var0.T = T0;          % convert to C
var0.P = Ptop/1e9;    % convert to GPa
var0.f = 1;           % in wt

[phs0,cal0]  =  meltmodel(var0,cal,'E');

m0  = phs0.f; x0 = 1-m0;
cm0 = phs0.cl;

% update phase oxide compositions
oxdm0 = cm0*cal.oxds;

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
        topshape = zeros(size(ZZ));
        botshape = zeros(size(ZZ));
    case 1  % top only
        topshape = exp( ( -ZZ+h/2)/dw);
        botshape = zeros(size(ZZ));
    case 2  % bot only
        topshape = zeros(size(ZZ));
        botshape = exp(-(D-ZZ-h/2)/dw);
    case 3  % top/bot only
        topshape = exp( ( -ZZ+h/2)/dw);
        botshape = exp(-(D-ZZ-h/2)/dw);
end
topshape = max(0,min(1,topshape));
topshape([1 end],:) = topshape([2 end-1],:);
topshape(:,[1 end]) = topshape(:,[2 end-1]);

botshape = max(0,min(1,botshape));
botshape([1 end],:) = botshape([2 end-1],:);
botshape(:,[1 end]) = botshape(:,[2 end-1]);

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
c   =  zeros(cal.nc,Nz,Nx); cxq = c; cmq = c;
for i=1:cal.nc
    c(i,:,:)   =  c0(i) + (cl(i)-c0(i)) .* (1+erf((ZZ/D-zlay)/wlay_c))/2 + dc(i).*rp;  if bndinit && ~isnan(cwall); c(i,:,:) = c(i,:,:) + (cwall-c(i,:,:)).*bndshape; end % major component
end

it  =  it0 + (it1-it0) .* (1+erf((ZZ/D-zlay)/wlay_c))/2 + dit.*rp;  if bndinit && ~isnan(itwall); it  = it  + (itwall-it ).*bndshape; end % incompatible trace element
ct  =  ct0 + (ct1-ct0) .* (1+erf((ZZ/D-zlay)/wlay_c))/2 + dct.*rp;  if bndinit && ~isnan(ctwall); ct  = ct  + (ctwall-ct ).*bndshape; end % compatible trace element
si  =  si0 + (si1-si0) .* (1+erf((ZZ/D-zlay)/wlay_c))/2 + dsi.*rp;  if bndinit && ~isnan(siwall); si  = si  + (siwall-si ).*bndshape; end % stable isotope ratio
rip =  ri0 + (ri1-ri0) .* (1+erf((ZZ/D-zlay)/wlay_c))/2 + dri.*rp;  if bndinit && ~isnan(riwall); rip = rip + (riwall-rip).*bndshape; end % radiogenic isotope parent
rid =  rip.*HLRID./HLRIP;                                           % radiogenic isotope daughter

U   =  zeros(size((XX(:,1:end-1)+XX(:,2:end))));  Ui = U;  res_U = 0.*U;
W   =  zeros(size((XX(1:end-1,:)+XX(2:end,:))));  Wi = W;  res_W = 0.*W; wf = 0.*W; wc = 0.*W;
P   =  zeros(size(XX));  Pi = P;  res_P = 0.*P;  meanQ = 0;  
SS  = [W(:);U(:);P(:)];

% initialise auxiliary fields
Wx  = W;  Ux  = U;
Wm  = W;  Um  = U;

eIIref =  1e-6;  
Div_V  =  0.*P;  Div_rhoV = 0.*P;  Div_rhoVo = Div_rhoV;
exx    =  0.*P;  ezz = 0.*P;  exz = zeros(Nz-1,Nx-1);  eII = 0.*P;  
txx    =  0.*P;  tzz = 0.*P;  txz = zeros(Nz-1,Nx-1);  tII = 0.*P; 
VolSrc =  0.*P;  MassErr = 0;  drhodt = 0.*P;  drhodto = 0.*P;

rhoref = rhom0(1);
rhoo   = rhoref.*ones(size(P));
Pt     = rhoref.*g0.*ZZ + Ptop;  
if Nz<=10; Pt = mean(mean(Pt(2:end-1,2:end-1))).*ones(size(Pt)); end
T   =  (Tp+273.15).*exp(aT./rhoref./cP.*Pt); % convert to [K]

% get volume fractions and bulk density
step   =  0;
theta  =  1.0;
res = 1;  tol = 1e-15;  x = ones(size(T))./10;
while res > tol
    xi = x;
    
    for i = 1:cal.nc
        var.c(:,i) = reshape(c(i,:,:),Nx*Nz,1);  % in wt
    end
    var.T = T(:)-273.15;  % convert to C
    var.P = Pt(:)/1e9;    % convert to GPa
    var.f = 1-x(:);       % in wt
    
    [phs,cal]  =  meltmodel(var,cal,'E');
    
    mq  = max(TINY,min(1-TINY,reshape(phs.f ,Nz,Nx))); xq = 1-mq; x = xq; m = mq;
    for i = 1:cal.nc
        cxq(i,:,:) = reshape(phs.cs(:,i),Nz,Nx); cx = cxq;
        cmq(i,:,:) = reshape(phs.cl(:,i),Nz,Nx); cm = cmq;
    end

    update;
    
    rhoref  = mean(mean(rho(2:end-1,2:end-1)));
    Pt      = Ptop + rhoref.*g0.*ZZ;
    if Nz<=10; Pt = mean(mean(Pt(2:end-1,2:end-1))); end
    rhoo  = rho;

    T    =  (Tp+273.15).*exp(aT./rhoref./cP.*Pt);

    res  = norm(x(:)-xi(:),2)./sqrt(length(x(:)));
end
dto   = dt;
Pto   = Pt;

% get geochemical phase compositions
itm  = it./(m + x.*KIT); itx = it./(m./KIT + x);
ctm  = ct./(m + x.*KCT); ctx = ct./(m./KCT + x);
ripm = rip./(m + x.*KRIP); ripx = rip./(m./KRIP + x);
ridm = rid./(m + x.*KRID); ridx = rid./(m./KRID + x);
  
% get bulk enthalpy, silica, volatile content densities
X  = rho.*x;
S  = rho.*(cP.*log(T/(T0+273.15)) + x.*Dsx - aT./rhoref.*(Pt-Ptop));  
S0 = rho.*(cP.*log((T0+273.15)) - aT./rhoref.*(Ptop)); 
s  = S./rho;
C  = 0.*c;
for i = 1:cal.nc; C(i,:,:) = rho.*(m.*squeeze(cm(i,:,:)) + x.*squeeze(cx(i,:,:))); end

sumC = squeeze(sum(C));
for i = 1:cal.nc; c(i,:,:) = squeeze(C(i,:,:))./sumC; end
for i = 1:cal.nc; C(i,:,:) = squeeze(c(i,:,:)).*rho; end

% get phase entropies
sm = S./rho - x.*Dsx;
sx = sm + Dsx;

% get geochemical content densities
IT  = rho.*(m.*itm + x.*itx);
CT  = rho.*(m.*ctm + x.*ctx);
SI  = rho.*(m.*si  + x.*si);
RIP = rho.*(m.*ripm + x.*ripx);
RID = rho.*(m.*ridm + x.*ridx);

% initialise reaction/decay rates
Gx = 0.*x; 
dcy_rip = rho.*rip./HLRIP.*log(2);
dcy_rid = rho.*rid./HLRID.*log(2);

% initialise auxiliary variables 
dTdt   = 0.*T;  diff_T = 0.*T; diss_h = 0.*T; bndS = zeros(size(ZZ));
dSdt   = 0.*T;
dCdt   = 0.*c;  diff_c = 0.*T;
dXdt   = 0.*x;  diff_x = 0.*x;  
dSIdt  = 0.*si;  diff_si  = 0.*SI;
dITdt  = 0.*IT;  diff_it  = 0.*IT;
dCTdt  = 0.*CT;  diff_ct  = 0.*CT;
dRIPdt = 0.*RIP; diff_rip  = 0.*RIP;
dRIDdt = 0.*RID; diff_rid  = 0.*RID;

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
        load(name,'U','W','P','Pt','x','m','chi','mu','S','C','T','X','c','cm','cx','IT','CT','SI','RIP','RID','it','ct','si','rip','rid','dSdt','dCdt','dXdt','dITdt','dCTdt','dSIdt','dRIDdt','dRIPdt','Gx','rho','eta','exx','ezz','exz','txx','tzz','txz','eII','tII','dt','time','step','hist','VolSrc','wx','wm');

        xq = x;
        SOL = [W(:);U(:);P(:)];
        dcy_rip = rho.*rip./HLRIP.*log(2);
        dcy_rid = rho.*rid./HLRID.*log(2);
        Pto = Pt; etao = eta; rhoo = rho; Div_rhoVo = Div_rhoV;
        update; output;
        time  = time+dt;
        step  = step+1;
        theta = 0.5;
    else % continuation file does not exist, start from scratch
        fprintf('\n   !!! restart file does not exist !!! \n   => starting run from scratch %s \n\n',name);
    end
else
    % complete, plot, and save initial condition
    update;
    fluidmech;
    history;
    output;
    step = step+1;
end