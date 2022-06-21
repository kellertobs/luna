load ocean;  % load custom colormap

% minimum cutoff phase, component fractions
TINY     =  1e-16;               

% get coordinate arrays
X         = -h/2:h:L+h/2;
Z         = -h/2:h:D+h/2;
[XX,ZZ]   = meshgrid(X,Z);
Xfc       = (X(1:end-1)+X(2:end))./2;
Zfc       = (Z(1:end-1)+Z(2:end))./2;
[XXu,ZZu] = meshgrid(Xfc,Z);
[XXw,ZZw] = meshgrid(X,Zfc);

Nx = length(X);
Nz = length(Z);

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
Wf  = W;  Uf  = U; 
Wx  = W;  Ux  = U;
Wm  = W;  Um  = U;

eIIref =  1e-6;  
Div_V  =  0.*P;  Div_rhoV = 0.*P;  Div_rhoVo = Div_rhoV;
exx    =  0.*P;  ezz = 0.*P;  exz = zeros(Nz-1,Nx-1);  eII = 0.*P;  
txx    =  0.*P;  tzz = 0.*P;  txz = zeros(Nz-1,Nx-1);  tII = 0.*P; 
VolSrc =  0.*P;  MassErr = 0;  drhodt = 0.*P;  drhodto = 0.*P;

if ~react;  Dsx = 0;  Dsf = 0;  end
rhoref = rhom0(1);
rhoo   = rhoref.*ones(size(P));
Pt     = rhoref.*g0.*ZZ + Ptop;  
if Nz<=10; Pt = mean(mean(Pt(2:end-1,2:end-1))).*ones(size(Pt)); end
T   =  (Tp+273.15).*exp(aT./rhoref./cP.*Pt); % convert to [K]

% get volume fractions and bulk density
ALPHA  =  1.0;
THETA  =  1.0;
res = 1;  tol = 1e-15;  x = ones(size(T))./10;
while res > tol
    xi = x;
    
    for i = 1:cal.nc
        var.c(:,i) = reshape(c(i,:,:),Nx*Nz,1);  % in wt
    end
    var.T = T(:)-273.15;  % convert to C
    var.P = Pt(:)/1e9;    % convert to GPa
    var.f = 1-x(:);       % in wt
    
    [phs,cal]  =  meltmodel(var,cal,'E'); % cs and cl component prop in each phase VS T 
    
    mq  = reshape(phs.f ,Nz,Nx); xq = 1-mq; x = xq; m = mq;
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
S = rho.*(cP.*log(T/(T0+273.15)) + x.*Dsx - aT./rhoref.*Pt);
C = 0.*c;
for i = 1:cal.nc; C(i,:,:) = rho.*(m.*squeeze(cm(i,:,:)) + x.*squeeze(cx(i,:,:))); end

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
dxdt   = 0.*x;  diff_x = 0.*x;  
dSIdt  = 0.*si;  diff_si  = 0.*SI;
dITdt  = 0.*IT;  diff_it  = 0.*IT;
dCTdt  = 0.*CT;  diff_ct  = 0.*CT;
dRIPdt = 0.*RIP; diff_rip  = 0.*RIP;
dRIDdt = 0.*RID; diff_rid  = 0.*RID;

% initialise timing and iterative parameters
step    =  0;
time    =  0;
iter    =  0;
hist    = [];
dsumMdt = 0;
dsumSdt = 0;
dsumCdt = zeros(1,cal.nc);

% overwrite fields from file if restarting run
if restart
    if     restart < 0  % restart from last continuation frame
        name = [opdir,'/',runID,'/',runID,'_cont'];
    elseif restart > 0  % restart from specified continuation frame
        name = [opdir,'/',runID,'/',runID,'_',num2str(restart)];
    end
    load(name,'U','W','P','Pt','x','m','chi','mu','S','C','T','c','cm','cx','IT','CT','SI','RIP','RID','it','ct','si','rip','rid','dTdt','dCdt','dITdt','dCTdt','dSIdt','dxdt','Gx','rho','eta','exx','ezz','exz','txx','tzz','txz','eII','tII','dt','time','step','hist','VolSrc','wx','wm');
    
    xq = x; 
    SS = [W(:);U(:);P(:)];
    dcy_rip = rho.*rip./HLRIP.*log(2);
    dcy_rid = rho.*rid./HLRID.*log(2);
    Pto = Pt; etao = eta; rhoo = rho; Div_rhoVo = Div_rhoV;
    update; output;
    time  = time+dt;
    step  = step+1;
else
    % update coefficients and run initial fluidmech solve
    if ~bnchm
        update; 
        fluidmech; 
        history; 
        output;
    end
    time  = time+dt;
    step  = step+1;
end