%% *****  THERMO-CHEMICAL EVOLUTION  **************************************

% store previous iterationcal.ncmp
Si = S; Ci = C; Xi = X;

%***  update heat content (entropy) density

% heat advection
advn_S = - advect(M.*sm,Um(2:end-1,:),Wm(:,2:end-1),h,{ADVN,''},[1,2],BCA) ...  % melt  advection
         - advect(X.*sx,Ux(2:end-1,:),Wx(:,2:end-1),h,{ADVN,''},[1,2],BCA);     % fluid advection

diff_S = diffus(T,ks,h,[1,2],BCD);

% heat dissipation
diss_h = diss ./ T;

% boundary layers
bnd_T = zeros(size(T));
if ~isnan(Ttop); bnd_T = bnd_T + ((Ttop+273.15)-T)./tau_T .* topshape; end % impose top boundary layer
if ~isnan(Tbot); bnd_T = bnd_T + ((Tbot+273.15)-T)./tau_T .* botshape; end % impose bot boundary layer
bnd_S = rho.*cP.*bnd_T./T;

% total rate of change
dSdt  = advn_S + diff_S + diss_h + bnd_S;

% residual of entropy evolution
res_S = (a1*S-a2*So-a3*Soo)/dt - (b1*dSdt + b2*dSdto + b3*dSdtoo);

% semi-implicit update of bulk entropy density
S = (a2*So+a3*Soo + (b1*dSdt + b2*dSdto + b3*dSdtoo)*dt)/a1;

% convert entropy desnity to temperature
T = (T0+273.15)*exp((S - X.*Dsx)./RHO./cP + Adbt./cP.*(Pt-Ptop));

% major component advection
advn_C = - advect(M.*cm,Um(2:end-1,:),Wm(:,2:end-1),h,{ADVN,''},[1,2],BCA) ...  % melt  advection
         - advect(X.*cx,Ux(2:end-1,:),Wx(:,2:end-1),h,{ADVN,''},[1,2],BCA);     % solid advection

% total rate of change
dCdt = advn_C;                                            

% residual of major component evolution
res_C = (a1*C-a2*Co-a3*Coo)/dt - (b1*dCdt + b2*dCdto + b3*dCdtoo);

% semi-implicit update of major component density
C = (a2*Co+a3*Coo + (b1*dCdt + b2*dCdto + b3*dCdtoo)*dt)/a1;

% apply minimum/maximum bounds
C = max(0, C );

% convert major component density to concentration
c = C./RHO;
c = c./sum(c,3);


%% *****  UPDATE PHASE EQUILIBRIUM  ***************************************

% update local phase equilibrium
var.c = reshape(c,Nz*Nx,cal.ncmp);  % in wt
var.T = reshape(T,Nz*Nx,1)-273.15;  % convert to C
var.P = reshape(Pt,Nz*Nx,1)/1e9;    % convert to GPa
var.f = 1-reshape(xq,Nz*Nx,1);      % in wt

[phs,cal]  =  meltmodel(var,cal,'E'); % cs and cl component prop in each phase VS T

mq = max(0,min(1,reshape(phs.f,Nz,Nx)));
xq = 1-mq;

cxq = reshape(phs.cs,Nz,Nx,cal.ncmp);
cmq = reshape(phs.cl,Nz,Nx,cal.ncmp);

cmq = cmq./sum(cmq,3);
cxq = cxq./sum(cxq,3);


%% *****  UPDATE PHASE PROPORTIONS  ***************************************

% phase advection rates
advn_X   = - advect(X,Ux(2:end-1,:),Wx(:,2:end-1),h,{ADVN,''},[1,2],BCA);
advn_M   = - advect(M,Um(2:end-1,:),Wm(:,2:end-1),h,{ADVN,''},[1,2],BCA);
advn_rho = advn_X+advn_M;

% phase mass transfer rates
Gx = (Gx + (xq.*RHO-X)./max(tau_r,4*dt))/2;
Gm = (Gm + (mq.*RHO-M)./max(tau_r,4*dt))/2;

% total rates of change
dXdt   = advn_X + Gx;
dMdt   = advn_M + Gm;

% residual of phase density evolution
res_X = (a1*X-a2*Xo-a3*Xoo)/dt - (b1*dXdt + b2*dXdto + b3*dXdtoo);
res_M = (a1*M-a2*Mo-a3*Moo)/dt - (b1*dMdt + b2*dMdto + b3*dMdtoo);

% semi-implicit update of phase fraction densities
X = (a2*Xo+a3*Xoo + (b1*dXdt + b2*dXdto + b3*dXdtoo)*dt)/a1;
M = (a2*Mo+a3*Moo + (b1*dMdt + b2*dMdto + b3*dMdtoo)*dt)/a1;

% apply minimum bound
X   = max(0, X );
M   = max(0, M );

% get dynamically evolving mixture density 
RHO = X+M;

%***  update phase fractions and component concentrations

% update phase fractions
x = X./RHO;
m = M./RHO;

% update phase entropies
sm = (S - X.*Dsx)./RHO;
sx = sm + Dsx;

% update major component phase composition
Kc  = cxq./cmq;
res = 1; tol = 1e-9;
while res>tol
    Kci = Kc;
    cm  = c./(m + x.*Kc); cm = cm./sum(cm,3);
    cx  = c./(m./Kc + x); cx = cx./sum(cx,3);
    Kc  = cx./cm;
    res = norm(Kci-Kc,'fro')./norm(Kc,'fro');
end
