%%*****  UPDATE PARAMETERS & AUXILIARY FIELDS  ****************************

% update phase oxide compositions
c_oxd  = reshape(reshape(c ,Nz*Nx,cal.nc)*cal.oxd,Nz,Nx,cal.nc);
cm_oxd = reshape(reshape(cm,Nz*Nx,cal.nc)*cal.oxd,Nz,Nx,cal.nc);
cx_oxd = reshape(reshape(cx,Nz*Nx,cal.nc)*cal.oxd,Nz,Nx,cal.nc);

% update phase mineral end-member compositions
c_mnr(:,:,cal.mnr_for) = squeeze(c(:,:,cal.for)).*100;
c_mnr(:,:,cal.mnr_fay) = squeeze(c(:,:,cal.fay)).*100;
c_mnr(:,:,cal.mnr_px1) = squeeze(c(:,:,cal.opx)).*100;
c_mnr(:,:,cal.mnr_px2) = squeeze(c(:,:,cal.cpx)).*100;
c_mnr(:,:,cal.mnr_px3) = squeeze(c(:,:,cal.eut).*cal.eut_mnr(cal.mnr_px3)).*100;
c_mnr(:,:,cal.mnr_ant) = squeeze(c(:,:,cal.ant)).*100;
c_mnr(:,:,cal.mnr_alb) = squeeze(c(:,:,cal.eut).*cal.eut_mnr(cal.mnr_alb)).*100;
c_mnr(:,:,cal.mnr_spn) = squeeze(c(:,:,cal.eut).*cal.eut_mnr(cal.mnr_spn)).*100;
c_mnr(:,:,cal.mnr_qtz) = squeeze(c(:,:,cal.eut).*cal.eut_mnr(cal.mnr_qtz)).*100;

cx_mnr(:,:,cal.mnr_for) = squeeze(cx(:,:,cal.for)).*100;
cx_mnr(:,:,cal.mnr_fay) = squeeze(cx(:,:,cal.fay)).*100;
cx_mnr(:,:,cal.mnr_px1) = squeeze(cx(:,:,cal.opx)).*100;
cx_mnr(:,:,cal.mnr_px2) = squeeze(cx(:,:,cal.cpx)).*100;
cx_mnr(:,:,cal.mnr_px3) = squeeze(cx(:,:,cal.eut).*cal.eut_mnr(cal.mnr_px3)).*100;
cx_mnr(:,:,cal.mnr_ant) = squeeze(cx(:,:,cal.ant)).*100;
cx_mnr(:,:,cal.mnr_alb) = squeeze(cx(:,:,cal.eut).*cal.eut_mnr(cal.mnr_alb)).*100;
cx_mnr(:,:,cal.mnr_spn) = squeeze(cx(:,:,cal.eut).*cal.eut_mnr(cal.mnr_spn)).*100;
cx_mnr(:,:,cal.mnr_qtz) = squeeze(cx(:,:,cal.eut).*cal.eut_mnr(cal.mnr_qtz)).*100;

cm_mnr(:,:,cal.mnr_for) = squeeze(cm(:,:,cal.for)).*100;
cm_mnr(:,:,cal.mnr_fay) = squeeze(cm(:,:,cal.fay)).*100;
cm_mnr(:,:,cal.mnr_px1) = squeeze(cm(:,:,cal.opx)).*100;
cm_mnr(:,:,cal.mnr_px2) = squeeze(cm(:,:,cal.cpx)).*100;
cm_mnr(:,:,cal.mnr_px3) = squeeze(cm(:,:,cal.eut).*cal.eut_mnr(cal.mnr_px3)).*100;
cm_mnr(:,:,cal.mnr_ant) = squeeze(cm(:,:,cal.ant)).*100;
cm_mnr(:,:,cal.mnr_alb) = squeeze(cm(:,:,cal.eut).*cal.eut_mnr(cal.mnr_alb)).*100;
cm_mnr(:,:,cal.mnr_spn) = squeeze(cm(:,:,cal.eut).*cal.eut_mnr(cal.mnr_spn)).*100;
cm_mnr(:,:,cal.mnr_qtz) = squeeze(cm(:,:,cal.eut).*cal.eut_mnr(cal.mnr_qtz)).*100;

% update phase densities
rhom = squeeze(sum(permute(cm,[3,1,2])./rhom0.')).^-1 .* (1 - aT.*(T-T0-273.15) + bPm.*(Pt-Ptop));
rhox = squeeze(sum(permute(cx,[3,1,2])./rhox0.')).^-1 .* (1 - aT.*(T-T0-273.15) + bPx.*(Pt-Ptop));

% convert weight to volume fraction, update bulk density
rho   = 1./(m./rhom + x./rhox);

rhofz = (rho(1:end-1,:)+rho(2:end,:))/2;

chi   = x.*rho./rhox;
mu    = m.*rho./rhom;

% determine adaptive max viscosity / min segregation coefficient
if exist('eta','var'); etamax = 1e+6.*min(eta (:)); else; etamax = 1e18; end

% update effective viscosity
wtm      = zeros(Nz*Nx,12);
wtm(:,1) = reshape(cm_oxd(:,:,1),Nz*Nx,1); % SiO2
wtm(:,3) = reshape(cm_oxd(:,:,2),Nz*Nx,1); % Al2O3
wtm(:,4) = reshape(cm_oxd(:,:,3),Nz*Nx,1); % FeO
wtm(:,6) = reshape(cm_oxd(:,:,4),Nz*Nx,1); % MgO
wtm(:,7) = reshape(cm_oxd(:,:,5),Nz*Nx,1); % CaO
wtm(:,8) = reshape(cm_oxd(:,:,6),Nz*Nx,1); % Na2O
etam  = max(1e-2,reshape(grdmodel08(wtm,T(:)-273.15),Nz,Nx));
etax  = etax0.* ones(size(x));                                             % constant crystal viscosity

% get permission weights
kv = permute(cat(3,etax,etam),[3,1,2]);
Mv = permute(repmat(kv,1,1,1,2),[4,1,2,3])./permute(repmat(kv,1,1,1,2),[1,4,2,3]);
 
ff = max(TINY,min(1-TINY,permute(cat(3,chi,mu),[3,1,2])));
FF = permute(repmat(ff,1,1,1,2),[4,1,2,3]);
Sf = (FF./BB).^(1./CC);  Sf = Sf./sum(Sf,2);
Xf = sum(AA.*Sf,2).*FF + (1-sum(AA.*Sf,2)).*Sf;

% get momentum and volume permissions
thtv = squeeze(prod(Mv.^Xf,2));

% get effective viscosity
eta  = squeeze(sum(ff.*kv.*thtv,1));
eta  = (1./etamax + 1./(etareg*eta)).^-1;% + etabnd(1);
eta([1 end],:) = eta([2 end-1],:);  
eta(:,[1 end]) = eta(:,[2 end-1]);
etaco  = (eta(1:end-1,1:end-1)+eta(2:end,1:end-1) ...                      % effective viscosity in cell corners
       +  eta(1:end-1,2:end  )+eta(2:end,2:end  ))./4;

% get segregation coefficients
dd   = permute(cat(3,d0*ones(size(mu)),d0*(1-mu)),[3,1,2]);
Ksgr = ((1-ff)./dd.^2.*(sgrreg.*kv.*thtv)).^-1 + TINY^2;

Ksgr_x = squeeze(Ksgr(1,:,:));  Ksgr_x([1 end],:) = Ksgr_x([2 end-1],:);  Ksgr_x(:,[1 end]) = Ksgr_x(:,[2 end-1]);
Ksgr_m = squeeze(Ksgr(2,:,:));  Ksgr_m([1 end],:) = Ksgr_m([2 end-1],:);  Ksgr_m(:,[1 end]) = Ksgr_m(:,[2 end-1]);

% update velocity divergence
Div_V(2:end-1,2:end-1) = ddz(W(:,2:end-1),h) ...                           % get velocity divergence
                       + ddx(U(2:end-1,:),h);
Div_V([1 end],:) = Div_V([2 end-1],:);                                     % apply boundary conditions
Div_V(:,[1 end]) = Div_V(:,[2 end-1]);

% update strain rates
exx(:,2:end-1) = diff(U,1,2)./h - Div_V(:,2:end-1)./2;                     % x-normal strain rate
exx([1 end],:) = exx([2 end-1],:);                                         % apply boundary conditions
exx(:,[1 end]) = exx(:,[2 end-1]);
ezz(2:end-1,:) = diff(W,1,1)./h - Div_V(2:end-1,:)./2;                     % z-normal strain rate
ezz([1 end],:) = ezz([2 end-1],:);                                         % apply boundary conditions
ezz(:,[1 end]) = ezz(:,[2 end-1]);
exz            = 1/2.*(diff(U,1,1)./h+diff(W,1,2)./h);                     % shear strain rate

% update stresses
txx = eta .* exx;                                                          % x-normal stress
tzz = eta .* ezz;                                                          % z-normal stress
txz = etaco .* exz;                                                        % xz-shear stress

% update tensor magnitudes
eII(2:end-1,2:end-1) = 1e-16 + (0.5.*(exx(2:end-1,2:end-1).^2 + ezz(2:end-1,2:end-1).^2 ...  % get strain rate magnitude
                             + 2.*(exz(1:end-1,1:end-1).^2.*exz(2:end,1:end-1).^2.*exz(1:end-1,2:end).^2.*exz(2:end,2:end).^2).^0.25)).^0.5;
eII(:,[1 end]) = eII(:,[2 end-1]);
eII([1 end],:) = eII([2 end-1],:);

tII(2:end-1,2:end-1) = 1e-16 + (0.5.*(txx(2:end-1,2:end-1).^2 + tzz(2:end-1,2:end-1).^2 ...  % get stress magnitude
                             + 2.*(txz(1:end-1,1:end-1).^2.*txz(2:end,1:end-1).^2.*txz(1:end-1,2:end).^2.*txz(2:end,2:end).^2).^0.25)).^0.5;  
tII(:,[1 end]) = tII(:,[2 end-1]);
tII([1 end],:) = tII([2 end-1],:);

% update phase segregation speeds
wm = ((rhom(1:end-1,:)+rhom(2:end,:))/2-(rho(1:end-1,:)+rho(2:end,:))/2).*g0.*2./(1./Ksgr_m(1:end-1,:)+1./Ksgr_m(2:end,:)); % melt segregation speed
wm([1 end],:) = 0;
wm(:,[1 end]) = wm(:,[2 end-1]);

wx = ((rhox(1:end-1,:)+rhox(2:end,:))/2-(rho(1:end-1,:)+rho(2:end,:))/2).*g0.*2./(1./Ksgr_x(1:end-1,:)+1./Ksgr_x(2:end,:)); % solid segregation speed
wx([1 end],:) = 0;
wx(:,[1 end]) = wx(:,[2 end-1]);

% diffusion parameters
if Nx==3
    [~,grdrhoz] = gradient(rho,h);
    Vel = abs(grdrhoz*h*g0*h^2./eta);                                      % estimate convective speed in 1D case
end
kW  = Vel/10*h/10;                                                         % convection fluctuation diffusivity
kwx = abs((rhox-rho).*g0.*Ksgr_x*d0*10);                                   % segregation fluctuation diffusivity
kx  = chi.*(kwx + kW);                                                     % solid fraction diffusion 
kT  = kT0 + (chi.*kwx + kW).*rho.*cP;                                      % heat diffusion
ks  = kT./T;                                                               % entropy diffusion

% heat dissipation (entropy production) rate
if Nz==3 && Nx==3  
    diss = 0.*T(inz,inx);  % no dissipation in 0-D mode (no diffusion, no shear deformation, no segregation)
else
    [grdTx,grdTz] = gradient(T,h);
    diss = ks(inz,inx).*(grdTz(inz,inx).^2 + grdTx(inz,inx).^2) ...
        + exx(inz,inx).*txx(inz,inx) + ezz(inz,inx).*tzz(inz,inx) ...
        + 2.*(exz(1:end-1,1:end-1)+exz(2:end,1:end-1)+exz(1:end-1,2:end)+exz(2:end,2:end))./4 ...
           .*(txz(1:end-1,1:end-1)+txz(2:end,1:end-1)+txz(1:end-1,2:end)+txz(2:end,2:end))./4 ...
        +  mu(inz,inx)./Ksgr_m(inz,inx) .* ((wm(inz,inx)+wm(inz,inx))./2).^2 ...
        + chi(inz,inx)./Ksgr_x(inz,inx) .* ((wx(inz,inx)+wx(inz,inx))./2).^2;
end
                        
% update volume source
if step>0 && ~restart
    Div_rhoV = + advect(M(inz,inx),Um(inz,:),Wm(:,inx),h,{ADVN,''},[1,2],BCA) ...  % melt advection
               + advect(X(inz,inx),Ux(inz,:),Wx(:,inx),h,{ADVN,''},[1,2],BCA);     % xtal advection
    F_DivV   = (alpha1*rho(inz,inx) - alpha2*rhoo(inz,inx) - alpha3*rhooo(inz,inx))./dt + (beta1*Div_rhoV + beta2*Div_rhoVo + beta3*Div_rhoVoo);  % get residual of mixture mass conservation
    VolSrc   = Div_V(inz,inx) - F_DivV./rho(inz,inx);  % correct volume source term by scaled residual
end

UBG    = - 0*mean(mean(VolSrc))./2 .* (L/2-XXu);
WBG    = - 2*mean(mean(VolSrc))./2 .* (D/2-ZZw);
