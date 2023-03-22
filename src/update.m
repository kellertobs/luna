%%*****  UPDATE PARAMETERS & AUXILIARY FIELDS  ****************************

% update phase oxide compositions
c_oxd  = reshape(reshape(c ,Nz*Nx,cal.ncmp)*cal.cmp_oxd,Nz,Nx,cal.noxd);
cm_oxd = reshape(reshape(cm,Nz*Nx,cal.ncmp)*cal.cmp_oxd,Nz,Nx,cal.noxd);
cx_oxd = reshape(reshape(cx,Nz*Nx,cal.ncmp)*cal.cmp_oxd,Nz,Nx,cal.noxd);

% update phase mineral end-member compositions
c_mnr  = reshape(reshape(c ,Nz*Nx,cal.ncmp)*cal.cmp_mnr,Nz,Nx,cal.nmnr);
cm_mnr = reshape(reshape(cm,Nz*Nx,cal.ncmp)*cal.cmp_mnr,Nz,Nx,cal.nmnr);
cx_mnr = reshape(reshape(cx,Nz*Nx,cal.ncmp)*cal.cmp_mnr,Nz,Nx,cal.nmnr);


% update phase densities
wtm      = zeros(Nz*Nx,9);
wtm(:,1) = reshape(cm_oxd(:,:,1),Nz*Nx,1); % SiO2
wtm(:,2) = reshape(cm_oxd(:,:,2),Nz*Nx,1); % SiO2
wtm(:,3) = reshape(cm_oxd(:,:,3),Nz*Nx,1); % Al2O3
wtm(:,4) = reshape(cm_oxd(:,:,4),Nz*Nx,1); % FeO
wtm(:,5) = reshape(cm_oxd(:,:,5),Nz*Nx,1); % MgO
wtm(:,6) = reshape(cm_oxd(:,:,6),Nz*Nx,1); % CaO
wtm(:,7) = reshape(cm_oxd(:,:,7),Nz*Nx,1); % Na2O
rhom     = reshape(DensityX(wtm,T0.*ones(size(T(:))),Ptop/1e8.*ones(size(Pt(:)))),Nz,Nx) .* (1 - aT.*(T-T0-273.15) + bPm.*(Pt-Ptop));
rhox     = reshape(sum(reshape(cx_mnr,Nz*Nx,cal.nmnr)./rhox0,2).^-1              ,Nz,Nx) .* (1 - aT.*(T-T0-273.15) + bPx.*(Pt-Ptop));

% convert weight to volume fraction, update bulk density
rho   = 1./(m./rhom + x./rhox);

chi   = x.*rho./rhox;
mu    = m.*rho./rhom;

if Nz==1; Pt = Ptop.*ones(size(Tp)); else
    rhofz       = (rho(1:end-1,:)+rho(2:end,:))/2;
    Pt(1,:)     = repmat(mean(rhofz(1,:),2).*g0.*h/2,1,Nx) + Ptop;
    Pt(2:end,:) = Pt(1,:) + repmat(cumsum(mean(rhofz,2).*g0.*h),1,Nx);
end

% determine adaptive max viscosity / min segregation coefficient
if exist('eta','var'); etamax = etacntr.*min(eta (:)); else; etamax = 1e18; end

% update effective viscosity
etam  = max(1e-2,reshape(Giordano08(wtm,T(:)-273.15),Nz,Nx));
etax  = etax0.* ones(size(x));

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
eta    = squeeze(sum(ff.*kv.*thtv,1)); if Nx==1; eta = eta.'; end
eta    = (etamax.^-0.5 + (etareg*eta).^-0.5).^-2;
etaco  = (eta([1,1:end],[1  ,1:end]).*eta([1:end,end],[1  ,1:end]) ...
       .* eta([1,1:end],[1:end,end]).*eta([1:end,end],[1:end,end])).^0.25;

% get segregation coefficients
dd   = permute(cat(3,d0*ones(size(mu)),d0*(1-mu)),[3,1,2]);
Ksgr = ((1-ff)./dd.^2.*(sgrreg.*kv.*thtv)).^-1 + TINY;

Ksgr_x = squeeze(Ksgr(1,:,:)); if Nx==1; Ksgr_x = Ksgr_x.'; end
Ksgr_m = squeeze(Ksgr(2,:,:)); if Nx==1; Ksgr_m = Ksgr_m.'; end

% diffusion parameters
if Nx==1 && Nz>1
    drhodz = ddz(rhofz([1,1:end,end],:),h);
    Vel = abs(drhodz*h*g0*h^2./eta);                                       % estimate convective speed in 1D case
end
kW  = Vel/10*h/10;                                                         % convection fluctuation diffusivity
kwx = abs((rhox-rho).*g0.*Ksgr_x*d0*10);                                   % segregation fluctuation diffusivity
kx  = chi.*mu.*(kwx + kW + mink);                                          % solid fraction diffusion 
kT  = kT0 + mu.*rho.*cP.*(chi.*kwx + kW + mink);                % heat diffusion
ks  = kT./T;                                                               % entropy diffusion

% update velocity divergence
Div_V = ddz(W(:,2:end-1),h) + ddx(U(2:end-1,:),h);                         % get velocity divergence
      
% update strain rates
exx = diff(U(2:end-1,:),1,2)./h - Div_V./2;                                % x-normal strain rate
ezz = diff(W(:,2:end-1),1,1)./h - Div_V./2;                                % z-normal strain rate
exz = 1/2.*(diff(U,1,1)./h+diff(W,1,2)./h);                                % shear strain rate

% update stresses
txx = eta   .* exx;                                                        % x-normal stress
tzz = eta   .* ezz;                                                        % z-normal stress
txz = etaco .* exz;                                                        % xz-shear stress

% update tensor magnitudes
eII = (0.5.*(exx.^2 + ezz.^2 ...
       + 2.*(exz(1:end-1,1:end-1).^2+exz(2:end,1:end-1).^2 ...
       +     exz(1:end-1,2:end  ).^2+exz(2:end,2:end  ).^2)/4)).^0.5 + TINY;

tII = (0.5.*(txx.^2 + tzz.^2 ...
       + 2.*(txz(1:end-1,1:end-1).^2+txz(2:end,1:end-1).^2 ...
       +     txz(1:end-1,2:end  ).^2+txz(2:end,2:end).^2)/4)).^0.5 + TINY;

% heat dissipation (entropy production) rate
if Nz==1 && Nx==1  
    diss = 0.*T;  % no dissipation in 0-D mode (no diffusion, no shear deformation, no segregation)
else
    [grdTx,grdTz] = gradient(T([1,1:end,end],[1,1:end,end]),h);
    diss = ks.*(grdTz(2:end-1,2:end-1).^2 + grdTx(2:end-1,2:end-1).^2) ...
         + exx.*txx + ezz.*tzz ...
         + 2.*(exz(1:end-1,1:end-1)+exz(2:end,1:end-1)+exz(1:end-1,2:end)+exz(2:end,2:end))./4 ...
            .*(txz(1:end-1,1:end-1)+txz(2:end,1:end-1)+txz(1:end-1,2:end)+txz(2:end,2:end))./4 ...
         +  mu./Ksgr_m .* ((wm(1:end-1,2:end-1)+wm(2:end,2:end-1))./2).^2 ...
         + chi./Ksgr_x .* ((wx(1:end-1,2:end-1)+wx(2:end,2:end-1))./2).^2;
end

