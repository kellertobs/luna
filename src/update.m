%%*****  UPDATE PARAMETERS & AUXILIARY FIELDS  ****************************

% update phase oxide compositions
c_oxd  = reshape(reshape(c ,Nz*Nx,cal.nc)*cal.oxd,Nz,Nx,cal.nc);
cm_oxd = reshape(reshape(cm,Nz*Nx,cal.nc)*cal.oxd,Nz,Nx,cal.nc);
cx_oxd = reshape(reshape(cx,Nz*Nx,cal.nc)*cal.oxd,Nz,Nx,cal.nc);

% update phase densities
rhom = squeeze(sum(permute(cm,[3,1,2])./rhom0.')).^-1 .* (1 - aT.*(T-T0-273.15) + bPm.*(Pt-Ptop));
rhox = squeeze(sum(permute(cx,[3,1,2])./rhox0.')).^-1 .* (1 - aT.*(T-T0-273.15) + bPx.*(Pt-Ptop));

% convert weight to volume fraction, update bulk density
rho   = 1./(m./rhom + x./rhox);  

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
etam  = reshape(grdmodel08(wtm,T(:)-273.15),Nz,Nx);
etax  = etax0.* ones(size(x));                                             % constant crystal viscosity

% get permission weights
kv = permute(cat(3,etax,etam),[3,1,2]);
Mv = permute(repmat(kv,1,1,1,2),[4,1,2,3])./permute(repmat(kv,1,1,1,2),[1,4,2,3]);
 
ff = max(1e-4,min(1-1e-4,permute(cat(3,chi,mu),[3,1,2])));
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
Csgr = ((1-ff)./dx.^2.*(sgrreg.*kv.*thtv)).^-1 + 1e-18;

Csgr_x = squeeze(Csgr(1,:,:));  Csgr_x([1 end],:) = Csgr_x([2 end-1],:);  Csgr_x(:,[1 end]) = Csgr_x(:,[2 end-1]);
Csgr_m = squeeze(Csgr(2,:,:));  Csgr_m([1 end],:) = Csgr_m([2 end-1],:);  Csgr_m(:,[1 end]) = Csgr_m(:,[2 end-1]);
Csgr_m = Csgr_m.*max(TINY,chi).^2; % dampen melt segregation at high melt fraction

% update velocity divergence
Div_V(2:end-1,2:end-1) = ddz(W(:,2:end-1),h) ...                           % get velocity divergence
                       + ddx(U(2:end-1,:),h);
Div_V([1 end],:) = Div_V([2 end-1],:);                                     % apply boundary conditions
Div_V(:,[1 end]) = Div_V(:,[2 end-1]);

% update strain rates
exx(:,2:end-1) = diff(U,1,2)./h - Div_V(:,2:end-1)./3;                     % x-normal strain rate
exx([1 end],:) = exx([2 end-1],:);                                         % apply boundary conditions
exx(:,[1 end]) = exx(:,[2 end-1]);
ezz(2:end-1,:) = diff(W,1,1)./h - Div_V(2:end-1,:)./3;                     % z-normal strain rate
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
% wx = ((rhox(1:end-1,:)+rhox(2:end,:))/2-(rho(1:end-1,:)+rho(2:end,:))/2).*g0.*2./(1./Csgr_x(1:end-1,:)+1./Csgr_x(2:end,:)); % melt segregation speed
wx = ((rhox(1:end-1,:)+rhox(2:end,:))/2-(rho(1:end-1,:)+rho(2:end,:))/2).*g0.*min(Csgr_x(1:end-1,:),Csgr_x(2:end,:)); % melt segregation speed
wx(1  ,:)     = 0;
wx(end,:)     = 0;
wx(:,[1 end]) = wx(:,[2 end-1]);

wm = 0.*((rhom(1:end-1,:)+rhom(2:end,:))/2-(rho(1:end-1,:)+rho(2:end,:))/2).*g0.*2./(1./Csgr_m(1:end-1,:)+1./Csgr_m(2:end,:)); % melt segregation speed
wm(1  ,:)     = 0;
wm(end,:)     = 0;
wm(:,[1 end]) = wm(:,[2 end-1]);

% diffusion parameters
kT = kT0 * dffreg;
ks = kT./T;
kc = 0.*rho.*abs((rhox-rho).*g0.*Csgr_x.*dx);           % chemical diffusion by fluctuation in crystal segregation speed
kx = abs((rhox-rho).*g0.*Csgr_x.*dx) * dffreg;
kx = kx.*min(kT./cP./rho./10,[],'all')./max(kx,[],'all');

% heat dissipation (entropy production) rate
[grdTx,grdTz] = gradient(T,h);
diss =  exx(2:end-1,2:end-1).*txx(2:end-1,2:end-1) ...
     +  ezz(2:end-1,2:end-1).*tzz(2:end-1,2:end-1) ...
     +  2.*(exz(1:end-1,1:end-1)+exz(2:end,1:end-1)+exz(1:end-1,2:end)+exz(2:end,2:end))./4 ...
         .*(txz(1:end-1,1:end-1)+txz(2:end,1:end-1)+txz(1:end-1,2:end)+txz(2:end,2:end))./4 ...
     +  chi(2:end-1,2:end-1)./Csgr_x(2:end-1,2:end-1) .* ((wx(1:end-1,2:end-1)+wx(2:end,2:end-1))./2).^2 ...
     +  mu (2:end-1,2:end-1)./Csgr_m(2:end-1,2:end-1) .* ((wm(1:end-1,2:end-1)+wm(2:end,2:end-1))./2).^2 ...
     +  ks(2:end-1,2:end-1).*(grdTz(2:end-1,2:end-1).^2 + grdTx(2:end-1,2:end-1).^2);
                        
% update volume source
Div_rhoV =  + advect(rho(inz,inx).*m(inz,inx),Um(inz,:)-U(inz,:),Wm(:,inx)-W(:,inx),h,{ADVN,''   },[1,2],BCA) ...
            + advect(rho(inz,inx).*x(inz,inx),Ux(inz,:)-U(inz,:),Wx(:,inx)-W(:,inx),h,{ADVN,''   },[1,2],BCA) ...
            + advect(rho(inz,inx)            ,          U(inz,:),          W(:,inx),h,{ADVN,'vdf'},[1,2],BCA);
if step>0; VolSrc(inz,inx) = -((rho(inz,inx)-rhoo(inz,inx))./dt + Div_rhoV)./rho(inz,inx); end
% if step>0; VolSrc(inz,inx) = -((rho(inz,inx)-rhoo(inz,inx))./dt + theta.*Div_rhoV + (1-theta).*Div_rhoVo)./rho(inz,inx); end

UBG    = - 1*mean(mean(VolSrc(inz,inx)))./2 .* (L/2-XXu);
WBG    = - 1*mean(mean(VolSrc(inz,inx)))./2 .* (D/2-ZZw);
