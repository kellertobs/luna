%%*****  UPDATE PARAMETERS & AUXILIARY FIELDS  ****************************

% update phase oxide compositions
oxd  = permute(reshape(reshape(permute(c ,[2,3,1]),Nz*Nx,cal.nc)*cal.oxds,Nz,Nx,cal.nc),[3,1,2]);
oxdm = permute(reshape(reshape(permute(cm,[2,3,1]),Nz*Nx,cal.nc)*cal.oxds,Nz,Nx,cal.nc),[3,1,2]);
oxdx = permute(reshape(reshape(permute(cx,[2,3,1]),Nz*Nx,cal.nc)*cal.oxds,Nz,Nx,cal.nc),[3,1,2]);

% update mineral system numbers
foNo = squeeze(c(1,:,:)./(c(1,:,:) + c(2,:,:)));  % fo# = fo /(fo +fay)
pxNo = squeeze(c(3,:,:)./(c(3,:,:) + c(4,:,:)));  % px# = opx/(opx+cpx)
anNo = squeeze(c(5,:,:)./(c(5,:,:) + c(6,:,:)));  % an# = an /(an +ab )

foNom = squeeze(cm(1,:,:)./(cm(1,:,:) + cm(2,:,:)));
pxNom = squeeze(cm(3,:,:)./(cm(3,:,:) + cm(4,:,:)));
anNom = squeeze(cm(5,:,:)./(cm(5,:,:) + cm(6,:,:)));

foNox = squeeze(cx(1,:,:)./(cx(1,:,:) + cx(2,:,:)));
pxNox = squeeze(cx(3,:,:)./(cx(3,:,:) + cx(4,:,:)));
anNox = squeeze(cx(5,:,:)./(cx(5,:,:) + cx(6,:,:)));

% update phase densities
rhom = squeeze(1./sum(cm./rhom0.')) .* (1 - aT.*(T-T0-273.15));
rhox = squeeze(1./sum(cx./rhox0.')) .* (1 - aT.*(T-T0-273.15));

% convert weight to volume fraction, update bulk density
rho   = 1./(m./rhom + x./rhox);  
rho([1 end],:) = rho([2 end-1],:);  
rho(:,[1 end]) = rho(:,[2 end-1]);

chi   = x.*rho./rhox;
mu    = m.*rho./rhom;

% update thermal properties
Ds    = x.*Dsx;

% update effective viscosity
wtm      = zeros(Nz*Nx,12);
wtm(:,1) = reshape(oxd(1,:,:),Nz*Nx,1); % SiO2
wtm(:,3) = reshape(oxd(2,:,:),Nz*Nx,1); % Al2O3
wtm(:,4) = reshape(oxd(3,:,:),Nz*Nx,1); % FeO
wtm(:,6) = reshape(oxd(4,:,:),Nz*Nx,1); % MgO
wtm(:,7) = reshape(oxd(5,:,:),Nz*Nx,1); % CaO
wtm(:,8) = reshape(oxd(6,:,:),Nz*Nx,1); % Na2O
etam  = reshape(grdmodel08(wtm,T(:)-273.15),Nz,Nx);
etax  = etax0.* ones(size(x));                                             % constant crysta viscosity

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
eta    = squeeze(sum(ff.*kv.*thtv,1));
eta    = (1./etamax + 1./eta).^-1 + etamin;
eta([1 end],:) = eta([2 end-1],:);  
eta(:,[1 end]) = eta(:,[2 end-1]);
etaco  = (eta(1:end-1,1:end-1)+eta(2:end,1:end-1) ...                      % effective viscosity in cell corners
       +  eta(1:end-1,2:end  )+eta(2:end,2:end  ))./4;

% get segregation coefficients
Csgr = ((1-ff)./[dx;dx].^2.*kv.*thtv).^-1;

Csgr_x = squeeze(Csgr(1,:,:)) + 1e-18;  Csgr_x([1 end],:) = Csgr_x([2 end-1],:);  Csgr_x(:,[1 end]) = Csgr_x(:,[2 end-1]);
Csgr_m = squeeze(Csgr(2,:,:)) + 1e-18;  Csgr_m([1 end],:) = Csgr_m([2 end-1],:);  Csgr_m(:,[1 end]) = Csgr_m(:,[2 end-1]);

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
txx = eta .* exx;                                                        % x-normal stress
tzz = eta .* ezz;                                                        % z-normal stress
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
wx = ((rhox(1:end-1,:)+rhox(2:end,:))/2-rhoref)*g0.*(Csgr_x(1:end-1,:).*Csgr_x(2:end,:)).^0.5; % crystal segregation speed
for i = 1:round(delta)
    wx(2:end-1,2:end-1) = wx(2:end-1,2:end-1) + diff(wx(:,2:end-1),2,1)./8 + diff(wx(2:end-1,:),2,2)./8;
    wx(1  ,:)     = min(1,1-top).*wx(1  ,:);
    wx(end,:)     = min(1,1-bot).*wx(end,:);
    wx(:,[1 end]) = -sds*wx(:,[2 end-1]);
end

wm = ((rhom(1:end-1,:)+rhom(2:end,:))/2-rhoref)*g0.*(Csgr_m(1:end-1,:).*Csgr_m(2:end,:)).^0.5.*((chi(1:end-1,:)+chi(2:end,:))./2).^2; % melt segregation speed
for i = 1:round(delta)
    wm(2:end-1,2:end-1) = wm(2:end-1,2:end-1) + diff(wm(:,2:end-1),2,1)./8 + diff(wm(2:end-1,:),2,2)./8;
    wm(1  ,:)     = min(1,1-top).*wm(1  ,:);
    wm(end,:)     = min(1,1-bot).*wm(end,:);
    wm(:,[1 end]) = -sds*wm(:,[2 end-1]);
end

% update volume source
Div_rhoV =  + advection(rho.*x,0.*U,wx,h,ADVN,'flx') ...
            + advection(rho.*m,0.*U,wm,h,ADVN,'flx') ...
            + advection(rho   ,   U,W ,h,ADVN,'flx');
VolSrc  = -((rho-rhoo)./dt + Div_rhoV - rho.*Div_V)./rho;

UBG    = - 2*mean(mean(VolSrc(2:end-1,2:end-1)))./2 .* (L/2-XXu);
WBG    = - 0*mean(mean(VolSrc(2:end-1,2:end-1)))./2 .* (D/2-ZZw);
