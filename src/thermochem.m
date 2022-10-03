%% *****  THERMO-CHEMICAL EVOLUTION  **************************************

% store previous iteration
Ti = T; xi = x;

% update temperature
if entr_mth  % use entropy equation to evolve heat
    
    advn_S = - advection(rho.*m.*sm,Um,Wm,h,ADVN,'flx') ...                % heat advection
             - advection(rho.*x.*sx,Ux,Wx,h,ADVN,'flx');
    
    qTz    = - (ks(1:end-1,:)+ks(2:end,:))/2 .* ddz(T,h);                  % heat diffusion z-flux
    qTx    = - (ks(:,1:end-1)+ks(:,2:end))/2 .* ddx(T,h);                  % heat diffusion x-flux
    diff_T(2:end-1,2:end-1) = (- ddz(qTz(:,2:end-1),h)  ...                % heat diffusion
                               - ddx(qTx(2:end-1,:),h));
    
    diss_h(2:end-1,2:end-1) = diss ./ T(2:end-1,2:end-1);
                        
    bndT = zeros(size(T));
    if ~isnan(Ttop); bndT = bndT + ((Ttop+273.15)-T)./tau_T .* topshape; end % impose top boundary layer
    if ~isnan(Tbot); bndT = bndT + ((Tbot+273.15)-T)./tau_T .* botshape; end % impose bot boundary layer
    bndS = rho.*cP.*bndT./T;
    
    dSdt = advn_S + diff_T + diss_h + bndS;                                % total rate of change
    
    S = So + (theta.*dSdt + (1-theta).*dSdto).*dt;                         % explicit update of entropy
    S([1 end],:) = S([2 end-1],:);                                         % apply zero flux boundary conditions
    S(:,[1 end]) = S(:,[2 end-1]);

    s = S./rho;
    
    T = (T0+273.15) * exp(S./rho./cP - x.*Dsx./cP + aT./rhoref./cP.*(Pt-Ptop));   % convert entropy to temperature
    
else  % use temperature equation to evolve heat

    advn_T = - advection(T,Ubar,Wbar,h,ADVN,'adv');                        % heat advection
    
    adbt_h =   T.*aT./rhoref./cP .* advection(Pt,Ubar,Wbar,h,ADVN,'adv');  % adiabatic heat
    
    latn_h = - T.*Dsx./rho./cP .* Gx;                                      % latent heat
                           
    qTz    = - kT .* ddz(T,h);                                             % heat diffusion z-flux
    qTx    = - kT .* ddx(T,h);                                             % heat diffusion x-flux
    diff_T(2:end-1,2:end-1) = (- ddz(qTz(:,2:end-1),h)  ...                % heat diffusion
                               - ddx(qTx(2:end-1,:),h)) ...
                               ./ rho(2:end-1,2:end-1)./cP;
    
    diss_h(2:end-1,2:end-1) = diss ./ rho(2:end-1,2:end-1)./cP;
                        
    bndT = zeros(size(T));
    if ~isnan(Ttop); bndT = bndT + ((Ttop+273.15)-T)./tau_T .* topshape; end % impose top boundary layer
    if ~isnan(Tbot); bndT = bndT + ((Tbot+273.15)-T)./tau_T .* botshape; end % impose bot boundary layer
    bndS = rho.*cP.*bndT./T;
    
    dTdt = advn_T + diff_T + adbt_h + latn_h + diss_h + bndT;              % total rate of change
    
    T = To + (theta.*dTdt + (1-theta).*dTdto).*dt;                         % explicit update of temperature
    T([1 end],:) = T([2 end-1],:);                                         % apply zero-flux boundary conditions
    T(:,[1 end]) = T(:,[2 end-1]);

    S = rho.*(cP.*log(T/(T0+273.15)) + x.*Dsx - aT./rhoref.*Pt);                    % convert temperature to entropy
    diss_h(2:end-1,2:end-1) = diss ./ T(2:end-1,2:end-1);                  % convert dissipation rate to entropy units

end

% update major component 
for i = 2:cal.nc
    advn_C = advection(rho.*m.*squeeze(cm(i,:,:)),Um,Wm,h,ADVN,'flx') ...
           + advection(rho.*x.*squeeze(cx(i,:,:)),Ux,Wx,h,ADVN,'flx');
    
    qcz   = - (kc(1:end-1,:)+kc(2:end,:))/2 .* ddz(squeeze(c(i,:,:)),h);   % major component diffusion z-flux
    qcx   = - (kc(:,1:end-1)+kc(:,2:end))/2 .* ddx(squeeze(c(i,:,:)),h);   % major component diffusion x-flux
    diff_c(2:end-1,2:end-1) = - ddz(qcz(:,2:end-1),h) ...                  % major component diffusion
                              - ddx(qcx(2:end-1,:),h);
    
    dCdt(i,:,:) = - advn_C + diff_c;                                       % total rate of change
    
    C(i,:,:) = Co(i,:,:) + (theta.*dCdt(i,:,:) + (1-theta).*dCdto(i,:,:)).*dt;  % explicit update of major component density
    C(i,:,:) = max(TINY,min(rho-TINY,squeeze(C(i,:,:))));
    C(i,[1 end],:) = C(i,[2 end-1],:);                                     % apply boundary conditions
    C(i,:,[1 end]) = C(i,:,[2 end-1]);
end

% convert component densities to concentrations
C(1,:,:) = max(TINY,min(rho-TINY, rho - squeeze(sum(C(2:end,:,:))) ));
for i = 1:cal.nc; c(i,:,:) = squeeze(C(i,:,:))./rho; end

% sumC = squeeze(sum(C));
% for i = 1:cal.nc; c(i,:,:) = squeeze(C(i,:,:))./sumC; end
% for i = 1:cal.nc; C(i,:,:) = squeeze(c(i,:,:)).*rho; end


%% *****  UPDATE PHASE PROPORTIONS  ***************************************

% update local phase equilibrium
for i = 1:cal.nc
    var.c(:,i) = reshape(c(i,:,:),Nx*Nz,1);  % in wt
end
var.T = T(:)-273.15;  % convert to C
var.P = Pt(:)/1e9;    % convert to GPa
var.f = 1-xq(:);      % in wt

[phs,cal]  =  meltmodel(var,cal,'E'); % cs and cl component prop in each phase VS T

mq = max(TINY,min(1-TINY,reshape(phs.f ,Nz,Nx)));
xq = 1-mq;
for i = 2:cal.nc
    cxq(i,:,:) = reshape(phs.cs(:,i),Nz,Nx);
    cmq(i,:,:) = reshape(phs.cl(:,i),Nz,Nx);
end
% for i = 1:cal.nc; cmq(i,:,:) = squeeze(cmq(i,:,:))./squeeze(sum(cmq)); end
% for i = 1:cal.nc; cxq(i,:,:) = squeeze(cxq(i,:,:))./squeeze(sum(cxq)); end
cmq(1,:,:) = 1 - squeeze(sum(cmq(2:end,:,:)));
cxq(1,:,:) = 1 - squeeze(sum(cxq(2:end,:,:)));

% update crystal fraction
if diseq  % quasi-equilibrium approach
    
    Gx = lambda.*Gx + (1-lambda) .* (xq-x).*rho./(4*dt);
        
    advn_X = - advection(rho.*x,Ux,Wx,h,ADVN,'flx');

    qxz   = - (kx(1:end-1,:)+kx(2:end,:))/2 .* ddz(x,h);                   % crystal fraction diffusion z-flux
    qxx   = - (kx(:,1:end-1)+kx(:,2:end))/2 .* ddx(x,h);                   % crystal fraction diffusion x-flux
    diff_x(2:end-1,2:end-1) = - ddz(qxz(:,2:end-1),h) ...                  % crystal fraction diffusion
                              - ddx(qxx(2:end-1,:),h);

    dXdt   = advn_X + diff_x + Gx;                                         % total rate of change
    
    X = Xo + (theta.*dXdt + (1-theta).*dXdto).*dt;                         % explicit update of crystal fraction
    X = min(rho.*(1-TINY),max(rho.*TINY,X));                               % enforce [0,1] limit
    X([1 end],:) = X([2 end-1],:);                                         % apply boundary conditions
    X(:,[1 end]) = X(:,[2 end-1]);
    
    x  = X./rho;
    
else  % equilibrium approach
    
    x  =  lambda.*x + (1-lambda).*xq;
    Gx = (rho.*x-rhoo.*xo)./dt + advection(rho.*x,Ux,Wx,h,ADVN,'flx');     % reconstruct crystallisation rate
    
end

% update melt fraction
m = min(1-TINY,max(TINY,1-x));

% update phase entropies & compositions
if step>0
    
    % phase entropies
    sm = s - x.*Dsx;                                                       
    sx = sm + Dsx;
    
    % phase compositions
    Kc = cxq./cmq;
    for i = 2:cal.nc
        cm(i,:,:) = squeeze(c(i,:,:))./(m + x.*squeeze(Kc(i,:,:)));
        cx(i,:,:) = squeeze(c(i,:,:))./(m./squeeze(Kc(i,:,:)) + x);
    end
%     for i = 1:cal.nc; cm(i,:,:) = squeeze(cm(i,:,:))./squeeze(sum(cm)); end
%     for i = 1:cal.nc; cx(i,:,:) = squeeze(cx(i,:,:))./squeeze(sum(cx)); end
    cm(1,:,:) = 1 - squeeze(sum(cm(2:end,:,:)));
    cx(1,:,:) = 1 - squeeze(sum(cx(2:end,:,:)));

end

% get residual of thermochemical equations from iterative update
resnorm_TC = norm( T - Ti,2)./(norm(T,2)+TINY) ...
           + norm((x - xi).*(x>10*TINY),2)./(norm(x,2)+TINY);
