%% *****  THERMO-CHEMICAL EVOLUTION  **************************************

% store previous iteration
Ti = T; xi = x;

% update temperature
if entr_mth  % use entropy equation to evolve heat
    
    sm = s - x.*Dsx;                                                       % phase entropies
    sx = sm + Dsx;
    
    advn_S = - advection(rho.*m.*sm,Um,Wm,h,ADVN,'flx') ...                % heat advection
             - advection(rho.*x.*sx,Ux,Wx,h,ADVN,'flx');
    
    qTz    = - (kT(1:end-1,:)+kT(2:end,:))/2 .* ddz(T,h);                  % heat diffusion z-flux
    qTx    = - (kT(:,1:end-1)+kT(:,2:end))/2 .* ddx(T,h);                  % heat diffusion x-flux
    diff_T(2:end-1,2:end-1) = (- ddz(qTz(:,2:end-1),h)  ...                % heat diffusion
                               - ddx(qTx(2:end-1,:),h));
    
    diss_h(2:end-1,2:end-1) = diss ./ T(2:end-1,2:end-1);
                        
    bndT = zeros(size(T));
    if ~isnan(Ttop); bndT = bndT + ((Ttop+273.15)-T)./tau_T .* topshape; end % impose top boundary layer
    if ~isnan(Tbot); bndT = bndT + ((Tbot+273.15)-T)./tau_T .* botshape; end % impose bot boundary layer
    bndS = rho.*cP.*bndT./T;
    
    dSdt = advn_S + diff_T + diss_h + bndS;                                % total rate of change
    
    S = So + (THETA.*dSdt + (1-THETA).*dSdto).*dt;                         % explicit update of entropy
    S([1 end],:) = S([2 end-1],:);                                         % apply zero flux boundary conditions
    S(:,[1 end]) = S(:,[2 end-1]);

    s = S./rho;
    
    T = (T0+273.15) * exp(S./rho./cP - x.*Dsx./cP + aT./rhoref./cP.*Pt);   % convert entropy to temperature
    
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
    
    T = To + (THETA.*dTdt + (1-THETA).*dTdto).*dt;                         % explicit update of temperature
    T([1 end],:) = T([2 end-1],:);                                         % apply zero-flux boundary conditions
    T(:,[1 end]) = T(:,[2 end-1]);

    S = rho.*(cP.*log(T/(T0+273.15)) + x.*Dsx - aT./rhoref.*Pt);                    % convert temperature to entropy
    diss_h(2:end-1,2:end-1) = diss ./ T(2:end-1,2:end-1);                  % convert dissipation rate to entropy units

end

% update major component 
for i = 2:cal.nc
    advn_C = advection(rho.*m.*squeeze(cm(i,:,:)),Um,Wm,h,ADVN,'flx') ...
           + advection(rho.*x.*squeeze(cx(i,:,:)),Ux,Wx,h,ADVN,'flx');
    
    qcz   = - kc.*(m(1:end-1,:)+m(2:end,:))/2 .* ddz(squeeze(c(i,:,:)),h); % major component diffusion z-flux
    qcx   = - kc.*(m(:,1:end-1)+m(:,2:end))/2 .* ddx(squeeze(c(i,:,:)),h); % major component diffusion x-flux
    diff_c(2:end-1,2:end-1) = - ddz(qcz(:,2:end-1),h) ...                  % major component diffusion
                              - ddx(qcx(2:end-1,:),h);
    
    dCdt(i,:,:) = - advn_C + diff_c;                                       % total rate of change
    
    C(i,:,:) = Co(i,:,:) + (THETA.*dCdt(i,:,:) + (1-THETA).*dCdto(i,:,:)).*dt;  % explicit update of major component density
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
if react
    for i = 1:cal.nc
        var.c(:,i) = reshape(c(i,:,:),Nx*Nz,1);  % in wt
    end
    var.T = T(:)-273.15;  % convert to C
    var.P = Pt(:)/1e9;    % convert to GPa
    var.f = 1-xq(:);      % in wt

    [phs,cal]  =  meltmodel(var,cal,'E'); % cs and cl component prop in each phase VS T 
    
    mq = reshape(phs.f ,Nz,Nx); 
    xq = 1-mq;
    for i = 2:cal.nc
        cxq(i,:,:) = reshape(phs.cs(:,i),Nz,Nx);
        cmq(i,:,:) = reshape(phs.cl(:,i),Nz,Nx);
    end
    cmq(1,:,:) = 1 - squeeze(sum(cmq(2:end,:,:)));
    cxq(1,:,:) = 1 - squeeze(sum(cxq(2:end,:,:)));
end

% update crystal fraction
if diseq || ~react  % disequilibrium approach
    
    if react
        Gx = ALPHA.*Gx + (1-ALPHA) .* (xq-x).*rho./max(4.*dt,tau_r);
    end
    
    advn_x = advection(rho.*x,Ux,Wx,h,ADVN,'flx');                         % get advection term
    
    dxdt   = - advn_x + Gx;                                                % total rate of change
    
    x = (rhoo.*xo + (THETA.*dxdt + (1-THETA).*dxdto).*dt)./rho;            % explicit update of crystal fraction
    x = min(1-TINY,max(TINY,x));                                         % enforce [0,1] limit
    x([1 end],:) = x([2 end-1],:);                                         % apply boundary conditions
    x(:,[1 end]) = x(:,[2 end-1]);
    
else  % equilibrium approach
    
    x  =  ALPHA.*x + (1-ALPHA).*xq;
    Gx = (rho.*x-rhoo.*xo)./dt + advection(rho.*x,Ux,Wx,h,ADVN,'flx');     % reconstruct crystallisation rate
    
end

% update melt fraction
m = min(1-TINY,max(TINY,1-x));

% update phase compositions
if react && step>0
    
    % major component
    Kc = cxq./cmq;
    for i = 2:cal.nc
        cm(i,:,:) = squeeze(c(i,:,:))./(m + x.*squeeze(Kc(i,:,:)));
        cx(i,:,:) = squeeze(c(i,:,:))./(m./squeeze(Kc(i,:,:)) + x);
    end
    cm(1,:,:) = 1 - squeeze(sum(cm(2:end,:,:)));
    cx(1,:,:) = 1 - squeeze(sum(cx(2:end,:,:)));

end

% get residual of thermochemical equations from iterative update
resnorm_TC = norm( T - Ti,2)./(norm(T,2)+TINY) ...
           + norm((x - xi).*(x>10*TINY),2)./(norm(x,2)+TINY);


%% ***** TRACE & ISOTOPE GEOCHEMISTRY  ************************************

% *****  Incompatible Trace Element  **************************************

% update incompatible trace element phase compositions
if react
    itm = it./(m + x.*KIT);
    itx = it./(m./KIT + x);
end

% update incompatible trace element composition
advn_IT = advection(rho.*m.*itm,Um,Wm,h,ADVN,'flx') ...
        + advection(rho.*x.*itx,Ux,Wx,h,ADVN,'flx');

qz   = - kc.*(m(1:end-1,:)+m(2:end,:))/2 .* ddz(it,h);
qx   = - kc.*(m(:,1:end-1)+m(:,2:end))/2 .* ddx(it,h);
diff_it(2:end-1,2:end-1) = - ddz(qz(:,2:end-1),h) ...                      % diffusion in melt
                           - ddx(qx(2:end-1,:),h);

dITdt = - advn_IT + diff_it;                                               % total rate of change

IT = ITo + (THETA.*dITdt + (1-THETA).*dITdto).*dt;                         % explicit update
IT = max(0+TINY, IT );
IT([1 end],:) = IT([2 end-1],:);                                           % boundary conditions
IT(:,[1 end]) = IT(:,[2 end-1]);


% *****  COMPATIBLE TRACE ELEMENT  ****************************************

% update compatible trace element phase compositions
if react
    ctm = ct./(m + x.*KCT);
    ctx = ct./(m./KCT + x);
end

% update compatible trace element composition
advn_CT = advection(rho.*m.*ctm,Um,Wm,h,ADVN,'flx') ...
        + advection(rho.*x.*ctx,Ux,Wx,h,ADVN,'flx');

qz   = - kc.*(m(1:end-1,:)+m(2:end,:))/2 .* ddz(ct,h);
qx   = - kc.*(m(:,1:end-1)+m(:,2:end))/2 .* ddx(ct,h);
diff_ct(2:end-1,2:end-1) = - ddz(qz(:,2:end-1),h) ...                      % diffusion in melt
                           - ddx(qx(2:end-1,:),h);

dCTdt = - advn_CT + diff_ct;                                               % total rate of change

CT = CTo + (THETA.*dCTdt + (1-THETA).*dCTdto).*dt;                         % explicit update
CT = max(0+TINY, CT );
CT([1 end],:) = CT([2 end-1],:);                                           % boundary conditions
CT(:,[1 end]) = CT(:,[2 end-1]);


% *****  STABLE ISOTOPE RATIO  ********************************************

% update stable isotope ratio in melt
advn_si = advection(rho.*m.*si,Um,Wm,h,ADVN,'flx') ...
        + advection(rho.*x.*si,Ux,Wx,h,ADVN,'flx');

qz   = - kc.*(m(1:end-1,:)+m(2:end,:))/2 .* ddz(si,h);
qx   = - kc.*(m(:,1:end-1)+m(:,2:end))/2 .* ddx(si,h);
diff_si(2:end-1,2:end-1) = - ddz(qz(:,2:end-1),h) ...                      % diffusion in melt
                           - ddx(qx(2:end-1,:),h);
                       
dSIdt = - advn_si + diff_si;                                               % total rate of change

SI = SIo + (THETA.*dSIdt + (1-THETA).*dSIdto).*dt;                         % explicit update
SI([1 end],:) = SI([2 end-1],:);                                           % boundary conditions
SI(:,[1 end]) = SI(:,[2 end-1]);


% *****  RADIOGENIC ISOTOPES  *********************************************

% decay rate of radiogenic isotope
dcy_rip = rho.*rip./HLRIP.*log(2);
dcy_rid = rho.*rid./HLRID.*log(2);

% update radiogenic parent isotope phase compositions
if react
    ripm = rip./(m + x.*KRIP);
    ripx = rip./(m./KRIP + x);
end

% update radiogenic parent isotope composition
advn_RIP = advection(rho.*m.*ripm,Um,Wm,h,ADVN,'flx') ...
         + advection(rho.*x.*ripx,Ux,Wx,h,ADVN,'flx');

qz   = - kc.*(m(1:end-1,:)+m(2:end,:))/2 .* ddz(rip,h);
qx   = - kc.*(m(:,1:end-1)+m(:,2:end))/2 .* ddx(rip,h);
diff_rip(2:end-1,2:end-1) = - ddz(qz(:,2:end-1),h) ...                     % diffusion in melt
                            - ddx(qx(2:end-1,:),h);

                                       % secular equilibrium!
dRIPdt = - advn_RIP + diff_rip - dcy_rip + dcy_rip;                        % total rate of change
                                       
RIP = RIPo + (THETA.*dRIPdt + (1-THETA).*dRIPdto).*dt;                     % explicit update
RIP = max(0+TINY, RIP );
RIP([1 end],:) = RIP([2 end-1],:);                                         % boundary conditions
RIP(:,[1 end]) = RIP(:,[2 end-1]);


% update radiogenic daughter isotope phase compositions
if react
    ridm = rid./(m + x.*KRID);
    ridx = rid./(m./KRID + x);
end

% update radiogenic daughter isotope composition
advn_RID = advection(rho.*m.*ridm,Um,Wm,h,ADVN,'flx') ...
         + advection(rho.*x.*ridx,Ux,Wx,h,ADVN,'flx');

qz   = - kc.*(m(1:end-1,:)+m(2:end,:))/2 .* ddz(rid,h);
qx   = - kc.*(m(:,1:end-1)+m(:,2:end))/2 .* ddx(rid,h);
diff_rid(2:end-1,2:end-1) = - ddz(qz(:,2:end-1),h) ...                     % diffusion in melt
                            - ddx(qx(2:end-1,:),h);
    
dRIDdt = - advn_RID + diff_rid - dcy_rid + dcy_rip;                        % total rate of change

RID = RIDo + (THETA.*dRIDdt + (1-THETA).*dRIDdto).*dt;                     % explicit update
RID = max(0+TINY, RID );
RID([1 end],:) = RID([2 end-1],:);                                         % boundary conditions
RID(:,[1 end]) = RID(:,[2 end-1]);

it  = IT./rho;
ct  = CT./rho;
si  = SI ./rho;
rip = RIP./rho;
rid = RID./rho;