%% *****  THERMO-CHEMICAL EVOLUTION  **************************************

% store previous iteration
Ti = T; xi = x;

% update temperature
advn_H = advection(rho.*m.*T.*(Cp + 0  ),Um,Wm,h,ADVN,'flx') ...
       + advection(rho.*x.*T.*(Cp + Dsx),Ux,Wx,h,ADVN,'flx');
                           
qTz    = - kT .* ddz(T,h);                     % heat diffusion z-flux
qTx    = - kT .* ddx(T,h);                     % heat diffusion x-flux
diff_T(2:end-1,2:end-1) = (- ddz(qTz(:,2:end-1),h) ...                     % heat diffusion
                           - ddx(qTx(2:end-1,:),h));
    
bndH = zeros(size(T));
if ~isnan(Ttop); bndH = bndH + rho.*Cp.*(Ttop+273.15-T)./tau_T .* topshape; end % impose top boundary layer
if ~isnan(Tbot); bndH = bndH + rho.*Cp.*(Tbot+273.15-T)./tau_T .* botshape; end % impose bot boundary layer

dHdt = - advn_H + diff_T + bndH;                                           % total rate of change
    
H = Ho + (THETA.*dHdt + (1-THETA).*dHdto).*dt;                             % explicit update of enthalpy
H([1 end],:) = H([2 end-1],:);                                             % apply boundary conditions
H(:,[1 end]) = H(:,[2 end-1]);    
    
% update major component
for i = 1:cal.nc
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

% convert enthalpy and component densities to temperature and concentrations
T = H./(rho.*(Cp + Ds));
% C(1,:,:) = max(TINY,min(rho-TINY, rho - squeeze(sum(C(2:end,:,:))) ));
sumC = squeeze(sum(C));
for i = 1:cal.nc; c(i,:,:) = squeeze(C(i,:,:))./sumC; end
for i = 1:cal.nc; C(i,:,:) = squeeze(c(i,:,:)).*rho; end


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
    for i = 1:cal.nc
        cxq(i,:,:) = reshape(phs.cs(:,i),Nz,Nx);
        cmq(i,:,:) = reshape(phs.cl(:,i),Nz,Nx);
    end
%     cmq(1,:,:) = 1 - squeeze(sum(cmq(2:end,:,:)));
%     cxq(1,:,:) = 1 - squeeze(sum(cxq(2:end,:,:)));
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
    for i = 1:cal.nc
        cm(i,:,:) = squeeze(c(i,:,:))./(m + x.*squeeze(Kc(i,:,:)));
        cx(i,:,:) = squeeze(c(i,:,:))./(m./squeeze(Kc(i,:,:)) + x);
    end
%     cm(1,:,:) = 1 - squeeze(sum(cm(2:end,:,:)));
%     cx(1,:,:) = 1 - squeeze(sum(cx(2:end,:,:)));

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