% ***** TRACE & ISOTOPE GEOCHEMISTRY  *************************************

% *****  Incompatible Trace Element  **************************************

% update incompatible trace element phase compositions
itm = it./(m + x.*KIT);
itx = it./(m./KIT + x);

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
ctm = ct./(m + x.*KCT);
ctx = ct./(m./KCT + x);

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

% update stable isotope in melt
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
ripm = rip./(m + x.*KRIP);
ripx = rip./(m./KRIP + x);

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
ridm = rid./(m + x.*KRID);
ridx = rid./(m./KRID + x);

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
si  = SI./rho;
rip = RIP./rho;
rid = RID./rho;