% ***** TRACE & ISOTOPE GEOCHEMISTRY  *************************************

% *****  Trace Elements  **************************************************

bnd_TE = zeros(Nz-2,Nx-2,cal.nte);
adv_TE = zeros(Nz-2,Nx-2,cal.nte);
Kte    = zeros(Nz,Nx,cal.nte);
for i = 1:cal.nte
    
    % update bulk partitioning coefficients
    for j=1:cal.nmnr; Kte(:,:,i) = Kte(:,:,i) + cal.Kte_cmp(i,j) .* cx_mnr(:,:,j)./100; end

    % update trace element phase compositions
    tem(:,:,i) = te(:,:,i)./(m + x.*Kte(:,:,i));
    tex(:,:,i) = te(:,:,i)./(m./Kte(:,:,i) + x);

    % get trace element advection
    adv_TE(:,:,i) = - advect(M(inz,inx).*tem(inz,inx,i),Um(inz,:),Wm(:,inx),h,{ADVN,''},[1,2],BCA) ...
                    - advect(X(inz,inx).*tex(inz,inx,i),Ux(inz,:),Wx(:,inx),h,{ADVN,''},[1,2],BCA);

end

% get total rate of change
dTEdt = adv_TE;

% update trace element concentrations
TE(inz,inx,:) = (a2*TEo(inz,inx,:) + a3*TEoo(inz,inx,:) + (b1*dTEdt + b2*dTEdto + b3*dTEdtoo)*dt)/a1;
TE = max(0, TE );                                                          % enforce min bound
TE([1 end],:,:) = TE([2 end-1],:,:);                                       % boundary conditions
TE(:,[1 end],:) = TE(:,[2 end-1],:);


% *****  Isotope Ratios  **************************************************

bnd_IR = zeros(Nz-2,Nx-2,cal.nir);
adv_IR = zeros(Nz-2,Nx-2,cal.nir);
for i = 1:cal.nir

    % update trace element phase compositions
    irm(:,:,i) = ir(:,:,i);
    irx(:,:,i) = ir(:,:,i);

    % get isotope ratio advection
    adv_IR(:,:,i) = - advect(M(inz,inx).*irm(inz,inx,i),Um(inz,:),Wm(:,inx),h,{ADVN,''},[1,2],BCA) ...
                    - advect(X(inz,inx).*irx(inz,inx,i),Ux(inz,:),Wx(:,inx),h,{ADVN,''},[1,2],BCA);

end

% get total rate of change
dIRdt = adv_IR;

% update isotope ratio concentrations
IR(inz,inx,:) = (a2*IRo(inz,inx,:) + a3*IRoo(inz,inx,:) + (b1*dIRdt + b2*dIRdto + b3*dIRdtoo)*dt)/a1;
IR([1 end],:,:) = IR([2 end-1],:,:);                                       % boundary conditions
IR(:,[1 end],:) = IR(:,[2 end-1],:);

% convert from mixture density to concentration
for i = 1:cal.nte; te(:,:,i) = TE(:,:,i)./rho; end
for i = 1:cal.nir; ir(:,:,i) = IR(:,:,i)./rho; end