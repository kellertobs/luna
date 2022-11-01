%% *****  THERMO-CHEMICAL EVOLUTION  **************************************

% store previous iteration
Ti = T; ci = c; xi = x;

% update temperature    
advn_S = - advect(rho(inz,inx).*m(inz,inx).*sm(inz,inx),Um(inz,:),Wm(:,inx),h,{ADVN,''},[1,2],BCA) ...  % heat advection
         - advect(rho(inz,inx).*x(inz,inx).*sx(inz,inx),Ux(inz,:),Wx(:,inx),h,{ADVN,''},[1,2],BCA);

qTz    = - (ks(1:end-1,:)+ks(2:end,:))./2 .* ddz(T,h);                     % heat diffusion
qTx    = - (ks(:,1:end-1)+ks(:,2:end))./2 .* ddx(T,h);
diff_T = (- ddz(qTz(:,inx),h)  ...
          - ddx(qTx(inz,:),h));

diss_h = diss ./ T(inz,inx);

bnd_T = zeros(size(T(inz,inx)));
if ~isnan(Ttop); bnd_T = bnd_T + ((Ttop+273.15)-T(inz,inx))./tau_T .* topshape; end % impose top boundary layer
if ~isnan(Tbot); bnd_T = bnd_T + ((Tbot+273.15)-T(inz,inx))./tau_T .* botshape; end % impose bot boundary layer
bnd_S = rho(inz,inx).*cP.*bnd_T./T(inz,inx);

dSdt = advn_S + diff_T + diss_h + bnd_S;                                   % total rate of change

S(inz,inx) = So(inz,inx) + (theta.*dSdt + (1-theta).*dSdto).*dt;           % explicit update of major component density
S([1 end],:) = S([2 end-1],:);                                             % apply zero flux boundary conditions
S(:,[1 end]) = S(:,[2 end-1]);

% convert entropy to temperature
T = (T0+273.15) * exp(S./rho./cP - x.*Dsx./cP + aT./rhoref./cP.*(Pt-Ptop));% convert entropy to temperature


% update major component 
advn_C = zeros(size(c(inz,inx,:)));
diff_C = zeros(size(c(inz,inx,:)));
for i = 1:cal.nc
    advn_C(:,:,i) = - advect(rho(inz,inx).*m(inz,inx).*cm(inz,inx,i),Um(inz,:),Wm(:,inx),h,{ADVN,''},[1,2],BCA) ...
                    - advect(rho(inz,inx).*x(inz,inx).*cx(inz,inx,i),Ux(inz,:),Wx(:,inx),h,{ADVN,''},[1,2],BCA);

    qCz = - (kc(1:end-1,:)+kc(2:end,:))./2 .* ddz(c(:,:,i),h);             % component diffusion
    qCx = - (kc(:,1:end-1)+kc(:,2:end))./2 .* ddx(c(:,:,i),h);
    diff_C(:,:,i) = (- ddz(qCz(:,inx),h)  ...
                     - ddx(qCx(inz,:),h));
end

% get total rate of change
dCdt = advn_C + diff_C;

% update trace element concentrations
C(inz,inx,:) = Co(inz,inx,:) + (theta.*dCdt + (1-theta).*dCdto).*dt;       % explicit update
C = max(0,min(rho, C ));                                                   % enforce min bound
C([1 end],:,:) = C([2 end-1],:,:);                                         % boundary conditions
C(:,[1 end],:) = C(:,[2 end-1],:);

% convert component densities to concentrations
C(:,:,1) = max(0,min(rho, rho - sum(C(:,:,2:end),3) ));
for i = 1:cal.nc; c(:,:,i) = C(:,:,i)./rho; end


%% *****  UPDATE PHASE PROPORTIONS  ***************************************

% update local phase equilibrium
for i = 1:cal.nc
    var.c(:,i) = reshape(c(:,:,i),Nx*Nz,1);  % in wt
end
var.T = T(:)-273.15;  % convert to C
var.P = Pt(:)/1e9;    % convert to GPa
var.f = 1-xq(:);      % in wt

[phs,cal]  =  meltmodel(var,cal,'E'); % cs and cl component prop in each phase VS T

mq = max(TINY,min(1-TINY,reshape(phs.f ,Nz,Nx)));
xq = 1-mq;

for i = 2:cal.nc
    cxq(:,:,i) = reshape(phs.cs(:,i),Nz,Nx);
    cmq(:,:,i) = reshape(phs.cl(:,i),Nz,Nx);
end
cmq(:,:,1) = 1 - sum(cmq(:,:,2:end),3);
cxq(:,:,1) = 1 - sum(cxq(:,:,2:end),3);

% update crystal fraction
if diseq % quasi-equilibrium approach
    
    Gx = lambda.*Gx + (1-lambda) .* (xq-x).*rho./(4*dt);
        
    advn_X = - advect(rho(inz,inx).*x(inz,inx),Ux(inz,:),Wx(:,inx),h,{ADVN,''},[1,2],BCA);

    dXdt   = advn_X + Gx(inz,inx);                                         % total rate of change
    
    X(inz,inx) = Xo(inz,inx) + (theta.*dXdt + (1-theta).*dXdto).*dt;       % explicit update of crystal fraction
    X = min(rho.*(1-TINY),max(rho.*TINY,X));                               % enforce [0,1] limit
    X([1 end],:) = X([2 end-1],:);                                         % apply boundary conditions
    X(:,[1 end]) = X(:,[2 end-1]);

else
    
    X  =  lambda.*X + (1-lambda).*xq.*rho;
    Gx(inz,inx) = (X(inz,inx)-Xo(inz,inx))./dt + advect(rho(inz,inx).*x(inz,inx),Ux(inz,:),Wx(:,inx),h,{ADVN,''},[1,2],BCA);     % reconstruct crystallisation rate
    
end

% update phase fractions
x = X./rho;
m = max(0,min(1,1-x));

% phase entropies
sm = S./rho - x.*Dsx;
sx = sm + Dsx;

% phase compositions
Kc = cxq./cmq;
for i = 2:cal.nc
    % update trace element phase compositions
    cm(:,:,i)  = c(:,:,i)./(m + x.*Kc(:,:,i) );
    cx(:,:,i)  = c(:,:,i)./(m./Kc(:,:,i)  + x);
end
cm(:,:,1) = 1 - sum(cm(:,:,2:end),3);
cx(:,:,1) = 1 - sum(cx(:,:,2:end),3);

% get residual of thermochemical equations from iterative update
resnorm_TC = norm( T(:) - Ti(:),2)./(norm(T(:),2)+TINY) ...
           + norm( c(:) - ci(:),2)./(norm(c(:),2)+TINY) ...
           + norm((x(:) - xi(:)).*(x(:)>10*TINY),2)./(norm(x(:),2)+TINY);
