%% *****  THERMO-CHEMICAL EVOLUTION  **************************************

% store previous iterationcal.ncmp
Si = S; Ci = C; Xi = X;

% update temperature    
advn_S = - advect(M(inz,inx).*sm(inz,inx),Um(inz,:),Wm(:,inx),h,{ADVN,''},[1,2],BCA) ...  % heat advection
         - advect(X(inz,inx).*sx(inz,inx),Ux(inz,:),Wx(:,inx),h,{ADVN,''},[1,2],BCA);

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

% semi-implicit update of bulk entropy density
S(inz,inx) = (a2*So(inz,inx) + a3*Soo(inz,inx) + (b1*dSdt + b2*dSdto + b3*dSdtoo)*dt)/a1; % semi-implicit update of bulk entropy density
S([1 end],:) = S([2 end-1],:);                                             % apply zero flux boundary conditions
S(:,[1 end]) = S(:,[2 end-1]);

% convert entropy to temperature
T = (T0+273.15) * exp((S - X.*Dsx)./rho./cP + Adbt./cP.*(Pt-Ptop));        % convert entropy to temperature


% update major component 
advn_C = zeros(size(c(inz,inx,:)));
diff_C = zeros(size(c(inz,inx,:)));
for i = 2:cal.ncmp
    advn_C(:,:,i) = - advect(M(inz,inx).*cm(inz,inx,i),Um(inz,:),Wm(:,inx),h,{ADVN,''},[1,2],BCA) ...  % melt advection
                    - advect(X(inz,inx).*cx(inz,inx,i),Ux(inz,:),Wx(:,inx),h,{ADVN,''},[1,2],BCA);     % xtal advection
end

% get total rate of change
dCdt = advn_C;

% semi-implicit update of major component density
C(inz,inx,:) = (a2*Co(inz,inx,:) + a3*Coo(inz,inx,:) + (b1*dCdt + b2*dCdto + b3*dCdtoo)*dt)/a1;
C = max(0, C );                                                            % enforce min bound
C([1 end],:,:) = C([2 end-1],:,:);                                         % boundary conditions
C(:,[1 end],:) = C(:,[2 end-1],:);

% convert component densities to concentrations
C(:,:,1) = max(0,min(rho, rho - sum(C(:,:,2:end),3) ));
for i = 1:cal.ncmp; c(:,:,i) = C(:,:,i)./rho; end


%% *****  UPDATE PHASE PROPORTIONS  ***************************************

% update local phase equilibrium
var.c = reshape(c(inz,inx,:),(Nx-2)*(Nz-2),cal.ncmp);  % in wt
var.T = reshape(T(inz,inx),(Nx-2)*(Nz-2),1)-273.15;  % convert to C
var.P = reshape(Pt(inz,inx),(Nx-2)*(Nz-2),1)/1e9;    % convert to GPa
var.f = 1-reshape(xq(inz,inx),(Nx-2)*(Nz-2),1);      % in wt

[phs,cal]  =  meltmodel(var,cal,'E'); % cs and cl component prop in each phase VS T

mq(inz,inx) = max(0,min(1,reshape(phs.f,(Nz-2),(Nx-2))));
mq([1 end],:) = mq([2 end-1],:);
mq(:,[1 end]) = mq(:,[2 end-1]);

xq = 1-mq;

for i = 2:cal.ncmp
    cxq(inz,inx,i) = reshape(phs.cs(:,i),(Nz-2),(Nx-2));
    cmq(inz,inx,i) = reshape(phs.cl(:,i),(Nz-2),(Nx-2));
end
cmq(:,:,1) = 1 - sum(cmq(:,:,2:end),3);
cxq(:,:,1) = 1 - sum(cxq(:,:,2:end),3);

cxq([1 end],:,:) = cxq([2 end-1],:,:);
cxq(:,[1 end],:) = cxq(:,[2 end-1],:);
cmq([1 end],:,:) = cmq([2 end-1],:,:);
cmq(:,[1 end],:) = cmq(:,[2 end-1],:);

% update crystal fraction
Gx = lambda * Gx + (1-lambda) * (rho(inz,inx).*xq(inz,inx)-X(inz,inx))./(3*dt);

advn_X = - advect(X(inz,inx),Ux(inz,:),Wx(:,inx),h,{ADVN,''},[1,2],BCA);   % xtal advection

dXdt   = advn_X + Gx;                                                      % total rate of change

% semi-implicit update of crystal fraction
X(inz,inx) = (a2*Xo(inz,inx) + a3*Xoo(inz,inx) + (b1*dXdt + b2*dXdto + b3*dXdtoo)*dt)/a1;
X = max(0, X );                                                            % enforce min bound
X([1 end],:) = X([2 end-1],:);                                             % apply boundary conditions
X(:,[1 end]) = X(:,[2 end-1]);

M = rho-X;

% update phase fractions
x = X./rho;
m = max(0,min(1,1-x));

% phase entropies
sm = (S - X.*Dsx)./rho;
sx = sm + Dsx;

% phase compositions
Kc = cxq./cmq;
for i = 2:cal.ncmp
    % update trace element phase compositions
    cm(:,:,i)  = c(:,:,i)./(m + x.*Kc(:,:,i) );
    cx(:,:,i)  = c(:,:,i)./(m./Kc(:,:,i)  + x);
end
cm(:,:,1) = 1 - sum(cm(:,:,2:end),3);
cx(:,:,1) = 1 - sum(cx(:,:,2:end),3);

% get residual of thermochemical equations from iterative update
resnorm_TC = norm( S(inz,inx) - Si(inz,inx),'fro')./(norm(S(inz,inx),'fro')+TINY) ...
           + norm( C(inz,inx) - Ci(inz,inx),'fro')./(norm(C(inz,inx),'fro')+TINY) ...
           + norm((X(inz,inx) - Xi(inz,inx)).*(x(inz,inx)>10*TINY),2)./(norm(X(inz,inx),2)+TINY);
