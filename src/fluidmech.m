%% assemble coefficients for matrix velocity diagonal and right-hand side
IIL  = [];       % equation indeces into L
JJL  = [];       % variable indeces into L
AAL  = [];       % coefficients for L
IIR  = [];       % equation indeces into R
AAR  = [];       % forcing entries for R


% assemble coefficients of z-stress divergence
    
% left boundary
ii = MapW(:,1); jj1 = ii; jj2 = MapW(:,2);
aa = zeros(size(ii));
IIL = [IIL; ii(:)]; JJL = [JJL; jj1(:)];   AAL = [AAL; aa(:)+1];
IIL = [IIL; ii(:)]; JJL = [JJL; jj2(:)];   AAL = [AAL; aa(:)+sds.*double(Nx<=10)];
IIR = [IIR; ii(:)]; AAR = [AAR; aa(:)];

% right boundary
ii = MapW(:,end); jj1 = ii; jj2 = MapW(:,end-1);
aa = zeros(size(ii));
IIL = [IIL; ii(:)]; JJL = [JJL; jj1(:)];   AAL = [AAL; aa(:)+1];
IIL = [IIL; ii(:)]; JJL = [JJL; jj2(:)];   AAL = [AAL; aa(:)+sds.*double(Nx<=10)];
IIR = [IIR; ii(:)]; AAR = [AAR; aa(:)];

% top boundary
ii = MapW(1,2:end-1); jj = ii;
aa = zeros(size(ii));
IIL = [IIL; ii(:)]; JJL = [JJL; jj(:)];   AAL = [AAL; aa(:)+1];
aa = zeros(size(ii)) + WBG(1,2:end-1);
IIR = [IIR; ii(:)]; AAR = [AAR; aa(:)];

% bottom boundary
ii = MapW(end,2:end-1); jj = ii;
aa = zeros(size(ii));
IIL = [IIL; ii(:)]; JJL = [JJL; jj(:)];   AAL = [AAL; aa(:)+1];
aa = zeros(size(ii)) + WBG(end,2:end-1);
IIR = [IIR; ii(:)]; AAR = [AAR; aa(:)];


% internal points
ii    = MapW(2:end-1,2:end-1);
EtaC1 = etaco(2:end-1,1:end-1); EtaC2 = etaco(2:end-1,2:end  );
EtaP1 = eta  (2:end-2,2:end-1); EtaP2 = eta  (3:end-1,2:end-1);

% coefficients multiplying z-velocities W
%             top          ||         bottom          ||           left            ||          right
jj1 = MapW(1:end-2,2:end-1); jj2 = MapW(3:end,2:end-1); jj3 = MapW(2:end-1,1:end-2); jj4 = MapW(2:end-1,3:end);

aa = - 2/3*(EtaP1+EtaP2)/h^2 - 1/2*(EtaC1+EtaC2)/h^2;
IIL = [IIL; ii(:)]; JJL = [JJL;  ii(:)];   AAL = [AAL; aa(:)           ];      % W on stencil centre
IIL = [IIL; ii(:)]; JJL = [JJL; jj1(:)];   AAL = [AAL; 2/3*EtaP1(:)/h^2];      % W one above
IIL = [IIL; ii(:)]; JJL = [JJL; jj2(:)];   AAL = [AAL; 2/3*EtaP2(:)/h^2];      % W one below
IIL = [IIL; ii(:)]; JJL = [JJL; jj3(:)];   AAL = [AAL; 1/2*EtaC1(:)/h^2];      % W one to the left
IIL = [IIL; ii(:)]; JJL = [JJL; jj4(:)];   AAL = [AAL; 1/2*EtaC2(:)/h^2];      % W one to the right

% what shall we do with a drunken sailor...
aa = -ddz(rho(2:end-1,2:end-1),h).*g0.*dt/2;
IIL = [IIL; ii(:)]; JJL = [JJL;  ii(:)];   AAL = [AAL; aa(:)];
    
% coefficients multiplying x-velocities U
%         top left         ||        bottom left          ||       top right       ||       bottom right
jj1 = MapU(2:end-2,1:end-1); jj2 = MapU(3:end-1,1:end-1); jj3 = MapU(2:end-2,2:end); jj4 = MapU(3:end-1,2:end);

IIL = [IIL; ii(:)]; JJL = [JJL; jj1(:)];   AAL = [AAL; (1/2*EtaC1(:)-1/3*EtaP1(:))/h^2];  % W one to the top and left
IIL = [IIL; ii(:)]; JJL = [JJL; jj2(:)];   AAL = [AAL;-(1/2*EtaC1(:)-1/3*EtaP2(:))/h^2];  % W one to the bottom and left
IIL = [IIL; ii(:)]; JJL = [JJL; jj3(:)];   AAL = [AAL;-(1/2*EtaC2(:)-1/3*EtaP1(:))/h^2];  % W one to the top and right
IIL = [IIL; ii(:)]; JJL = [JJL; jj4(:)];   AAL = [AAL; (1/2*EtaC2(:)-1/3*EtaP2(:))/h^2];  % W one to the bottom and right


% z-RHS vector

rr = - ((rho(2:end-2,2:end-1)+rho(3:end-1,2:end-1))/2 - rhoref) .* g0;
if bnchm; rr = rr + src_W_mms(2:end-1,2:end-1); end

IIR = [IIR; ii(:)];  AAR = [AAR; rr(:)];


%  assemble coefficients of x-stress divergence

% top boundary
ii = MapU(1,:); jj1 = ii; jj2 = MapU(2,:);
aa = zeros(size(ii));
IIL = [IIL; ii(:)]; JJL = [JJL; jj1(:)];   AAL = [AAL; aa(:)+1];
IIL = [IIL; ii(:)]; JJL = [JJL; jj2(:)];   AAL = [AAL; aa(:)+top];
IIR = [IIR; ii(:)]; AAR = [AAR; aa(:)];

% bottom boundary
ii = MapU(end,:); jj1 = ii; jj2 = MapU(end-1,:);
aa = zeros(size(ii));
IIL = [IIL; ii(:)]; JJL = [JJL; jj1(:)];   AAL = [AAL; aa(:)+1];
IIL = [IIL; ii(:)]; JJL = [JJL; jj2(:)];   AAL = [AAL; aa(:)+bot];
IIR = [IIR; ii(:)]; AAR = [AAR; aa(:)];

% left side boundary
ii = MapU(2:end-1,1); jj = ii;
aa = zeros(size(ii));
IIL = [IIL; ii(:)]; JJL = [JJL; jj(:)];   AAL = [AAL; aa(:)+1];
aa = zeros(size(ii)) + UBG(2:end-1,1);
IIR = [IIR; ii(:)]; AAR = [AAR; aa(:)];

% right side boundary
ii = MapU(2:end-1,end); jj = ii;
aa = zeros(size(ii));
IIL = [IIL; ii(:)]; JJL = [JJL; jj(:)];   AAL = [AAL; aa(:)+1];
aa = zeros(size(ii)) + UBG(2:end-1,end);
IIR = [IIR; ii(:)]; AAR = [AAR; aa(:)];


% internal points
ii    = MapU(2:end-1,2:end-1);
EtaC1 = etaco(1:end-1,2:end-1); EtaC2 = etaco(2:end  ,2:end-1);
EtaP1 = eta  (2:end-1,2:end-2); EtaP2 = eta  (2:end-1,3:end-1);

% coefficients multiplying x-velocities U
%            left          ||          right          ||           top             ||          bottom
jj1 = MapU(2:end-1,1:end-2); jj2 = MapU(2:end-1,3:end); jj3 = MapU(1:end-2,2:end-1); jj4 = MapU(3:end,2:end-1);

aa = - 2/3*(EtaP1+EtaP2)/h^2 - 1/2*(EtaC1+EtaC2)/h^2;
IIL = [IIL; ii(:)]; JJL = [JJL;  ii(:)];   AAL = [AAL; aa(:)           ];      % U on stencil centre
IIL = [IIL; ii(:)]; JJL = [JJL; jj1(:)];   AAL = [AAL; 2/3*EtaP1(:)/h^2];      % U one to the left
IIL = [IIL; ii(:)]; JJL = [JJL; jj2(:)];   AAL = [AAL; 2/3*EtaP2(:)/h^2];      % U one to the right
IIL = [IIL; ii(:)]; JJL = [JJL; jj3(:)];   AAL = [AAL; 1/2*EtaC1(:)/h^2];      % U one above
IIL = [IIL; ii(:)]; JJL = [JJL; jj4(:)];   AAL = [AAL; 1/2*EtaC2(:)/h^2];      % U one below

% coefficients multiplying z-velocities W
%         top left         ||        top right          ||       bottom left       ||       bottom right
jj1 = MapW(1:end-1,2:end-2); jj2 = MapW(1:end-1,3:end-1); jj3 = MapW(2:end,2:end-2); jj4 = MapW(2:end,3:end-1);

IIL = [IIL; ii(:)]; JJL = [JJL; jj1(:)];   AAL = [AAL; (1/2*EtaC1(:)-1/3*EtaP1(:))/h^2];  % W one to the top and left
IIL = [IIL; ii(:)]; JJL = [JJL; jj2(:)];   AAL = [AAL;-(1/2*EtaC1(:)-1/3*EtaP2(:))/h^2];  % W one to the top and right
IIL = [IIL; ii(:)]; JJL = [JJL; jj3(:)];   AAL = [AAL;-(1/2*EtaC2(:)-1/3*EtaP1(:))/h^2];  % W one to the bottom and left
IIL = [IIL; ii(:)]; JJL = [JJL; jj4(:)];   AAL = [AAL; (1/2*EtaC2(:)-1/3*EtaP2(:))/h^2];  % W one to the bottom and right


% x-RHS vector
rr = zeros(size(ii));
if bnchm; rr = rr + src_U_mms(2:end-1,2:end-1); end

IIR = [IIR; ii(:)];  AAR = [AAR; rr(:)];


% assemble coefficient matrix & right-hand side vector
KV = sparse(IIL,JJL,AAL,NW+NU,NW+NU);
RV = sparse(IIR,ones(size(IIR)),AAR);



%% assemble coefficients for gradient operator
if ~exist('GG','var') || bnchm
    IIL  = [];       % equation indeces into A
    JJL  = [];       % variable indeces into A
    AAL  = [];       % coefficients for A
    
    
    % coefficients for z-gradient
    ii = MapW(2:end-1,2:end-1);
    
    %         top              ||          bottom
    jj1 = MapP(2:end-2,2:end-1); jj2 = MapP(3:end-1,2:end-1);
    
    aa = zeros(size(ii));
    IIL = [IIL; ii(:)]; JJL = [JJL; jj1(:)];   AAL = [AAL; aa(:)-1/h];     % one to the top
    IIL = [IIL; ii(:)]; JJL = [JJL; jj2(:)];   AAL = [AAL; aa(:)+1/h];     % one to the bottom
    
    
    % coefficients for x-gradient
    ii = MapU(2:end-1,2:end-1);
    
    %         left             ||           right
    jj1 = MapP(2:end-1,2:end-2); jj2 = MapP(2:end-1,3:end-1);
    
    aa = zeros(size(ii));
    IIL = [IIL; ii(:)]; JJL = [JJL; jj1(:)];   AAL = [AAL; aa(:)-1/h];     % one to the left
    IIL = [IIL; ii(:)]; JJL = [JJL; jj2(:)];   AAL = [AAL; aa(:)+1/h];     % one to the right
    
    
    % assemble coefficient matrix
    GG = sparse(IIL,JJL,AAL,NW+NU,NP);
end


%% assemble coefficients for divergence operator
if ~exist('DD','var') || bnchm
    IIL  = [];       % equation indeces into A
    JJL  = [];       % variable indeces into A
    AAL  = [];       % coefficients for A
    
    %internal points
    ii = MapP(2:end-1,2:end-1);
    
    % coefficients multiplying velocities U, W
    %          left U          ||           right U       ||           top W           ||          bottom W
    jj1 = MapU(2:end-1,1:end-1); jj2 = MapU(2:end-1,2:end); jj3 = MapW(1:end-1,2:end-1); jj4 = MapW(2:end,2:end-1);
    
    aa = zeros(size(ii));
    IIL = [IIL; ii(:)]; JJL = [JJL; jj1(:)];   AAL = [AAL; aa(:)-1/h];  % U one to the left
    IIL = [IIL; ii(:)]; JJL = [JJL; jj2(:)];   AAL = [AAL; aa(:)+1/h];  % U one to the right
    IIL = [IIL; ii(:)]; JJL = [JJL; jj3(:)];   AAL = [AAL; aa(:)-1/h];  % W one above
    IIL = [IIL; ii(:)]; JJL = [JJL; jj4(:)];   AAL = [AAL; aa(:)+1/h];  % W one below
    
    % assemble coefficient matrix
    DD = sparse(IIL,JJL,AAL,NP,NW+NU);
end


%% assemble coefficients for matrix pressure diagonal and right-hand side
if ~exist('KP','var') || bnchm
    IIL  = [];       % equation indeces into A
    JJL  = [];       % variable indeces into A
    AAL  = [];       % coefficients for A
    
    % boundary points
    ii  = [MapP(1,:).'; MapP(end  ,:).']; % top & bottom
    jj1 = ii;
    jj2 = [MapP(2,:).'; MapP(end-1,:).'];
    
    aa = zeros(size(ii));
    IIL = [IIL; ii(:)]; JJL = [JJL; jj1(:)];   AAL = [AAL; aa(:)+1];
    IIL = [IIL; ii(:)]; JJL = [JJL; jj2(:)];   AAL = [AAL; aa(:)-1];
    
    ii  = [MapP(:,1); MapP(:,end  )]; % left & right
    jj1 = ii;
    jj2 = [MapP(:,2); MapP(:,end-1)];
    
    aa = zeros(size(ii));
    IIL = [IIL; ii(:)]; JJL = [JJL; jj1(:)];   AAL = [AAL; aa(:)+1];
    IIL = [IIL; ii(:)]; JJL = [JJL; jj2(:)];   AAL = [AAL; aa(:)-1];
    
    % internal points
    ii = MapP(2:end-1,2:end-1);
    
    % coefficients multiplying matrix pressure P
    aa = zeros(size(ii));
    IIL = [IIL; ii(:)]; JJL = [JJL; ii(:)];    AAL = [AAL; aa(:)];  % P on stencil centre
    
    % assemble coefficient matrix
    KP = sparse(IIL,JJL,AAL,NP,NP);
end

% RHS
IIR  = [];       % equation indeces into R
AAR  = [];       % forcing entries for R

ii = MapP(2:end-1,2:end-1);

rr = - VolSrc(2:end-1,2:end-1);
if bnchm; rr = rr + src_P_mms(2:end-1,2:end-1); end

IIR = [IIR; ii(:)]; AAR = [AAR; rr(:)];

% assemble right-hand side vector
RP = sparse(IIR,ones(size(IIR)),AAR,NP,1);

% get pressure scaling factor
Pscale = sqrt(geomean(eta(:))/h^2);

nzp = round((Nz-2)/2)+1;
nxp = round((Nx-2)/2)+1;
DD(MapP(nzp,nxp),:) = 0;
KP(MapP(nzp,nxp),:) = 0;
KP(MapP(nzp,nxp),MapP(nzp,nxp)) = 1;
RP(MapP(nzp,nxp),:) = 0;
if bnchm; RP(MapP(nzp,nxp),:) = P_mms(nzp,nxp); end


%% assemble global coefficient matrix and right-hand side vector
LL = [ KV  -GG  ; ...
      -DD   KP ];

RR = [RV; RP];


%% get residual
% get non-linear residual
FF         = LL*SS - RR;
resnorm_VP = norm(FF(:),2)./(norm(RR(:),2)+TINY);

% map residual vector to 2D arrays
res_W  = full(reshape(FF(MapW(:))        ,(Nz-1), Nx   ));  % z-velocity
res_U  = full(reshape(FF(MapU(:))        , Nz   ,(Nx-1)));  % x-velocity
res_P  = full(reshape(FF(MapP(:)+(NW+NU)), Nz   , Nx   ));  % dynamic pressure
        

%% Solve linear system of equations for vx, vz, P
SC = sqrt(abs(diag(LL)));
SC = diag(sparse(1./(SC+1)));

LL = SC*LL*SC;
RR = SC*RR;

SS = SC*(LL\RR);  % update solution


% map solution vector to 2D arrays
W  = full(reshape(SS(MapW(:))        ,(Nz-1), Nx   ));                     % z-velocity
U  = full(reshape(SS(MapU(:))        , Nz   ,(Nx-1)));                     % x-velocity
P  = full(reshape(SS(MapP(:)+(NW+NU)), Nz   , Nx   ));                     % dynamic pressure
% Pt = P + rhoref.*g0.*ZZ + Ptop;                                            % total pressure

% update phase velocities
Wx   = W + wx;                                                             % xtl z-velocity
Ux   = U + 0.;                                                             % xtl x-velocity
Wm   = W + wm;                                                             % mlt z-velocity
Um   = U + 0.;                                                             % mlt x-velocity

Wbar = (m(1:end-1,:)+m(2:end,:))/2 .* Wm ...
     + (x(1:end-1,:)+x(2:end,:))/2 .* Wx;
Ubar = (m(:,1:end-1)+m(:,2:end))/2 .* Um ...
     + (x(:,1:end-1)+x(:,2:end))/2 .* Ux; 

 
%% update time step
dtk = min((h/2)^2./max([kT(:)./rho(:)./cP;kc./rho(:)]))/2;                    % diffusive time step size
dta = CFL*min(min(h/2/max(abs([Ux(:);Wx(:);Um(:);Wm(:)]+1e-16)))); % advective time step size
dt  = min([2*dto,dtmax,min(dtk,dta)]);                                     % physical time step size
