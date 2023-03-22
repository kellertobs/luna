% record run history

dsumMdtoo = dsumMdto; dsumMdto = dsumMdt;
dsumSdtoo = dsumSdto; dsumSdto = dsumSdt;
dsumCdtoo = dsumCdto; dsumCdto = dsumCdt;

stp = max(1,step+1);

% record model time
hist.time(stp) = time;

% record total mass, heat, component mass in model (assume hy = 1, unit length in third dimension)
hist.sumM(stp) = sum(sum(rho*h*h*1));  % [kg]
hist.sumS(stp) = sum(sum((S+S0)*h*h*1));  % [J/K]
for i = 1:cal.ncmp; hist.sumC(stp,i) = sum(sum(  squeeze(C(:,:,i))*h*h*1));  end% [kg]

% record expected rates of change by volume change and imposed boundaries layers
dsumMdt = sum(  X(1,:).*Wx(1,2:end-1)*h*1) - sum(X(end,:).*Wx(end,2:end-1)*h*1) ...
        + sum(  X(:,1).*Ux(2:end-1,1)*h*1) - sum(X(:,end).*Ux(2:end-1,end)*h*1) ...
        + sum(  M(1,:).*Wm(1,2:end-1)*h*1) - sum(M(end,:).*Wm(end,2:end-1)*h*1) ...
        + sum(  M(:,1).*Um(2:end-1,1)*h*1) - sum(M(:,end).*Um(2:end-1,end)*h*1);  % [kg/s]
dsumSdt = sum(sum(bnd_S*h*h*1)) + sum(sum(diss_h*h*h*1)) ...
        + sum(  X(1,:).*sx(1,:).*Wx(1,2:end-1)*h*1) - sum(X(end,:).*sx(end,:).*Wx(end,2:end-1)*h*1) ...
        + sum(  X(:,1).*sx(:,1).*Ux(2:end-1,1)*h*1) - sum(X(:,end).*sx(:,end).*Ux(2:end-1,end)*h*1) ...
        + sum(  M(1,:).*sm(1,:).*Wm(1,2:end-1)*h*1) - sum(M(end,:).*sm(end,:).*Wm(end,2:end-1)*h*1) ...
        + sum(  M(:,1).*sm(:,1).*Um(2:end-1,1)*h*1) - sum(M(:,end).*sm(:,end).*Um(2:end-1,end)*h*1);  % [J /s]
dsumCdt = squeeze( ...
        + sum(  X(1,:).*cx(1,:,:).*Wx(1,2:end-1)*h*1,2) - sum(X(end,:).*cx(end,:,:).*Wx(end,2:end-1)*h*1,2) ...
        + sum(  X(:,1).*cx(:,1,:).*Ux(2:end-1,1)*h*1,1) - sum(X(:,end).*cx(:,end,:).*Ux(2:end-1,end)*h*1,1) ...
        + sum(  M(1,:).*cm(1,:,:).*Wm(1,2:end-1)*h*1,2) - sum(M(end,:).*cm(end,:,:).*Wm(end,2:end-1)*h*1,2) ...
        + sum(  M(:,1).*cm(:,1,:).*Um(2:end-1,1)*h*1,1) - sum(M(:,end).*cm(:,end,:).*Um(2:end-1,end)*h*1,1)).';  % [kg/s]

if step>=1; hist.DM(stp  ) = (a2*hist.DM(stp-1  ) + a3*hist.DM(max(1,stp-2)  ) + (b1*dsumMdt + b2*dsumMdto + b3*dsumMdtoo)*dt)/a1; else; hist.DM(stp  ) = 0; end  % [kg]
if step>=1; hist.DS(stp  ) = (a2*hist.DS(stp-1  ) + a3*hist.DS(max(1,stp-2)  ) + (b1*dsumSdt + b2*dsumSdto + b3*dsumSdtoo)*dt)/a1; else; hist.DS(stp  ) = 0; end  % [kg]
if step>=1; hist.DC(stp,:) = (a2*hist.DC(stp-1,:) + a3*hist.DC(max(1,stp-2),:) + (b1*dsumCdt + b2*dsumCdto + b3*dsumCdtoo)*dt)/a1; else; hist.DC(stp,:) = zeros(1,cal.ncmp); end  % [kg]

% record conservation error of mass M, heat S, major component C, volatile component V
hist.EM(stp)   = (hist.sumM(stp)   - hist.DM(stp))  ./hist.sumM(1)   - 1;  % [kg/kg]
hist.ES(stp)   = (hist.sumS(stp)   - hist.DS(stp))  ./hist.sumS(1)   - 1;  % [JK/JK]
hist.EC(stp,:) = (hist.sumC(stp,:) - hist.DC(stp,:))./hist.sumC(1,:) - 1;  % [kg/kg]

% record variable and coefficient diagnostics
hist.W(stp,1) = min(    -W(:,2:end-1),[],'all');
hist.W(stp,2) = mean(abs(W(:,2:end-1))  ,'all');
hist.W(stp,3) = max(    -W(:,2:end-1),[],'all');

hist.U(stp,1) = min(     U(2:end-1,:),[],'all');
hist.U(stp,2) = mean(abs(U(2:end-1,:))  ,'all');
hist.U(stp,3) = max(     U(2:end-1,:),[],'all');

hist.P(stp,1) = min(     P(2:end-1,2:end-1),[],'all');
hist.P(stp,2) = mean(abs(P(2:end-1,2:end-1))  ,'all');
hist.P(stp,3) = max(     P(2:end-1,2:end-1),[],'all');

hist.x(stp,1) = min( x,[],'all');
hist.x(stp,2) = mean(x   ,'all');
hist.x(stp,3) = max( x,[],'all');

hist.m(stp,1) = min( m,[],'all');
hist.m(stp,2) = mean(m   ,'all');
hist.m(stp,3) = max( m,[],'all');

hist.chi(stp,1) = min( chi,[],'all');
hist.chi(stp,2) = mean(chi   ,'all');
hist.chi(stp,3) = max( chi,[],'all');

hist.mu(stp,1) = min( mu,[],'all');
hist.mu(stp,2) = mean(mu   ,'all');
hist.mu(stp,3) = max( mu,[],'all');

hist.T(stp,1) = min( T,[],'all');
hist.T(stp,2) = mean(T   ,'all');
hist.T(stp,3) = max( T,[],'all');

for i = 1:cal.ncmp
    hist.c(stp,1,i) = min( c(:,:,i),[],'all');
    hist.c(stp,2,i) = mean(c(:,:,i)   ,'all');
    hist.c(stp,3,i) = max( c(:,:,i),[],'all');
end
for i=1:cal.noxd
    hist.c_oxd(stp,1,i) = min(min(c_oxd(:,:,i)));
    hist.c_oxd(stp,2,i) = mean(mean(c_oxd(:,:,i)));
    hist.c_oxd(stp,3,i) = max(max(c_oxd(:,:,i)));
end

indx     = x>1e-6;
indx_cmp = repmat(x>1e-6,1,1,cal.ncmp);
indx_oxd = repmat(x>1e-6,1,1,cal.noxd);
if any(indx(:)>0)
    for i=1:cal.ncmp
        hist.cx(stp,1,i) = min(min(cx(indx_cmp(:,:,i))));
        hist.cx(stp,2,i) = sum(sum(cx(:,:,i).*x.*rho))./sum(sum(x.*rho));
        hist.cx(stp,3,i) = max(max(cx(indx_cmp(:,:,i))));
    end
    for i=1:cal.noxd
        hist.cx_oxd(stp,1,i) = min(min(cx_oxd(indx_oxd(:,:,i))));
        hist.cx_oxd(stp,2,i) = sum(sum(cx_oxd(:,:,i).*x.*rho))./sum(sum(x.*rho));
        hist.cx_oxd(stp,3,i) = max(max(cx_oxd(indx_oxd(:,:,i))));
    end

    hist.rhox(stp,1) = min(min(rhox(indx)));
    hist.rhox(stp,2) = sum(sum(rhox.*x.*rho))./sum(sum(x.*rho));
    hist.rhox(stp,3) = max(max(rhox(indx)));
else
    hist.cx(stp,1:3,1:cal.ncmp) = NaN;
    hist.cx_oxd(stp,1:3,1:cal.noxd) = NaN;
    hist.rhox(stp,1:3) = NaN;
end

indm     = m>1e-6;
indm_cmp = repmat(m>1e-6,1,1,cal.ncmp);
indm_oxd = repmat(m>1e-6,1,1,cal.noxd);
if any(indm(:)>0)
    for i=1:cal.ncmp
        hist.cm(stp,1,i) = min(min(cm(indm_cmp(:,:,i))));
        hist.cm(stp,2,i) = sum(sum(cm(:,:,i).*m.*rho))./sum(sum(m.*rho));
        hist.cm(stp,3,i) = max(max(cm(indm_cmp(:,:,i))));
    end
    for i=1:cal.noxd
        hist.cm_oxd(stp,1,i) = min(min(cm_oxd(indm_oxd(:,:,i))));
        hist.cm_oxd(stp,2,i) = sum(sum(cm_oxd(:,:,i).*m.*rho))./sum(sum(m.*rho));
        hist.cm_oxd(stp,3,i) = max(max(cm_oxd(indm_oxd(:,:,i))));
    end

    hist.rhom(stp,1) = min(min(rhom(indm)));
    hist.rhom(stp,2) = sum(sum(rhom.*m))./sum(sum(m));
    hist.rhom(stp,3) = max(max(rhom(indm)));

    hist.etam(stp,1) = min(min(etam(indm)));
    hist.etam(stp,2) = sum(sum(etam.*m))./sum(sum(m));
    hist.etam(stp,3) = max(max(etam(indm)));
else
    hist.cm(stp,1:3,1:cal.ncmp) = NaN;
    hist.cm_oxd(stp,1:3,1:cal.noxd) = NaN;
    hist.rhom(stp,1:3) = NaN;
    hist.etam(stp,1:3) = NaN;
end

hist.Gx(stp,1) = min( Gx,[],'all');
hist.Gx(stp,2) = mean(Gx   ,'all');
hist.Gx(stp,3) = max( Gx,[],'all');

hist.dV(stp,1) = min( VolSrc,[],'all');
hist.dV(stp,2) = mean(VolSrc   ,'all');
hist.dV(stp,3) = max( VolSrc,[],'all');

hist.rho(stp,1) = min( rho,[],'all');
hist.rho(stp,2) = mean(rho   ,'all');
hist.rho(stp,3) = max( rho,[],'all');

hist.eta(stp,1) = min(    eta,[],'all');
hist.eta(stp,2) = geomean(eta   ,'all');
hist.eta(stp,3) = max(    eta,[],'all');