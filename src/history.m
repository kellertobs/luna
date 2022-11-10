% record run history

dsumMdto = dsumMdt;
dsumSdto = dsumSdt;
dsumCdto = dsumCdt;

stp = max(1,step-1);

% record model time
hist.time(stp) = time;

% record total mass, heat, component mass in model (assume hy = 1, unit length in third dimension)
hist.sumM(stp) = sum(sum(rho(2:end-1,2:end-1)*h*h*1));  % [kg]
hist.sumS(stp) = sum(sum((S(2:end-1,2:end-1)+S0(2:end-1,2:end-1))*h*h*1));  % [J/K]
for i = 1:cal.nc; hist.sumC(stp,i) = sum(sum(  squeeze(C(2:end-1,2:end-1,i))*h*h*1));  end% [kg]

% record expected rates of change by volume change and imposed boundaries layers
dsumMdt = sum(rho(2,2:end-1).*W(1,2:end-1)*h*1) - sum(rho(end-1,2:end-1).*W(end,2:end-1)*h*1) ...
        + sum(rho(2:end-1,2).*U(2:end-1,1)*h*1) - sum(rho(2:end-1,end-1).*U(2:end-1,end)*h*1);  % [kg/s]
dsumSdt = sum(sum(bnd_S*h*h*1)) + sum(sum(diss_h*h*h*1))...
        + sum(  S(2,2:end-1).*W(1,2:end-1)*h*1) - sum(  S(end-1,2:end-1).*W(end,2:end-1)*h*1) ...
        + sum(  S(2:end-1,2).*U(2:end-1,1)*h*1) - sum(  S(2:end-1,end-1).*U(2:end-1,end)*h*1);  % [J/K/s]
for i = 1:cal.nc
    dsumCdt(i) = sum(  C(2,2:end-1,i).*W(1,2:end-1)*h*1) - sum(  C(end-1,2:end-1,i).*W(end,2:end-1)*h*1) ...
               + sum(  C(2:end-1,2,i).*U(2:end-1,1)*h*1) - sum(  C(2:end-1,end-1,i).*U(2:end-1,end)*h*1);  % [kg/s]
end

if stp>1; hist.DM(stp)   = hist.DM(stp-1)   +        dsumMdt                      .*dt; else; hist.DM(stp) = 0; end  % [kg]
if stp>1; hist.DS(stp)   = hist.DS(stp-1)   + (theta*dsumSdt + (1-theta)*dsumSdto).*dt; else; hist.DS(stp) = 0; end  % [J/K]
if stp>1; hist.DC(stp,:) = hist.DC(stp-1,:) + (theta*dsumCdt + (1-theta)*dsumCdto).*dt; else; hist.DC(stp,:) = zeros(1,cal.nc); end  % [kg]

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

hist.x(stp,1) = min( x(2:end-1,2:end-1),[],'all');
hist.x(stp,2) = mean(x(2:end-1,2:end-1)   ,'all');
hist.x(stp,3) = max( x(2:end-1,2:end-1),[],'all');

hist.m(stp,1) = min( m(2:end-1,2:end-1),[],'all');
hist.m(stp,2) = mean(m(2:end-1,2:end-1)   ,'all');
hist.m(stp,3) = max( m(2:end-1,2:end-1),[],'all');

hist.chi(stp,1) = min( chi(2:end-1,2:end-1),[],'all');
hist.chi(stp,2) = mean(chi(2:end-1,2:end-1)   ,'all');
hist.chi(stp,3) = max( chi(2:end-1,2:end-1),[],'all');

hist.mu(stp,1) = min( mu(2:end-1,2:end-1),[],'all');
hist.mu(stp,2) = mean(mu(2:end-1,2:end-1)   ,'all');
hist.mu(stp,3) = max( mu(2:end-1,2:end-1),[],'all');

hist.T(stp,1) = min( T(2:end-1,2:end-1),[],'all');
hist.T(stp,2) = mean(T(2:end-1,2:end-1)   ,'all');
hist.T(stp,3) = max( T(2:end-1,2:end-1),[],'all');

for i = 1:cal.nc
    hist.c(stp,1,i) = min( c(2:end-1,2:end-1,i),[],'all');
    hist.c(stp,2,i) = mean(c(2:end-1,2:end-1,i)   ,'all');
    hist.c(stp,3,i) = max( c(2:end-1,2:end-1,i),[],'all');
    hist.c_oxd(stp,1,i) = min( c_oxd(2:end-1,2:end-1,i),[],'all');
    hist.c_oxd(stp,2,i) = mean(c_oxd(2:end-1,2:end-1,i)   ,'all');
    hist.c_oxd(stp,3,i) = max( c_oxd(2:end-1,2:end-1,i),[],'all');
end

indx = repmat(x>1e-6,1,1,cal.nc);
if any(indx(:)>0)
    for i = 1:cal.nc
        hist.cx(stp,1,i) = min(min(cx(indx(2:end-1,2:end-1,i))));
        hist.cx(stp,2,i) = sum(sum(cx(2:end-1,2:end-1,i).*x(2:end-1,2:end-1).*rho(2:end-1,2:end-1)))./sum(sum(x(2:end-1,2:end-1).*rho(2:end-1,2:end-1)));
        hist.cx(stp,3,i) = max(max(cx(indx(2:end-1,2:end-1,i))));

        hist.cx_oxd(stp,1,i) = min(min(cx_oxd(indx(2:end-1,2:end-1,i))));
        hist.cx_oxd(stp,2,i) = sum(sum(cx_oxd(2:end-1,2:end-1,i).*x(2:end-1,2:end-1).*rho(2:end-1,2:end-1)))./sum(sum(x(2:end-1,2:end-1).*rho(2:end-1,2:end-1)));
        hist.cx_oxd(stp,3,i) = max(max(cx_oxd(indx(2:end-1,2:end-1,i))));
    end
    hist.rhox(stp,1) = min(min(rhox(indx(2:end-1,2:end-1,1))));
    hist.rhox(stp,2) = sum(sum(rhox(2:end-1,2:end-1).*x(2:end-1,2:end-1).*rho(2:end-1,2:end-1)))./sum(sum(x(2:end-1,2:end-1).*rho(2:end-1,2:end-1)));
    hist.rhox(stp,3) = max(max(rhox(indx(2:end-1,2:end-1,1))));
else
    for i = 1:cal.nc
        hist.cx(stp,1:3,i) = NaN;
        hist.cx_oxd(stp,1:3,i) = NaN;
    end
    hist.rhox(stp,1:3) = NaN;
end

indm = repmat(m>1e-6,1,1,cal.nc);
if any(indm(:)>0)
    for i = 1:cal.nc
        hist.cm(stp,1,i) = min(min(cm(indm(2:end-1,2:end-1,i))));
        hist.cm(stp,2,i) = sum(sum(cm(2:end-1,2:end-1,i).*m(2:end-1,2:end-1).*rho(2:end-1,2:end-1)))./sum(sum(m(2:end-1,2:end-1).*rho(2:end-1,2:end-1)));
        hist.cm(stp,3,i) = max(max(cm(indm(2:end-1,2:end-1,i))));

        hist.cm_oxd(stp,1,i) = min(min(cm_oxd(indm(2:end-1,2:end-1,i))));
        hist.cm_oxd(stp,2,i) = sum(sum(cm_oxd(2:end-1,2:end-1,i).*m(2:end-1,2:end-1).*rho(2:end-1,2:end-1)))./sum(sum(m(2:end-1,2:end-1).*rho(2:end-1,2:end-1)));
        hist.cm_oxd(stp,3,i) = max(max(cm_oxd(indm(2:end-1,2:end-1,i))));
    end
    hist.rhom(stp,1) = min(min(rhom(indm(2:end-1,2:end-1,1))));
    hist.rhom(stp,2) = sum(sum(rhom(2:end-1,2:end-1).*m(2:end-1,2:end-1)))./sum(sum(m(2:end-1,2:end-1)));
    hist.rhom(stp,3) = max(max(rhom(indm(2:end-1,2:end-1,1))));

    hist.etam(stp,1) = min(min(etam(indm(2:end-1,2:end-1,1))));
    hist.etam(stp,2) = sum(sum(etam(2:end-1,2:end-1).*m(2:end-1,2:end-1)))./sum(sum(m(2:end-1,2:end-1)));
    hist.etam(stp,3) = max(max(etam(indm(2:end-1,2:end-1,1))));
else
    for i = 1:cal.nc
        hist.cm(stp,1:3,i) = NaN;
        hist.cm_oxd(stp,1:3,i) = NaN;
    end
    hist.rhom(stp,1:3) = NaN;
    hist.etam(stp,1:3) = NaN;
end

hist.Gx(stp,1) = min( Gx,[],'all');
hist.Gx(stp,2) = mean(Gx   ,'all');
hist.Gx(stp,3) = max( Gx,[],'all');

hist.dV(stp,1) = min( VolSrc,[],'all');
hist.dV(stp,2) = mean(VolSrc   ,'all');
hist.dV(stp,3) = max( VolSrc,[],'all');

hist.rho(stp,1) = min( rho(2:end-1,2:end-1),[],'all');
hist.rho(stp,2) = mean(rho(2:end-1,2:end-1)   ,'all');
hist.rho(stp,3) = max( rho(2:end-1,2:end-1),[],'all');

hist.eta(stp,1) = min(    eta(2:end-1,2:end-1),[],'all');
hist.eta(stp,2) = geomean(eta(2:end-1,2:end-1)   ,'all');
hist.eta(stp,3) = max(    eta(2:end-1,2:end-1),[],'all');