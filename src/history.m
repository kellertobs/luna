% record run history

dsumMdto = dsumMdt;
dsumHdto = dsumHdt;
dsumCdto = dsumCdt;

stp = max(1,step);

% record model time
hist.time(stp) = time;

% record total mass, heat, component mass in model (assume hy = 1, unit length in third dimension)
hist.sumM(stp) = sum(sum(rho(2:end-1,2:end-1)*h*h*1));  % [kg]
hist.sumH(stp) = sum(sum(  H(2:end-1,2:end-1)*h*h*1));  % [J]
for i = 1:cal.nc; hist.sumC(stp,i) = sum(sum(  squeeze(C(i,2:end-1,2:end-1))*h*h*1));  end% [kg]

% record expected rates of change by volume change and imposed boundaries layers
dsumMdt = sum(rho(2,2:end-1).*W(1,2:end-1)*h*1) - sum(rho(end-1,2:end-1).*W(end,2:end-1)*h*1) ...
        + sum(rho(2:end-1,2).*U(2:end-1,1)*h*1) - sum(rho(2:end-1,end-1).*U(2:end-1,end)*h*1);  % [kg/s]
dsumHdt = sum(sum(bndH(2:end-1,2:end-1)*h*h*1)) ...
        + sum(  H(2,2:end-1).*W(1,2:end-1)*h*1) - sum(  H(end-1,2:end-1).*W(end,2:end-1)*h*1) ...
        + sum(  H(2:end-1,2).*U(2:end-1,1)*h*1) - sum(  H(2:end-1,end-1).*U(2:end-1,end)*h*1);  % [J /s]
for i = 1:cal.nc
    dsumCdt(i) = sum(  squeeze(C(i,2,2:end-1)).'.*W(1,2:end-1)*h*1) - sum(  squeeze(C(i,end-1,2:end-1)).'.*W(end,2:end-1)*h*1) ...
               + sum(  squeeze(C(i,2:end-1,2)).'.*U(2:end-1,1)*h*1) - sum(  squeeze(C(i,2:end-1,end-1)).'.*U(2:end-1,end)*h*1);  % [kg/s]
end

% if step>1; hist.DM(stp) = hist.DM(stp-1) + (THETA*dsumMdt + (1-THETA)*dsumMdto).*dt; else; hist.DM(stp) = 0; end  % [kg]
if step>1; hist.DM(stp) = hist.DM(stp-1) +        dsumMdt                      .*dt; else; hist.DM(stp) = 0; end  % [kg]
if step>1; hist.DH(stp) = hist.DH(stp-1) + (THETA*dsumHdt + (1-THETA)*dsumHdto).*dt; else; hist.DH(stp) = 0; end  % [J ]
if step>1; hist.DC(stp,:) = hist.DC(stp-1,:) + (THETA*dsumCdt + (1-THETA)*dsumCdto).*dt; else; hist.DC(stp,:) = zeros(1,cal.nc); end  % [kg]

% record conservation error of mass M, heat H, major component C, volatile component V
hist.EM(stp) = (hist.sumM(stp) - hist.DM(stp))./hist.sumM(1) - 1;  % [kg/kg]
hist.EH(stp) = (hist.sumH(stp) - hist.DH(stp))./hist.sumH(1) - 1;  % [J /J ]
hist.EC(stp,:) = (hist.sumC(stp,:) - hist.DC(stp,:))./hist.sumC(1,:) - 1;  % [kg/kg]

% record variable and coefficient diagnostics
hist.W(stp,1) = min(min(-W(:,2:end-1)));
hist.W(stp,2) = mean(mean(abs(W(:,2:end-1))));
hist.W(stp,3) = max(max(-W(:,2:end-1)));

hist.U(stp,1) = min(min(U(2:end-1,:)));
hist.U(stp,2) = mean(mean(abs(U(2:end-1,:))));
hist.U(stp,3) = max(max(U(2:end-1,:)));

hist.P(stp,1) = min(min(P(2:end-1,2:end-1)));
hist.P(stp,2) = mean(mean(abs(P(2:end-1,2:end-1))));
hist.P(stp,3) = max(max(P(2:end-1,2:end-1)));

hist.x(stp,1) = min(min(x(2:end-1,2:end-1)));
hist.x(stp,2) = mean(mean(x(2:end-1,2:end-1)));
hist.x(stp,3) = max(max(x(2:end-1,2:end-1)));

hist.m(stp,1) = min(min(m(2:end-1,2:end-1)));
hist.m(stp,2) = mean(mean(m(2:end-1,2:end-1)));
hist.m(stp,3) = max(max(m(2:end-1,2:end-1)));

hist.chi(stp,1) = min(min(chi(2:end-1,2:end-1)));
hist.chi(stp,2) = mean(mean(chi(2:end-1,2:end-1)));
hist.chi(stp,3) = max(max(chi(2:end-1,2:end-1)));

hist.mu(stp,1) = min(min(mu(2:end-1,2:end-1)));
hist.mu(stp,2) = mean(mean(mu(2:end-1,2:end-1)));
hist.mu(stp,3) = max(max(mu(2:end-1,2:end-1)));

hist.T(stp,1) = min(min(T(2:end-1,2:end-1)));
hist.T(stp,2) = mean(mean(T(2:end-1,2:end-1)));
hist.T(stp,3) = max(max(T(2:end-1,2:end-1)));

for i = 1:cal.nc
    hist.c(i,stp,1) = min(min(c(i,2:end-1,2:end-1)));
    hist.c(i,stp,2) = mean(mean(c(i,2:end-1,2:end-1)));
    hist.c(i,stp,3) = max(max(c(i,2:end-1,2:end-1)));
    hist.oxd(i,stp,1) = min(min(oxd(i,2:end-1,2:end-1)));
    hist.oxd(i,stp,2) = mean(mean(oxd(i,2:end-1,2:end-1)));
    hist.oxd(i,stp,3) = max(max(oxd(i,2:end-1,2:end-1)));
end

indx = x>1e-6;
if any(indx(:)>0)
    for i = 1:cal.nc
        hist.cx(i,stp,1) = min(min(cx(i,indx(2:end-1,2:end-1))));
        hist.cx(i,stp,2) = sum(sum(squeeze(cx(i,2:end-1,2:end-1)).*x(2:end-1,2:end-1).*rho(2:end-1,2:end-1)))./sum(sum(x(2:end-1,2:end-1).*rho(2:end-1,2:end-1)));
        hist.cx(i,stp,3) = max(max(cx(i,indx(2:end-1,2:end-1))));
        
        hist.oxdx(i,stp,1) = min(min(oxdx(i,indx(2:end-1,2:end-1))));
        hist.oxdx(i,stp,2) = sum(sum(squeeze(oxdx(i,2:end-1,2:end-1)).*x(2:end-1,2:end-1).*rho(2:end-1,2:end-1)))./sum(sum(x(2:end-1,2:end-1).*rho(2:end-1,2:end-1)));
        hist.oxdx(i,stp,3) = max(max(oxdx(i,indx(2:end-1,2:end-1))));
    end
    
    hist.rhox(stp,1) = min(min(rhox(indx(2:end-1,2:end-1))));
    hist.rhox(stp,2) = sum(sum(rhox(2:end-1,2:end-1).*chi(2:end-1,2:end-1)))./sum(sum(chi(2:end-1,2:end-1)));
    hist.rhox(stp,3) = max(max(rhox(indx(2:end-1,2:end-1))));
else
    hist.cx(1:cal.nc,stp,1:3) = NaN;
    hist.oxdx(1:cal.nc,stp,1:3) = NaN;
    hist.rhox(stp,1:3) = NaN;
end

indm = m>1e-6;
if any(indm(:)>0)
    for i = 1:cal.nc
        hist.cm(i,stp,1) = min(min(cm(i,indm(2:end-1,2:end-1))));
        hist.cm(i,stp,2) = sum(sum(squeeze(cm(i,2:end-1,2:end-1)).*m(2:end-1,2:end-1).*rho(2:end-1,2:end-1)))./sum(sum(m(2:end-1,2:end-1).*rho(2:end-1,2:end-1)));
        hist.cm(i,stp,3) = max(max(cm(i,indm(2:end-1,2:end-1))));
        
        hist.oxdm(i,stp,1) = min(min(oxdm(i,indm(2:end-1,2:end-1))));
        hist.oxdm(i,stp,2) = sum(sum(squeeze(oxdm(i,2:end-1,2:end-1)).*m(2:end-1,2:end-1).*rho(2:end-1,2:end-1)))./sum(sum(m(2:end-1,2:end-1).*rho(2:end-1,2:end-1)));
        hist.oxdm(i,stp,3) = max(max(oxdm(i,indm(2:end-1,2:end-1))));
    end
    
    hist.rhom(stp,1) = min(min(rhom(indm(2:end-1,2:end-1))));
    hist.rhom(stp,2) = sum(sum(rhom(2:end-1,2:end-1).*mu(2:end-1,2:end-1)))./sum(sum(mu(2:end-1,2:end-1)));
    hist.rhom(stp,3) = max(max(rhom(indm(2:end-1,2:end-1))));
    
    hist.etam(stp,1) = min(min(etam(indm(2:end-1,2:end-1))));
    hist.etam(stp,2) = geomean(geomean(etam(indm(2:end-1,2:end-1))));
    hist.etam(stp,3) = max(max(etam(indm(2:end-1,2:end-1))));
else
    hist.cm(1:cal.nc,stp,1:3) = NaN;
    hist.oxdm(1:cal.nc,stp,1:3) = NaN;
    hist.rhom(stp,1:3) = NaN;
    hist.etam(stp,1:3) = NaN;
end

hist.Gx(stp,1) = min(min(Gx(2:end-1,2:end-1)));
hist.Gx(stp,2) = mean(mean(Gx(2:end-1,2:end-1)));
hist.Gx(stp,3) = max(max(Gx(2:end-1,2:end-1)));

hist.dV(stp,1) = min(min(VolSrc(2:end-1,2:end-1)));
hist.dV(stp,2) = mean(mean(VolSrc(2:end-1,2:end-1)));
hist.dV(stp,3) = max(max(VolSrc(2:end-1,2:end-1)));

hist.rho(stp,1) = min(min(rho(2:end-1,2:end-1)));
hist.rho(stp,2) = mean(mean(rho(2:end-1,2:end-1)));
hist.rho(stp,3) = max(max(rho(2:end-1,2:end-1)));

hist.eta(stp,1) = min(min(eta(2:end-1,2:end-1)));
hist.eta(stp,2) = geomean(geomean(eta(2:end-1,2:end-1)));
hist.eta(stp,3) = max(max(eta(2:end-1,2:end-1)));

hist.wx(stp,1) = min(min(-(chi(1:end-1,2:end-1)+chi(2:end,2:end-1))/2.*wx(:,2:end-1)));
hist.wx(stp,2) = mean(mean(abs((chi(1:end-1,2:end-1)+chi(2:end,2:end-1))/2.*wx(:,2:end-1))));
hist.wx(stp,3) = max(max(-(chi(1:end-1,2:end-1)+chi(2:end,2:end-1))/2.*wx(:,2:end-1)));

hist.it(stp,1) = min(min(it(2:end-1,2:end-1)));
hist.it(stp,2) = mean(mean(it(2:end-1,2:end-1)));
hist.it(stp,3) = max(max(it(2:end-1,2:end-1)));

hist.ct(stp,1) = min(min(ct(2:end-1,2:end-1)));
hist.ct(stp,2) = mean(mean(ct(2:end-1,2:end-1)));
hist.ct(stp,3) = max(max(ct(2:end-1,2:end-1)));

hist.si(stp,1) = min(min(si(2:end-1,2:end-1)));
hist.si(stp,2) = mean(mean(si(2:end-1,2:end-1)));
hist.si(stp,3) = max(max(si(2:end-1,2:end-1)));

hist.rip(stp,1) = min(min(rip(2:end-1,2:end-1)));
hist.rip(stp,2) = mean(mean(rip(2:end-1,2:end-1)));
hist.rip(stp,3) = max(max(rip(2:end-1,2:end-1)));

hist.rid(stp,1) = min(min(rid(2:end-1,2:end-1)));
hist.rid(stp,2) = mean(mean(rid(2:end-1,2:end-1)));
hist.rid(stp,3) = max(max(rid(2:end-1,2:end-1)));

% % fraction, composition, and temperature of eruptible magma suspension (mu>0.55)
% indmagma = max(0,min(1,(1+erf((mu-0.55)./0.05))/2));
% hist.Fmagma(stp) = sum(sum(rho(2:end-1,2:end-1).*indmagma(2:end-1,2:end-1).*h^2))./sum(sum(rho(2:end-1,2:end-1).*h^2));
% hist.Cmagma(stp) = sum(sum(rho(2:end-1,2:end-1).*indmagma(2:end-1,2:end-1).*c(2:end-1,2:end-1).*h^2))./sum(sum(rho(2:end-1,2:end-1).*indmagma(2:end-1,2:end-1).*h^2));
% hist.Tmagma(stp) = sum(sum(rho(2:end-1,2:end-1).*indmagma(2:end-1,2:end-1).*T(2:end-1,2:end-1).*h^2))./sum(sum(rho(2:end-1,2:end-1).*indmagma(2:end-1,2:end-1).*h^2));
% 
% % fraction, composition, and temperature of plutonic rock (mu<0.15)
% indpluton = max(0,min(1,(1+erf((chi-0.85)./0.05))/2));
% hist.Fpluton(stp) = sum(sum(rho(2:end-1,2:end-1).*indpluton(2:end-1,2:end-1).*h^2))./sum(sum(rho(2:end-1,2:end-1).*h^2));
% hist.Cpluton(stp) = sum(sum(rho(2:end-1,2:end-1).*indpluton(2:end-1,2:end-1).*c(2:end-1,2:end-1).*h^2))./sum(sum(rho(2:end-1,2:end-1).*indpluton(2:end-1,2:end-1).*h^2));
% hist.Tpluton(stp) = sum(sum(rho(2:end-1,2:end-1).*indpluton(2:end-1,2:end-1).*T(2:end-1,2:end-1).*h^2))./sum(sum(rho(2:end-1,2:end-1).*indpluton(2:end-1,2:end-1).*h^2));
% 
% % fraction, composition, and temperature of magma mush (0.15<mu<0.55)
% indmush = max(0,min(1,1-indmagma-indpluton));
% hist.Fmush(stp) = sum(sum(rho(2:end-1,2:end-1).*indmush(2:end-1,2:end-1).*h^2))./sum(sum(rho(2:end-1,2:end-1).*h^2));
% hist.Cmush(stp) = sum(sum(rho(2:end-1,2:end-1).*indmush(2:end-1,2:end-1).*c(2:end-1,2:end-1).*h^2))./sum(sum(rho(2:end-1,2:end-1).*indmush(2:end-1,2:end-1).*h^2));
% hist.Tmush(stp) = sum(sum(rho(2:end-1,2:end-1).*indmush(2:end-1,2:end-1).*T(2:end-1,2:end-1).*h^2))./sum(sum(rho(2:end-1,2:end-1).*indmush(2:end-1,2:end-1).*h^2));
% 
% % fraction, crystallinity, and temperature of felsic materials (c > (perCm_cphs1)/2)
% indfelsic = max(0,min(1,(1+erf((c-(perCm+cphs1)/2)./0.005))/2));
% hist.Ffelsic(stp) = sum(sum(rho(2:end-1,2:end-1).*indfelsic(2:end-1,2:end-1).*h^2))./sum(sum(rho(2:end-1,2:end-1).*h^2));
% hist.Xfelsic(stp) = sum(sum(rho(2:end-1,2:end-1).*indfelsic(2:end-1,2:end-1).*x(2:end-1,2:end-1).*h^2))./sum(sum(rho(2:end-1,2:end-1).*indfelsic(2:end-1,2:end-1).*h^2));
% hist.Tfelsic(stp) = sum(sum(rho(2:end-1,2:end-1).*indfelsic(2:end-1,2:end-1).*T(2:end-1,2:end-1).*h^2))./sum(sum(rho(2:end-1,2:end-1).*indfelsic(2:end-1,2:end-1).*h^2));
% 
% % fraction, crystallinity, and temperature of intermediate materials (perCm < c < (perCm_cphs1)/2)
% indinterm = max(0,min(1,(1+erf((c-perCm)./0.005))/2 .* (1-indfelsic)));
% hist.Finterm(stp) = sum(sum(rho(2:end-1,2:end-1).*indinterm(2:end-1,2:end-1).*h^2))./sum(sum(rho(2:end-1,2:end-1).*h^2));
% hist.Xinterm(stp) = sum(sum(rho(2:end-1,2:end-1).*indinterm(2:end-1,2:end-1).*x(2:end-1,2:end-1).*h^2))./sum(sum(rho(2:end-1,2:end-1).*indinterm(2:end-1,2:end-1).*h^2));
% hist.Tinterm(stp) = sum(sum(rho(2:end-1,2:end-1).*indinterm(2:end-1,2:end-1).*T(2:end-1,2:end-1).*h^2))./sum(sum(rho(2:end-1,2:end-1).*indinterm(2:end-1,2:end-1).*h^2));
% 
% % fraction, crystallinity, and temperature of mafic materials (perCx < c < perCm)
% indmafic = max(0,min(1,(1+erf((c-perCx)./0.005))/2 .* (1-indinterm-indfelsic)));
% hist.Fmafic(stp) = sum(sum(rho(2:end-1,2:end-1).*indmafic(2:end-1,2:end-1).*h^2))./sum(sum(rho(2:end-1,2:end-1).*h^2));
% hist.Xmafic(stp) = sum(sum(rho(2:end-1,2:end-1).*indmafic(2:end-1,2:end-1).*x(2:end-1,2:end-1).*h^2))./sum(sum(rho(2:end-1,2:end-1).*indmafic(2:end-1,2:end-1).*h^2));
% hist.Tmafic(stp) = sum(sum(rho(2:end-1,2:end-1).*indmafic(2:end-1,2:end-1).*T(2:end-1,2:end-1).*h^2))./sum(sum(rho(2:end-1,2:end-1).*indmafic(2:end-1,2:end-1).*h^2));
% 
% % fraction, crystallinity, and temperature of ultramafic materials (c < perCx)
% indultram = max(0,min(1,1-indmafic-indinterm-indfelsic));
% hist.Fultram(stp) = sum(sum(rho(2:end-1,2:end-1).*indultram(2:end-1,2:end-1).*h^2))./sum(sum(rho(2:end-1,2:end-1).*h^2));
% hist.Xultram(stp) = sum(sum(rho(2:end-1,2:end-1).*indultram(2:end-1,2:end-1).*x(2:end-1,2:end-1).*h^2))./sum(sum(rho(2:end-1,2:end-1).*indultram(2:end-1,2:end-1).*h^2));
% hist.Tultram(stp) = sum(sum(rho(2:end-1,2:end-1).*indultram(2:end-1,2:end-1).*T(2:end-1,2:end-1).*h^2))./sum(sum(rho(2:end-1,2:end-1).*indultram(2:end-1,2:end-1).*h^2));
% 
% % differentiation index
% nobnd = bndshape<1e-2;
% if any(nobnd(:))
%     hist.Rdiff(stp) = (max(max(c(nobnd(2:end-1,2:end-1))))-min(min(c(nobnd(2:end-1,2:end-1)))))./(cphs1-cphs0);
% end