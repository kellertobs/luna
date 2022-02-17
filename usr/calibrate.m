%% prep workspace
clear all; close all;
addpath('../src');

% prep plotting options
TX = {'Interpreter','latex'};
TL = {'TickLabelInterpreter','latex'};
FS = {'FontSize',16,20};
MS = {'MarkerSize',8};
LW = {'LineWidth',2};
CL = {'Color',[0.4660, 0.6740, 0.1880],[0.4940, 0.1840, 0.5560],[0, 0.4470, 0.7410],[0.8500, 0.3250, 0.0980],[0.8,0.8,0.8]};
TINY = 1e-16;

%% load experimental data from Schmidt & Kraettli (2020), Table 3
Load_SKTable3;
nphs = 6; olv=1; opx=2; cpx=3; plg=4; tis=5; mlt=6;
ncmp = 11; Si=1; Ti=2; Al=3; Cr=4; Fe=5; Mn=6; Mg=7; Ca=8; Na=9; K=10; P=11;

oxds(:,tis,:) = [];  fphs(:,tis) = []; % ignore Ti-sp phase

oxds(:,:,P ) = [];  % ignore P2O5
oxds(:,:,K ) = [];  % ignore K2O
oxds(:,:,Mn) = [];  % ignore MnO
oxds(:,:,Cr) = [];  % ignore Cr2O3
oxds(:,:,Ti) = [];  % ignore TiO2

oxds(:,opx,:) = (fphs(:,opx).*squeeze(oxds(:,opx,:)) + fphs(:,cpx).*squeeze(oxds(:,cpx,:)))./(fphs(:,opx) + fphs(:,cpx) + TINY);
fphs(:,opx  ) = fphs(:,opx) + fphs(:,cpx); 
oxds(:,cpx,:) = []; fphs(:,cpx) = [];

oxds = oxds./(sum(oxds,3)+TINY)*100;  % normalise remaining oxides to 100%

nphs = 4; olv=1; pxn=2; plg=3; mlt=4; blk=5;
ncmp = 6; Si=1; Al=2; Fe=3; Mg=4; Ca=5; Na=6;

%% load calibration
cal_magma6;

% plot Harker diagrams
figure(1); clf;
subplot(2,3,1);
plot(Temp(fphs(:,olv)>0),fphs(fphs(:,olv)>0,olv),'d',CL{[1,2]},MS{:},LW{1},1.5); axis xy tight; box on; hold on
plot(Temp(fphs(:,pxn)>0),fphs(fphs(:,pxn)>0,pxn),'v',CL{[1,3]},MS{:},LW{1},1.5); axis xy tight; box on; hold on
plot(Temp(fphs(:,plg)>0),fphs(fphs(:,plg)>0,plg),'s',CL{[1,4]},MS{:},LW{1},1.5); axis xy tight; box on; hold on
plot(Temp(fphs(:,mlt)>0),fphs(fphs(:,mlt)>0,mlt),'o',CL{[1,5]},MS{:},LW{1},1.5); axis xy tight; box on; hold on
set(gca,TL{:}); xlabel('$T [^\circ$C]',TX{:}); ylabel('$f_\mathrm{phs}$ [wt]',TX{:});

subplot(2,3,2);
plot(oxds(:,blk,Si),oxds(:,blk,Al),'k.',MS{1},12,LW{1},1.5); axis xy tight; box on; hold on
plot(oxds(fphs(:,olv)>0,olv,Si),oxds(fphs(:,olv)>0,olv,Al),'d',CL{[1,2]},MS{:},LW{1},1.5);
plot(oxds(fphs(:,pxn)>0,pxn,Si),oxds(fphs(:,pxn)>0,pxn,Al),'v',CL{[1,3]},MS{:},LW{1},1.5); 
plot(oxds(fphs(:,plg)>0,plg,Si),oxds(fphs(:,plg)>0,plg,Al),'s',CL{[1,4]},MS{:},LW{1},1.5); 
plot(oxds(fphs(:,mlt)>0,mlt,Si),oxds(fphs(:,mlt)>0,mlt,Al),'o',CL{[1,5]},MS{:},LW{1},1.5);
plot(cal.oxds(cal.fo ,Si),cal.oxds(cal.fo ,Al),'d',CL{[1,2]},MS{1},12,LW{:});
plot(cal.oxds(cal.fay,Si),cal.oxds(cal.fay,Al),'d',CL{[1,2]},MS{1},12,LW{:});
plot(cal.oxds(cal.opx,Si),cal.oxds(cal.opx,Al),'v',CL{[1,3]},MS{1},12,LW{:});
plot(cal.oxds(cal.cpx,Si),cal.oxds(cal.cpx,Al),'v',CL{[1,3]},MS{1},12,LW{:});
plot(cal.oxds(cal.an ,Si),cal.oxds(cal.an ,Al),'s',CL{[1,4]},MS{1},12,LW{:});
plot(cal.oxds(cal.ab ,Si),cal.oxds(cal.ab ,Al),'s',CL{[1,4]},MS{1},12,LW{:});
set(gca,TL{:}); xlabel('SiO$_2$ [wt\%]',TX{:}); ylabel('Al$_2$O$_3$ [wt\%]',TX{:});

subplot(2,3,3);
plot(oxds(:,blk,Si),oxds(:,blk,Fe),'k.',MS{1},12,LW{1},1.5); axis xy tight; box on; hold on
plot(oxds(fphs(:,olv)>0,olv,Si),oxds(fphs(:,olv)>0,olv,Fe),'d',CL{[1,2]},MS{:},LW{1},1.5);
plot(oxds(fphs(:,pxn)>0,pxn,Si),oxds(fphs(:,pxn)>0,pxn,Fe),'v',CL{[1,3]},MS{:},LW{1},1.5); 
plot(oxds(fphs(:,plg)>0,plg,Si),oxds(fphs(:,plg)>0,plg,Fe),'s',CL{[1,4]},MS{:},LW{1},1.5); 
plot(oxds(fphs(:,mlt)>0,mlt,Si),oxds(fphs(:,mlt)>0,mlt,Fe),'o',CL{[1,5]},MS{:},LW{1},1.5);
plot(cal.oxds(cal.fo ,Si),cal.oxds(cal.fo ,Fe),'d',CL{[1,2]},MS{1},12,LW{:});
plot(cal.oxds(cal.fay,Si),cal.oxds(cal.fay,Fe),'d',CL{[1,2]},MS{1},12,LW{:});
plot(cal.oxds(cal.opx,Si),cal.oxds(cal.opx,Fe),'v',CL{[1,3]},MS{1},12,LW{:});
plot(cal.oxds(cal.cpx,Si),cal.oxds(cal.cpx,Fe),'v',CL{[1,3]},MS{1},12,LW{:});
plot(cal.oxds(cal.an ,Si),cal.oxds(cal.an ,Fe),'s',CL{[1,4]},MS{1},12,LW{:});
plot(cal.oxds(cal.ab ,Si),cal.oxds(cal.ab ,Fe),'s',CL{[1,4]},MS{1},12,LW{:});
set(gca,TL{:}); xlabel('SiO$_2$ [wt\%]',TX{:}); ylabel('FeO [wt\%]',TX{:});

subplot(2,3,4);
plot(oxds(:,blk,Si),oxds(:,blk,Mg),'k.',MS{1},12,LW{1},1.5); axis xy tight; box on; hold on
plot(oxds(fphs(:,olv)>0,olv,Si),oxds(fphs(:,olv)>0,olv,Mg),'d',CL{[1,2]},MS{:},LW{1},1.5);
plot(oxds(fphs(:,pxn)>0,pxn,Si),oxds(fphs(:,pxn)>0,pxn,Mg),'v',CL{[1,3]},MS{:},LW{1},1.5); 
plot(oxds(fphs(:,plg)>0,plg,Si),oxds(fphs(:,plg)>0,plg,Mg),'s',CL{[1,4]},MS{:},LW{1},1.5); 
plot(oxds(fphs(:,mlt)>0,mlt,Si),oxds(fphs(:,mlt)>0,mlt,Mg),'o',CL{[1,5]},MS{:},LW{1},1.5);
plot(cal.oxds(cal.fo ,Si),cal.oxds(cal.fo ,Mg),'d',CL{[1,2]},MS{1},12,LW{:});
plot(cal.oxds(cal.fay,Si),cal.oxds(cal.fay,Mg),'d',CL{[1,2]},MS{1},12,LW{:});
plot(cal.oxds(cal.opx,Si),cal.oxds(cal.opx,Mg),'v',CL{[1,3]},MS{1},12,LW{:});
plot(cal.oxds(cal.cpx,Si),cal.oxds(cal.cpx,Mg),'v',CL{[1,3]},MS{1},12,LW{:});
plot(cal.oxds(cal.an ,Si),cal.oxds(cal.an ,Mg),'s',CL{[1,4]},MS{1},12,LW{:});
plot(cal.oxds(cal.ab ,Si),cal.oxds(cal.ab ,Mg),'s',CL{[1,4]},MS{1},12,LW{:});
set(gca,TL{:}); xlabel('SiO$_2$ [wt\%]',TX{:}); ylabel('MgO [wt\%]',TX{:});

subplot(2,3,5);
plot(oxds(:,blk,Si),oxds(:,blk,Ca),'k.',MS{1},12,LW{1},1.5); axis xy tight; box on; hold on
plot(oxds(fphs(:,olv)>0,olv,Si),oxds(fphs(:,olv)>0,olv,Ca),'d',CL{[1,2]},MS{:},LW{1},1.5);
plot(oxds(fphs(:,pxn)>0,pxn,Si),oxds(fphs(:,pxn)>0,pxn,Ca),'v',CL{[1,3]},MS{:},LW{1},1.5); 
plot(oxds(fphs(:,plg)>0,plg,Si),oxds(fphs(:,plg)>0,plg,Ca),'s',CL{[1,4]},MS{:},LW{1},1.5); 
plot(oxds(fphs(:,mlt)>0,mlt,Si),oxds(fphs(:,mlt)>0,mlt,Ca),'o',CL{[1,5]},MS{:},LW{1},1.5);
plot(cal.oxds(cal.fo ,Si),cal.oxds(cal.fo ,Ca),'d',CL{[1,2]},MS{1},12,LW{:});
plot(cal.oxds(cal.fay,Si),cal.oxds(cal.fay,Ca),'d',CL{[1,2]},MS{1},12,LW{:});
plot(cal.oxds(cal.opx,Si),cal.oxds(cal.opx,Ca),'v',CL{[1,3]},MS{1},12,LW{:});
plot(cal.oxds(cal.cpx,Si),cal.oxds(cal.cpx,Ca),'v',CL{[1,3]},MS{1},12,LW{:});
plot(cal.oxds(cal.an ,Si),cal.oxds(cal.an ,Ca),'s',CL{[1,4]},MS{1},12,LW{:});
plot(cal.oxds(cal.ab ,Si),cal.oxds(cal.ab ,Ca),'s',CL{[1,4]},MS{1},12,LW{:});
set(gca,TL{:}); xlabel('SiO$_2$ [wt\%]',TX{:}); ylabel('CaO [wt\%]',TX{:});

subplot(2,3,6);
plot(oxds(:,blk,Si),oxds(:,blk,Na),'k.',MS{1},12,LW{1},1.5); axis xy tight; box on; hold on
plot(oxds(fphs(:,olv)>0,olv,Si),oxds(fphs(:,olv)>0,olv,Na),'d',CL{[1,2]},MS{:},LW{1},1.5);
plot(oxds(fphs(:,pxn)>0,pxn,Si),oxds(fphs(:,pxn)>0,pxn,Na),'v',CL{[1,3]},MS{:},LW{1},1.5); 
plot(oxds(fphs(:,plg)>0,plg,Si),oxds(fphs(:,plg)>0,plg,Na),'s',CL{[1,4]},MS{:},LW{1},1.5); 
plot(oxds(fphs(:,mlt)>0,mlt,Si),oxds(fphs(:,mlt)>0,mlt,Na),'o',CL{[1,5]},MS{:},LW{1},1.5);
plot(cal.oxds(cal.fo ,Si),cal.oxds(cal.fo ,Na),'d',CL{[1,2]},MS{1},12,LW{:});
plot(cal.oxds(cal.fay,Si),cal.oxds(cal.fay,Na),'d',CL{[1,2]},MS{1},12,LW{:});
plot(cal.oxds(cal.opx,Si),cal.oxds(cal.opx,Na),'v',CL{[1,3]},MS{1},12,LW{:});
plot(cal.oxds(cal.cpx,Si),cal.oxds(cal.cpx,Na),'v',CL{[1,3]},MS{1},12,LW{:});
plot(cal.oxds(cal.an ,Si),cal.oxds(cal.an ,Na),'s',CL{[1,4]},MS{1},12,LW{:});
plot(cal.oxds(cal.ab ,Si),cal.oxds(cal.ab ,Na),'s',CL{[1,4]},MS{1},12,LW{:});
set(gca,TL{:}); xlabel('SiO$_2$ [wt\%]',TX{:}); ylabel('Na$_2$O [wt\%]',TX{:});


%% calibrate model against experiments

cal0  = cal;
msft  = 1e3;
bstft = msft;
tol   = 7.5e-2;

while msft > tol

red = 1;%bstft^0.1;

% draw randomised calibration coefficients
% set pure component melting points T_m^i at P=0
cal.T0(cal.fo )  =  cal0.T0(cal.fo ) + 10.*randn(1).*red;
cal.T0(cal.fay)  =  cal0.T0(cal.fay) + 10.*randn(1).*red;
cal.T0(cal.opx)  =  cal0.T0(cal.opx) + 10.*randn(1).*red;
cal.T0(cal.cpx)  =  cal0.T0(cal.cpx) + 10.*randn(1).*red;
cal.T0(cal.an )  =  cal0.T0(cal.an ) + 10.*randn(1).*red;
cal.T0(cal.ab )  =  cal0.T0(cal.ab ) + 10.*randn(1).*red;

% set first coeff. for P-dependence of T_m^i [GPa]
cal.A(cal.fo )   =  max(1,cal0.A(cal.fo ) + 0.2.*randn(1)).*red;
cal.A(cal.fay)   =  max(1,cal0.A(cal.fay) + 0.2.*randn(1)).*red;
cal.A(cal.opx)   =  max(1,cal0.A(cal.opx) + 0.2.*randn(1)).*red;
cal.A(cal.cpx)   =  max(1,cal0.A(cal.cpx) + 0.2.*randn(1)).*red;
cal.A(cal.an )   =  max(1,cal0.A(cal.an ) + 0.2.*randn(1)).*red;
cal.A(cal.ab )   =  max(1,cal0.A(cal.ab ) + 0.2.*randn(1)).*red;

% set second coeff. for P-dependence of T_m^i [1]
cal.B(cal.fo )   =  max(1,cal0.B(cal.fo ) + 0.02.*randn(1)).*red;
cal.B(cal.fay)   =  max(1,cal0.B(cal.fay) + 0.02.*randn(1)).*red;
cal.B(cal.opx)   =  max(1,cal0.B(cal.opx) + 0.02.*randn(1)).*red;
cal.B(cal.cpx)   =  max(1,cal0.B(cal.cpx) + 0.02.*randn(1)).*red;
cal.B(cal.an )   =  max(1,cal0.B(cal.an ) + 0.02.*randn(1)).*red;
cal.B(cal.ab )   =  max(1,cal0.B(cal.ab ) + 0.02.*randn(1)).*red;

% set coeff. for T-dependence of partition coefficients K^i [1/K]
cal.r(cal.fo )   =  max(5,cal0.r(cal.fo ) + 0.5.*randn(1)).*red;
cal.r(cal.fay)   =  max(5,cal0.r(cal.fay) + 0.5.*randn(1)).*red;
cal.r(cal.opx)   =  max(5,cal0.r(cal.opx) + 0.5.*randn(1)).*red;
cal.r(cal.cpx)   =  max(5,cal0.r(cal.cpx) + 0.5.*randn(1)).*red;
cal.r(cal.an )   =  max(5,cal0.r(cal.an ) + 0.5.*randn(1)).*red;
cal.r(cal.ab )   =  max(5,cal0.r(cal.ab ) + 0.5.*randn(1)).*red;

% select stages to include in parameters fitting
stages = [1,2,3,4,5,6,7,8,9,10];

% extract experimental data and compute model outcome for each stage
for stg = stages

    % get computed phase proportions, compositions
    c0       =  max(0,min(1,cal.oxds.'\squeeze(oxds(stg,blk,:))));  % find end-member proportions of first step
    var.c    =  c0.'./sum(c0);
    var.T    =  Temp(stg);
    var.P    =  Pres(stg);
    [phs,~]  =  meltmodel(var,cal,'E');

    % record model results for all stages
    mdl.f(stg   ) = phs.f;
    mdl.cs(stg,:) = phs.cs;
    mdl.cl(stg,:) = phs.cl;

    % exclude results where experiments do not have one of the phases
    hasolv(stg) = (fphs(stg,olv)>0);
    haspxn(stg) = (fphs(stg,pxn)>0);
    hasplg(stg) = (fphs(stg,plg)>0);

    % bring model result into form of experimental data
    mdl.fphs(stg,olv) = (phs.cs(cal.fo ) + phs.cs(cal.fay)) .* (1-phs.f) .* hasolv(stg);
    mdl.fphs(stg,pxn) = (phs.cs(cal.opx) + phs.cs(cal.cpx)) .* (1-phs.f) .* haspxn(stg);
    mdl.fphs(stg,plg) = (phs.cs(cal.an ) + phs.cs(cal.ab )) .* (1-phs.f) .* hasplg(stg);
    mdl.fphs(stg,mlt) = phs.f;

    mdl.oxds(stg,olv,:) = (phs.cs(cal.fo )*cal.oxds(cal.fo ,:) + phs.cs(cal.fay)*cal.oxds(cal.fay,:)) ./ (phs.cs(cal.fo )+phs.cs(cal.fay)+TINY) .* hasolv(stg);
    mdl.oxds(stg,pxn,:) = (phs.cs(cal.opx)*cal.oxds(cal.opx,:) + phs.cs(cal.cpx)*cal.oxds(cal.cpx,:)) ./ (phs.cs(cal.opx)+phs.cs(cal.cpx)+TINY) .* haspxn(stg);
    mdl.oxds(stg,plg,:) = (phs.cs(cal.an )*cal.oxds(cal.an ,:) + phs.cs(cal.ab )*cal.oxds(cal.ab ,:)) ./ (phs.cs(cal.an )+phs.cs(cal.ab )+TINY) .* hasplg(stg);
    mdl.oxds(stg,mlt,:) = phs.cl*cal.oxds;

    exp.Temp(stg,:  ) = Temp(stg);
    exp.fphs(stg,:  ) = fphs(stg,:);
    exp.oxds(stg,:,:) = oxds(stg,1:nphs,:);

end

% get misfit norm to characterise fit
msft = 1*norm(exp.fphs(:)-mdl.fphs(:),2)./sqrt(length(exp.fphs(:))) ...
     + 1*norm(exp.oxds(:)-mdl.oxds(:),2)./sqrt(length(exp.oxds(:)))/50;

% record if new bestfit found
if msft < 1.01*bstft
    cal0  = cal;
    bstft = msft;
    fprintf(1,'misfit = %1.4e\n',msft);

    % plot current fit
    figure(2); clf;
    subplot(2,3,1);
    plot(exp.Temp(hasolv),exp.fphs(hasolv,olv),'d',CL{[1,6]},MS{:},LW{1},1.5); axis xy tight; box on; hold on
    plot(exp.Temp(haspxn),exp.fphs(haspxn,pxn),'v',CL{[1,6]},MS{:},LW{1},1.5);
    plot(exp.Temp(hasplg),exp.fphs(hasplg,plg),'s',CL{[1,6]},MS{:},LW{1},1.5);
    plot(exp.Temp,exp.fphs(:,mlt),'o',CL{[1,6]},MS{:},LW{1},1.5);
    plot(exp.Temp(hasolv),mdl.fphs(hasolv,olv),'d',CL{[1,2]},MS{:},LW{1},1.5);
    plot(exp.Temp(haspxn),mdl.fphs(haspxn,pxn),'v',CL{[1,3]},MS{:},LW{1},1.5);
    plot(exp.Temp(hasplg),mdl.fphs(hasplg,plg),'s',CL{[1,4]},MS{:},LW{1},1.5);
    plot(exp.Temp,mdl.fphs(:,mlt),'o',CL{[1,5]},MS{:},LW{1},1.5);
    set(gca,TL{:}); xlabel('$T [^\circ$C]',TX{:}); ylabel('$f_\mathrm{phs}$ [wt]',TX{:});

    subplot(2,3,2);
    plot(exp.oxds(hasolv,olv,Si),exp.oxds(hasolv,olv,Al),'d',CL{[1,6]},MS{:},LW{1},1.5); axis xy tight; box on; hold on
    plot(exp.oxds(haspxn,pxn,Si),exp.oxds(haspxn,pxn,Al),'v',CL{[1,6]},MS{:},LW{1},1.5);
    plot(exp.oxds(hasplg,plg,Si),exp.oxds(hasplg,plg,Al),'s',CL{[1,6]},MS{:},LW{1},1.5);
    plot(exp.oxds(:,mlt,Si),exp.oxds(:,mlt,Al),'o',CL{[1,6]},MS{:},LW{1},1.5);
    plot(mdl.oxds(hasolv,olv,Si),mdl.oxds(hasolv,olv,Al),'d',CL{[1,2]},MS{:},LW{1},1.5);
    plot(mdl.oxds(haspxn,pxn,Si),mdl.oxds(haspxn,pxn,Al),'v',CL{[1,3]},MS{:},LW{1},1.5);
    plot(mdl.oxds(hasplg,plg,Si),mdl.oxds(hasplg,plg,Al),'s',CL{[1,4]},MS{:},LW{1},1.5);
    plot(mdl.oxds(:,mlt,Si),mdl.oxds(:,mlt,Al),'o',CL{[1,5]},MS{:},LW{1},1.5);
    set(gca,TL{:}); xlabel('SiO$_2$ [wt\%]',TX{:}); ylabel('Al$_2$O$_3$ [wt\%]',TX{:});

    subplot(2,3,3);
    plot(exp.oxds(hasolv,olv,Si),exp.oxds(hasolv,olv,Fe),'d',CL{[1,6]},MS{:},LW{1},1.5); axis xy tight; box on; hold on
    plot(exp.oxds(haspxn,pxn,Si),exp.oxds(haspxn,pxn,Fe),'v',CL{[1,6]},MS{:},LW{1},1.5);
    plot(exp.oxds(hasplg,plg,Si),exp.oxds(hasplg,plg,Fe),'s',CL{[1,6]},MS{:},LW{1},1.5);
    plot(exp.oxds(:,mlt,Si),exp.oxds(:,mlt,Fe),'o',CL{[1,6]},MS{:},LW{1},1.5);
    plot(mdl.oxds(hasolv,olv,Si),mdl.oxds(hasolv,olv,Fe),'d',CL{[1,2]},MS{:},LW{1},1.5);
    plot(mdl.oxds(haspxn,pxn,Si),mdl.oxds(haspxn,pxn,Fe),'v',CL{[1,3]},MS{:},LW{1},1.5);
    plot(mdl.oxds(hasplg,plg,Si),mdl.oxds(hasplg,plg,Fe),'s',CL{[1,4]},MS{:},LW{1},1.5);
    plot(mdl.oxds(:,mlt,Si),mdl.oxds(:,mlt,Fe),'o',CL{[1,5]},MS{:},LW{1},1.5);
    set(gca,TL{:}); xlabel('SiO$_2$ [wt\%]',TX{:}); ylabel('FeO [wt\%]',TX{:});

    subplot(2,3,4);
    plot(exp.oxds(hasolv,olv,Si),exp.oxds(hasolv,olv,Mg),'d',CL{[1,6]},MS{:},LW{1},1.5); axis xy tight; box on; hold on
    plot(exp.oxds(haspxn,pxn,Si),exp.oxds(haspxn,pxn,Mg),'v',CL{[1,6]},MS{:},LW{1},1.5);
    plot(exp.oxds(hasplg,plg,Si),exp.oxds(hasplg,plg,Mg),'s',CL{[1,6]},MS{:},LW{1},1.5);
    plot(exp.oxds(:,mlt,Si),exp.oxds(:,mlt,Mg),'o',CL{[1,6]},MS{:},LW{1},1.5);
    plot(mdl.oxds(hasolv,olv,Si),mdl.oxds(hasolv,olv,Mg),'d',CL{[1,2]},MS{:},LW{1},1.5);
    plot(mdl.oxds(haspxn,pxn,Si),mdl.oxds(haspxn,pxn,Mg),'v',CL{[1,3]},MS{:},LW{1},1.5);
    plot(mdl.oxds(hasplg,plg,Si),mdl.oxds(hasplg,plg,Mg),'s',CL{[1,4]},MS{:},LW{1},1.5);
    plot(mdl.oxds(:,mlt,Si),mdl.oxds(:,mlt,Mg),'o',CL{[1,5]},MS{:},LW{1},1.5);
    set(gca,TL{:}); xlabel('SiO$_2$ [wt\%]',TX{:}); ylabel('MgO [wt\%]',TX{:});

    subplot(2,3,5);
    plot(exp.oxds(hasolv,olv,Si),exp.oxds(hasolv,olv,Ca),'d',CL{[1,6]},MS{:},LW{1},1.5); axis xy tight; box on; hold on
    plot(exp.oxds(haspxn,pxn,Si),exp.oxds(haspxn,pxn,Ca),'v',CL{[1,6]},MS{:},LW{1},1.5);
    plot(exp.oxds(hasplg,plg,Si),exp.oxds(hasplg,plg,Ca),'s',CL{[1,6]},MS{:},LW{1},1.5);
    plot(exp.oxds(:,mlt,Si),exp.oxds(:,mlt,Ca),'o',CL{[1,6]},MS{:},LW{1},1.5);
    plot(mdl.oxds(hasolv,olv,Si),mdl.oxds(hasolv,olv,Ca),'d',CL{[1,2]},MS{:},LW{1},1.5);
    plot(mdl.oxds(haspxn,pxn,Si),mdl.oxds(haspxn,pxn,Ca),'v',CL{[1,3]},MS{:},LW{1},1.5);
    plot(mdl.oxds(hasplg,plg,Si),mdl.oxds(hasplg,plg,Ca),'s',CL{[1,4]},MS{:},LW{1},1.5);
    plot(mdl.oxds(:,mlt,Si),mdl.oxds(:,mlt,Ca),'o',CL{[1,5]},MS{:},LW{1},1.5);
    set(gca,TL{:}); xlabel('SiO$_2$ [wt\%]',TX{:}); ylabel('CaO [wt\%]',TX{:});

    subplot(2,3,6);
    plot(exp.oxds(hasolv,olv,Si),exp.oxds(hasolv,olv,Na),'d',CL{[1,6]},MS{:},LW{1},1.5); axis xy tight; box on; hold on
    plot(exp.oxds(haspxn,pxn,Si),exp.oxds(haspxn,pxn,Na),'v',CL{[1,6]},MS{:},LW{1},1.5);
    plot(exp.oxds(hasplg,plg,Si),exp.oxds(hasplg,plg,Na),'s',CL{[1,6]},MS{:},LW{1},1.5);
    plot(exp.oxds(:,mlt,Si),exp.oxds(:,mlt,Na),'o',CL{[1,6]},MS{:},LW{1},1.5);
    plot(mdl.oxds(hasolv,olv,Si),mdl.oxds(hasolv,olv,Na),'d',CL{[1,2]},MS{:},LW{1},1.5);
    plot(mdl.oxds(haspxn,pxn,Si),mdl.oxds(haspxn,pxn,Na),'v',CL{[1,3]},MS{:},LW{1},1.5);
    plot(mdl.oxds(hasplg,plg,Si),mdl.oxds(hasplg,plg,Na),'s',CL{[1,4]},MS{:},LW{1},1.5);
    plot(mdl.oxds(:,mlt,Si),mdl.oxds(:,mlt,Na),'o',CL{[1,5]},MS{:},LW{1},1.5);
    set(gca,TL{:}); xlabel('SiO$_2$ [wt\%]',TX{:}); ylabel('Na$_2$O [wt\%]',TX{:});

end

end

%% plot best fit calibration
% set bulk composition as ['fo','fa','opx','px2','cpx,'an','ab']
c0 = max(0,min(1,cal.oxds.'\squeeze(oxds(1,blk,:))));  % find end-member proportions of first step
c0 = c0.'./sum(c0);

% plot partition coefficients (T-dependent)
var.c = ones(100,1)*c0;
var.T = linspace(1000,1800,100).';
var.P = linspace(1.5,1.5,100).';

[~,cal]  =  meltmodel(var,cal,'K');

figure(3); clf;
semilogy(var.T,cal.K,LW{:}); axis xy tight; box on;
legend(cal.CompStr{:},FS{[1,2]},TX{:})
set(gca,TL{:},FS{[1,2]})
xlabel('$T \ [^\circ$C]',TX{:},FS{[1,3]})
ylabel('log$_{10} \ K$',TX{:},FS{[1,3]})


% plot partition coefficients (P-dependent)
var.c = ones(100,1)*c0;
var.T = linspace(1300,1300,100).';
var.P = linspace(0,3,100).';

[~,cal]  =  meltmodel(var,cal,'K');

figure(4); clf;
semilogx(cal.K,var.P,LW{:}); axis ij tight; box on;
legend(cal.CompStr{:},FS{[1,2]},TX{:})
set(gca,TL{:},FS{[1,2]})
xlabel('log$_{10} \ K$',TX{:},FS{[1,3]})
ylabel('$P$ [GPa]',TX{:},FS{[1,3]})

% plot melting points, solidus, liquidus (P-dependent)
var.c = ones(100,1)*c0;
var.T = linspace(1300,1300,100).';
var.P = linspace(0,3,100).';

[~,cal]  =  meltmodel(var,cal,'T');

figure(5); clf;
plot(cal.Tm,var.P,LW{:}); axis ij tight; box on; hold on
plot(cal.Tsol,var.P,'k--',LW{1},3); axis ij tight; box on;
plot(cal.Tliq,var.P,'k-.',LW{1},3); axis ij tight; box on;
legend(cal.CompStr{:},'$T_\mathrm{sol}$','$T_\mathrm{liq}$',FS{[1,2]},TX{:})
set(gca,TL{:},FS{[1,2]})
xlabel('$T \ [^\circ$C]',TX{:},FS{[1,3]})
ylabel('$P$ [GPa]',TX{:},FS{[1,3]})

% plot melt fraction, phase compositions (T-dependent)
var.c = ones(1000,1)*c0;
var.T = linspace(1300,1300,1000).';
var.P = linspace(1.5,1.5,1000).';

[~,cal]  =  meltmodel(var,cal,'T');
var.T = linspace(mean(cal.Tsol),mean(cal.Tliq),1000).';

[phs,~]  =  meltmodel(var,cal,'E');

figure(6); clf;
plot(var.T,1-phs.f,LW{:}); axis xy tight; box on; hold on;
plot(var.T,phs.f,LW{:}); axis xy tight; box on;
legend('xtals','melt',FS{[1,2]},TX{:},'location','east')
set(gca,TL{:},FS{[1,2]})
xlabel('$T \ [^\circ$C]',TX{:},FS{[1,3]})
ylabel('$f$ [wt]',TX{:},FS{[1,3]})

figure(7); clf;
plot(var.T,phs.cs,LW{:}); axis xy tight; box on;
legend(cal.CompStr{:},FS{[1,2]},TX{:})
set(gca,TL{:},FS{[1,2]})
xlabel('$T \ [^\circ$C]',TX{:},FS{[1,3]})
ylabel('$c^s$ [wt]',TX{:},FS{[1,3]})

figure(8); clf;
plot(var.T,phs.cl,LW{:}); axis xy tight; box on;
legend(cal.CompStr{:},FS{[1,2]},TX{:})
set(gca,TL{:},FS{[1,2]})
xlabel('$T \ [^\circ$C]',TX{:},FS{[1,3]})
ylabel('$c^\ell$ [wt]',TX{:},FS{[1,3]})

figure(9); clf;
plot(var.T,phs.cl*cal.oxds,LW{:}); axis xy tight; box on;
legend('SiO$_2$','Al$_2$O$_3$','FeO','MgO','CaO','Na$_2$O',FS{[1,2]},TX{:})
set(gca,TL{:},FS{[1,2]})
xlabel('$T \ [^\circ$C]',TX{:},FS{[1,3]})
ylabel('melt composition [wt\%]',TX{:},FS{[1,3]})
