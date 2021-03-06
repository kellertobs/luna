%% prep workspace
clear all; close all;
addpath('../src');
addpath('../cal');

% prep plotting options
TX = {'Interpreter','latex'};
TL = {'TickLabelInterpreter','latex'};
FS = {'FontSize',16,20};
MS = {'MarkerSize',8};
LW = {'LineWidth',2};
CL = {'Color',[0.4660, 0.6740, 0.1880],[0.4940, 0.1840, 0.5560],[0, 0.4470, 0.7410],[0.8500, 0.3250, 0.0980],[0.8,0.8,0.8]};
TINY = 1e-16;

%% catmip parameters - adjust this

% catmip options
Niter   = 500;            % number of iterations per catmip temperature
Nstep   = 100;            % number of steps in MCMC in catmip
pllopt  = 8;              % whether to run in parallel, number of workers

rng(10); % for reproducibility during testing, comment this for actual run

%% load calibration
cal_catmip

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

% assign errors for phases and oxides - please review
factor             = 2; % to adjust sigmas, maybe too small
sig_fphs           = factor*0.020.*ones(size(fphs));
sig_oxds           = zeros(size(oxds));
sig_oxds(oxds>=20) = factor*max(0.02*oxds(oxds>=20), 1.00);
sig_oxds(oxds< 20) = factor*max(0.10*oxds(oxds< 20), 0.20);
sig_oxds(:,pxn,:)  = 2*sig_oxds(:,pxn,:);

%% prepare information for the parameter fitting

% select stages to include in parameters fitting
stages = [1,2,3,4,5,6,7,8,9,10];

for stg = stages
    
    % exclude results where experiments do not have one of the phases
    hasolv(stg) = (fphs(stg,olv)>0);
    haspxn(stg) = (fphs(stg,pxn)>0);
    hasplg(stg) = (fphs(stg,plg)>0);
    
    exp.Temp(stg,:  ) = Temp(stg);
    exp.fphs(stg,:  ) = fphs(stg,:);
    exp.oxds(stg,:,:) = oxds(stg,1:nphs,:);
    
    sig.fphs(stg,:  ) = sig_fphs(stg,:);
    sig.oxds(stg,:,:) = sig_oxds(stg,1:nphs,:);
end

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
drawnow;

%% define and check functions for mcmc

cal0 = cal;

prsampfunc= @(Ni)    priorsamp(bnds.mat, Ni);
priorfunc = @(model) prior(model, bnds.mat);
likefunc  = @(model) likefrommodel(model, cal0, oxds, Temp, Pres, stages, hasolv, haspxn, hasplg, exp, sig);
dhatfunc  = @(model) runmodel(model, cal0, oxds, Temp, Pres, stages, hasolv, haspxn, hasplg);

%% run catmip

[mout, Pout, ~, ~, allmodels] = catmip(priorfunc,prsampfunc,likefunc,...
    'Niter',Niter,'Nsteps',Nstep,'Parallel',logical(pllopt),'Ncores',pllopt);

%% plot catmip results

PlotTemperingSteps(allmodels,bnds.mat,vname);
mMAP = plotpdfs(mout, Pout, vname, bnds.mat);
% [mMAP] = plotcorner(mout, Pout, [], bnds.mat, 1, 1, vname);

%% plot catmip histograms and predicted data from best-fit model

cal  = model2cal(cal0, mMAP)
mdl  = dhatfunc(mMAP);

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



%% plot best fit calibration

% plot fo-fay phase diagram
var.c = zeros(100,cal.nc);
var.c(:,cal.fo ) = linspace(0,1,100);
var.c(:,cal.fay) = 1-var.c(:,cal.fo);
var.T = 1300.*ones(100,1);
var.P =  0.1.*ones(100,1);

[~,cal]  =  meltmodel(var,cal,'T');

figure();
plot(var.c(:,cal.fay),cal.Tsol,'k-',var.c(:,cal.fay),cal.Tliq,'k-',LW{:}); axis xy tight; box on;
set(gca,TL{:},FS{[1,2]})
ylabel('$T \ [^\circ$C]',TX{:},FS{[1,3]})
xlabel('fay [wt]',TX{:},FS{[1,3]})

% plot opx-cpx phase diagram
var.c = zeros(100,cal.nc);
var.c(:,cal.opx) = linspace(0,1,100);
var.c(:,cal.cpx) = 1-var.c(:,cal.opx);
var.T = 1300.*ones(100,1);
var.P =  0.1.*ones(100,1);

[~,cal]  =  meltmodel(var,cal,'T');

figure();
plot(var.c(:,cal.cpx),cal.Tsol,'k-',var.c(:,cal.cpx),cal.Tliq,'k-',LW{:}); axis xy tight; box on;
set(gca,TL{:},FS{[1,2]})
ylabel('$T \ [^\circ$C]',TX{:},FS{[1,3]})
xlabel('cpx [wt]',TX{:},FS{[1,3]})

% plot an-ab phase diagram
var.c = zeros(100,cal.nc);
var.c(:,cal.an) = linspace(0,1,100);
var.c(:,cal.ab) = 1-var.c(:,cal.an);
var.T = 1300.*ones(100,1);
var.P =  0.1.*ones(100,1);

[~,cal]  =  meltmodel(var,cal,'T');

figure();
plot(var.c(:,cal.ab),cal.Tsol,'k-',var.c(:,cal.ab),cal.Tliq,'k-',LW{:}); axis xy tight; box on;
set(gca,TL{:},FS{[1,2]})
ylabel('$T \ [^\circ$C]',TX{:},FS{[1,3]})
xlabel('ab [wt]',TX{:},FS{[1,3]})

% set bulk composition as ['fo','fa','opx','px2','cpx,'an','ab']
c0 = [0.35;0.10;0.23;0.13;0.15;0.04];  % find end-member proportions of first step
c0 = c0.'./sum(c0);

% plot partition coefficients (T-dependent)
var.c = ones(100,1)*c0;
var.T = linspace(1000,1800,100).';
var.P = linspace(1.5,1.5,100).';

[~,cal]  =  meltmodel(var,cal,'K');

figure(); clf;
semilogy(var.T,cal.K,LW{:}); axis xy tight; box on;
legend(cal.CompStr{:},FS{[1,2]},TX{:})
set(gca,TL{:},FS{[1,2]})
xlabel('$T \ [^\circ$C]',TX{:},FS{[1,3]})
ylabel('log$_{10} \ K$',TX{:},FS{[1,3]})


% plot partition coefficients (P-dependent)
var.c = ones(100,1)*c0;
var.T = linspace(1300,1300,100).';
var.P = linspace(0,4,100).';

[~,cal]  =  meltmodel(var,cal,'K');

figure(); clf;
semilogx(cal.K,var.P,LW{:}); axis ij tight; box on;
legend(cal.CompStr{:},FS{[1,2]},TX{:})
set(gca,TL{:},FS{[1,2]})
xlabel('log$_{10} \ K$',TX{:},FS{[1,3]})
ylabel('$P$ [GPa]',TX{:},FS{[1,3]})

% plot melting points, solidus, liquidus (P-dependent)
var.c = ones(100,1)*c0;
var.T = linspace(1300,1300,100).';
var.P = linspace(0,4,100).';

[~,cal]  =  meltmodel(var,cal,'T');

figure(); clf;
plot(cal.Tm,var.P,LW{:}); axis ij tight; box on; hold on
plot(cal.Tsol,var.P,'k--',LW{1},3); axis ij tight; box on;
plot(cal.Tliq,var.P,'k-.',LW{1},3); axis ij tight; box on;
legend(cal.CompStr{:},'$T_\mathrm{sol}$','$T_\mathrm{liq}$',FS{[1,2]},TX{:})
set(gca,TL{:},FS{[1,2]})
xlabel('$T \ [^\circ$C]',TX{:},FS{[1,3]})
ylabel('$P$ [GPa]',TX{:},FS{[1,3]})

% plot melt fraction, phase compositions (T-dependent)
var.c = ones(100,1)*c0;
var.T = linspace(1300,1300,100).';
var.P = linspace(0.1,0.1,100).';

[~,cal]  =  meltmodel(var,cal,'T');
var.T = linspace(mean(cal.Tsol),mean(cal.Tliq),100).';

[phs,~]  =  meltmodel(var,cal,'E');

figure(); clf;
plot(var.T,1-phs.f,LW{:}); axis xy tight; box on; hold on;
plot(var.T,phs.f,LW{:}); axis xy tight; box on;
legend('xtals','melt',FS{[1,2]},TX{:},'location','east')
set(gca,TL{:},FS{[1,2]})
xlabel('$T \ [^\circ$C]',TX{:},FS{[1,3]})
ylabel('$f$ [wt]',TX{:},FS{[1,3]})

figure(); clf;
plot(var.T,phs.cs,LW{:}); axis xy tight; box on;
legend(cal.CompStr{:},FS{[1,2]},TX{:})
set(gca,TL{:},FS{[1,2]})
xlabel('$T \ [^\circ$C]',TX{:},FS{[1,3]})
ylabel('$c^s$ [wt]',TX{:},FS{[1,3]})

figure(); clf;
plot(var.T,phs.cl,LW{:}); axis xy tight; box on;
legend(cal.CompStr{:},FS{[1,2]},TX{:})
set(gca,TL{:},FS{[1,2]})
xlabel('$T \ [^\circ$C]',TX{:},FS{[1,3]})
ylabel('$c^\ell$ [wt]',TX{:},FS{[1,3]})

figure(); clf;
plot(var.T,phs.cl*cal.oxds,LW{:}); axis xy tight; box on;
legend('SiO$_2$','Al$_2$O$_3$','FeO','MgO','CaO','Na$_2$O',FS{[1,2]},TX{:})
set(gca,TL{:},FS{[1,2]})
xlabel('$T \ [^\circ$C]',TX{:},FS{[1,3]})
ylabel('melt composition [wt\%]',TX{:},FS{[1,3]})
