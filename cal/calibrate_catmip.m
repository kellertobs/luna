%% prep workspace
clear all; close all;
addpath('../src');
addpath('util/')

% prep plotting options
TX = {'Interpreter','latex'};
TL = {'TickLabelInterpreter','latex'};
FS = {'FontSize',16,20};
MS = {'MarkerSize',8};
LW = {'LineWidth',2};
CL = {'Color',[0.4660, 0.6740, 0.1880],[0.4940, 0.1840, 0.5560],[0, 0.4470, 0.7410],[0.300, 0.2500, 0.1500],[0.80, 0.60, 0.50],[0.8500, 0.3250, 0.0980],[0.8,0.8,0.8]};
TINY = 1e-16;

%% catmip parameters - adjust this

% catmip options
Niter   = 500;            % number of iterations per catmip temperature
Nstep   = 100;            % number of steps in MCMC in catmip
pllopt  = 8;              % whether to run in parallel, number of workers

rng(10); % for reproducibility during testing, comment this for actual run

%% load calibration
cal_catmip;
cal0 = cal;

%% set up data for calibration

% load experimental data from Schmidt & Kraettli (2020), Table 3
Load_SKTable3;
nphs = 6; olv=1; opx=2; cpx=3; plg=4; spn=5; qtz=6; mlt=7; blk=8;
ncmp = 11; Si=1; Ti=2; Al=3; Cr=4; Fe=5; Mn=6; Mg=7; Ca=8; Na=9; K=10; P=11;

oxd(:,:,P ) = [];  % ignore P2O5
oxd(:,:,K ) = [];  % ignore K2O
oxd(:,:,Mn) = [];  % ignore MnO
oxd(:,:,Cr) = [];  % ignore Cr2O3

% combine opx and cpx
oxd(:,opx,:) = (phs(:,opx).*squeeze(oxd(:,opx,:)) + phs(:,cpx).*squeeze(oxd(:,cpx,:)))./(phs(:,opx) + phs(:,cpx) + TINY);
phs(:,opx  ) = phs(:,opx) + phs(:,cpx);
oxd(:,cpx,:) = []; phs(:,cpx) = [];

oxd = oxd./(sum(oxd,3)+TINY)*100;  % normalise remaining oxides to 100%

nphs = 6; olv=1; pxn=2; plg=3; spn=4; qtz=5; mlt=6; blk=7;
ncmp = 6; Si=1; Ti=2; Al=3; Fe=4; Mg=5; Ca=6; Na=7;

% assign errors for phases and oxides - please review
factor           = 2; % to adjust sigmas, maybe too small
sig_phs         = factor*0.020.*ones(size(phs));
sig_oxd          = zeros(size(oxd));
sig_oxd(oxd>=20) = factor*max(0.02*oxd(oxd>=20), 1.00);
sig_oxd(oxd< 20) = factor*max(0.10*oxd(oxd< 20), 0.20);
sig_oxd(:,pxn,:) = 2*sig_oxd(:,pxn,:);


% load data for solidus and liquidus
Psl = (0.1:0.2:4)'; % all across pressure space
[Tsol, Tliq] = solidusliquidus('johnson2021', Psl);
Tsl_sigma = 10;

% weights for datasets - [SK2020, Tsol/Tliq]
wgts = [10; 1];

%% prepare information for the parameter fitting

% select stages to include in parameters fitting
stages = [1,2,3,4,5,6,7,8,9,10];

for stg = stages
    
    % exclude results where experiments do not have one of the phases
    hasolv(stg) = (phs(stg,olv)>0);
    haspxn(stg) = (phs(stg,pxn)>0);
    hasplg(stg) = (phs(stg,plg)>0);
    hasspn(stg) = (phs(stg,spn)>0);
    hasqtz(stg) = (phs(stg,qtz)>0);

    exp.Tmp(stg,:  ) = Tmp(stg);
    exp.phs(stg,:  ) = phs(stg,:);
    exp.oxd(stg,:,:) = oxd(stg,1:nphs,:);
    
    sig.phs(stg,:  ) = sig_phs(stg,:);
    sig.oxd(stg,:,:) = sig_oxd(stg,1:nphs,:);

    % fit component fractions to bulk compositions
    c0(stg,:)  =  fitc0(squeeze(oxd(stg,blk,:)).',cal.oxd);  % find end-member proportions of first step

end

cal.c0 = c0(1,:);

% plot Harker diagrams
figure(1); clf;
plot(Tmp(phs(:,olv)>0),phs(phs(:,olv)>0,olv),'d',CL{[1,2]},MS{:},LW{1},1.5); axis xy tight; box on; hold on
plot(Tmp(phs(:,pxn)>0),phs(phs(:,pxn)>0,pxn),'v',CL{[1,3]},MS{:},LW{1},1.5);
plot(Tmp(phs(:,plg)>0),phs(phs(:,plg)>0,plg),'s',CL{[1,4]},MS{:},LW{1},1.5);
plot(Tmp(phs(:,spn)>0),phs(phs(:,spn)>0,spn),'^',CL{[1,5]},MS{:},LW{1},1.5);
plot(Tmp(phs(:,qtz)>0),phs(phs(:,qtz)>0,qtz),'x',CL{[1,6]},MS{:},LW{1},1.5);
plot(Tmp(phs(:,mlt)>0),phs(phs(:,mlt)>0,mlt),'o',CL{[1,7]},MS{:},LW{1},1.5);
set(gca,TL{:}); xlabel('$T [^\circ$C]',TX{:}); ylabel('$f_\mathrm{phs}$ [wt]',TX{:});

figure(2); clf;
subplot(2,3,1);
plot(oxd(:,blk,Si),oxd(:,blk,Ti),'k.',MS{1},12,LW{1},1.5); axis xy tight; box on; hold on
plot(oxd(phs(:,olv)>0,olv,Si),oxd(phs(:,olv)>0,olv,Ti),'d',CL{[1,2]},MS{:},LW{1},1.5);
plot(oxd(phs(:,pxn)>0,pxn,Si),oxd(phs(:,pxn)>0,pxn,Ti),'v',CL{[1,3]},MS{:},LW{1},1.5);
plot(oxd(phs(:,plg)>0,plg,Si),oxd(phs(:,plg)>0,plg,Ti),'s',CL{[1,4]},MS{:},LW{1},1.5);
plot(oxd(phs(:,mlt)>0,mlt,Si),oxd(phs(:,mlt)>0,mlt,Ti),'o',CL{[1,7]},MS{:},LW{1},1.5);
plot(cal.mnr_oxd(cal.mnr_for,Si),cal.mnr_oxd(cal.mnr_for,Ti),'d',CL{[1,2]},MS{1},12,LW{:});
plot(cal.mnr_oxd(cal.mnr_fay,Si),cal.mnr_oxd(cal.mnr_fay,Ti),'d',CL{[1,2]},MS{1},12,LW{:});
plot(cal.mnr_oxd(cal.mnr_px1,Si),cal.mnr_oxd(cal.mnr_px1,Ti),'v',CL{[1,3]},MS{1},12,LW{:});
plot(cal.mnr_oxd(cal.mnr_px2,Si),cal.mnr_oxd(cal.mnr_px2,Ti),'v',CL{[1,3]},MS{1},12,LW{:});
plot(cal.mnr_oxd(cal.mnr_px3,Si),cal.mnr_oxd(cal.mnr_px3,Ti),'v',CL{[1,3]},MS{1},12,LW{:});
plot(cal.mnr_oxd(cal.mnr_ant,Si),cal.mnr_oxd(cal.mnr_ant,Ti),'s',CL{[1,4]},MS{1},12,LW{:});
plot(cal.mnr_oxd(cal.mnr_alb,Si),cal.mnr_oxd(cal.mnr_alb,Ti),'s',CL{[1,4]},MS{1},12,LW{:});
plot(cal.oxd(cal.eut,Si),cal.oxd(cal.eut,Ti),'h',CL{[1,7]},MS{1},12,LW{:});
set(gca,TL{:}); xlabel('SiO$_2$ [wt\%]',TX{:}); ylabel('TiO$_2$ [wt\%]',TX{:});

subplot(2,3,2);
plot(oxd(:,blk,Si),oxd(:,blk,Al),'k.',MS{1},12,LW{1},1.5); axis xy tight; box on; hold on
plot(oxd(phs(:,olv)>0,olv,Si),oxd(phs(:,olv)>0,olv,Al),'d',CL{[1,2]},MS{:},LW{1},1.5);
plot(oxd(phs(:,pxn)>0,pxn,Si),oxd(phs(:,pxn)>0,pxn,Al),'v',CL{[1,3]},MS{:},LW{1},1.5);
plot(oxd(phs(:,plg)>0,plg,Si),oxd(phs(:,plg)>0,plg,Al),'s',CL{[1,4]},MS{:},LW{1},1.5);
plot(oxd(phs(:,mlt)>0,mlt,Si),oxd(phs(:,mlt)>0,mlt,Al),'o',CL{[1,7]},MS{:},LW{1},1.5);
plot(cal.mnr_oxd(cal.mnr_for,Si),cal.mnr_oxd(cal.mnr_for,Al),'d',CL{[1,2]},MS{1},12,LW{:});
plot(cal.mnr_oxd(cal.mnr_fay,Si),cal.mnr_oxd(cal.mnr_fay,Al),'d',CL{[1,2]},MS{1},12,LW{:});
plot(cal.mnr_oxd(cal.mnr_px1,Si),cal.mnr_oxd(cal.mnr_px1,Al),'v',CL{[1,3]},MS{1},12,LW{:});
plot(cal.mnr_oxd(cal.mnr_px2,Si),cal.mnr_oxd(cal.mnr_px2,Al),'v',CL{[1,3]},MS{1},12,LW{:});
plot(cal.mnr_oxd(cal.mnr_px3,Si),cal.mnr_oxd(cal.mnr_px3,Al),'v',CL{[1,3]},MS{1},12,LW{:});
plot(cal.mnr_oxd(cal.mnr_ant,Si),cal.mnr_oxd(cal.mnr_ant,Al),'s',CL{[1,4]},MS{1},12,LW{:});
plot(cal.mnr_oxd(cal.mnr_alb,Si),cal.mnr_oxd(cal.mnr_alb,Al),'s',CL{[1,4]},MS{1},12,LW{:});
plot(cal.oxd(cal.eut,Si),cal.oxd(cal.eut,Al),'h',CL{[1,7]},MS{1},12,LW{:});
set(gca,TL{:}); xlabel('SiO$_2$ [wt\%]',TX{:}); ylabel('Al$_2$O$_3$ [wt\%]',TX{:});

subplot(2,3,3);
plot(oxd(:,blk,Si),oxd(:,blk,Fe),'k.',MS{1},12,LW{1},1.5); axis xy tight; box on; hold on
plot(oxd(phs(:,olv)>0,olv,Si),oxd(phs(:,olv)>0,olv,Fe),'d',CL{[1,2]},MS{:},LW{1},1.5);
plot(oxd(phs(:,pxn)>0,pxn,Si),oxd(phs(:,pxn)>0,pxn,Fe),'v',CL{[1,3]},MS{:},LW{1},1.5);
plot(oxd(phs(:,plg)>0,plg,Si),oxd(phs(:,plg)>0,plg,Fe),'s',CL{[1,4]},MS{:},LW{1},1.5);
plot(oxd(phs(:,mlt)>0,mlt,Si),oxd(phs(:,mlt)>0,mlt,Fe),'o',CL{[1,7]},MS{:},LW{1},1.5);
plot(cal.mnr_oxd(cal.mnr_for,Si),cal.mnr_oxd(cal.mnr_for,Fe),'d',CL{[1,2]},MS{1},12,LW{:});
plot(cal.mnr_oxd(cal.mnr_fay,Si),cal.mnr_oxd(cal.mnr_fay,Fe),'d',CL{[1,2]},MS{1},12,LW{:});
plot(cal.mnr_oxd(cal.mnr_px1,Si),cal.mnr_oxd(cal.mnr_px1,Fe),'v',CL{[1,3]},MS{1},12,LW{:});
plot(cal.mnr_oxd(cal.mnr_px2,Si),cal.mnr_oxd(cal.mnr_px2,Fe),'v',CL{[1,3]},MS{1},12,LW{:});
plot(cal.mnr_oxd(cal.mnr_px3,Si),cal.mnr_oxd(cal.mnr_px3,Fe),'v',CL{[1,3]},MS{1},12,LW{:});
plot(cal.mnr_oxd(cal.mnr_ant,Si),cal.mnr_oxd(cal.mnr_ant,Fe),'s',CL{[1,4]},MS{1},12,LW{:});
plot(cal.mnr_oxd(cal.mnr_alb,Si),cal.mnr_oxd(cal.mnr_alb,Fe),'s',CL{[1,4]},MS{1},12,LW{:});
plot(cal.oxd(cal.eut,Si),cal.oxd(cal.eut,Fe),'h',CL{[1,7]},MS{1},12,LW{:});
set(gca,TL{:}); xlabel('SiO$_2$ [wt\%]',TX{:}); ylabel('FeO [wt\%]',TX{:});

subplot(2,3,4);
plot(oxd(:,blk,Si),oxd(:,blk,Mg),'k.',MS{1},12,LW{1},1.5); axis xy tight; box on; hold on
plot(oxd(phs(:,olv)>0,olv,Si),oxd(phs(:,olv)>0,olv,Mg),'d',CL{[1,2]},MS{:},LW{1},1.5);
plot(oxd(phs(:,pxn)>0,pxn,Si),oxd(phs(:,pxn)>0,pxn,Mg),'v',CL{[1,3]},MS{:},LW{1},1.5);
plot(oxd(phs(:,plg)>0,plg,Si),oxd(phs(:,plg)>0,plg,Mg),'s',CL{[1,4]},MS{:},LW{1},1.5);
plot(oxd(phs(:,mlt)>0,mlt,Si),oxd(phs(:,mlt)>0,mlt,Mg),'o',CL{[1,7]},MS{:},LW{1},1.5);
plot(cal.mnr_oxd(cal.mnr_for,Si),cal.mnr_oxd(cal.mnr_for,Mg),'d',CL{[1,2]},MS{1},12,LW{:});
plot(cal.mnr_oxd(cal.mnr_fay,Si),cal.mnr_oxd(cal.mnr_fay,Mg),'d',CL{[1,2]},MS{1},12,LW{:});
plot(cal.mnr_oxd(cal.mnr_px1,Si),cal.mnr_oxd(cal.mnr_px1,Mg),'v',CL{[1,3]},MS{1},12,LW{:});
plot(cal.mnr_oxd(cal.mnr_px2,Si),cal.mnr_oxd(cal.mnr_px2,Mg),'v',CL{[1,3]},MS{1},12,LW{:});
plot(cal.mnr_oxd(cal.mnr_px3,Si),cal.mnr_oxd(cal.mnr_px3,Mg),'v',CL{[1,3]},MS{1},12,LW{:});
plot(cal.mnr_oxd(cal.mnr_ant,Si),cal.mnr_oxd(cal.mnr_ant,Mg),'s',CL{[1,4]},MS{1},12,LW{:});
plot(cal.mnr_oxd(cal.mnr_alb,Si),cal.mnr_oxd(cal.mnr_alb,Mg),'s',CL{[1,4]},MS{1},12,LW{:});
plot(cal.oxd(cal.eut,Si),cal.oxd(cal.eut,Mg),'h',CL{[1,7]},MS{1},12,LW{:});
set(gca,TL{:}); xlabel('SiO$_2$ [wt\%]',TX{:}); ylabel('MgO [wt\%]',TX{:});

subplot(2,3,5);
plot(oxd(:,blk,Si),oxd(:,blk,Ca),'k.',MS{1},12,LW{1},1.5); axis xy tight; box on; hold on
plot(oxd(phs(:,olv)>0,olv,Si),oxd(phs(:,olv)>0,olv,Ca),'d',CL{[1,2]},MS{:},LW{1},1.5);
plot(oxd(phs(:,pxn)>0,pxn,Si),oxd(phs(:,pxn)>0,pxn,Ca),'v',CL{[1,3]},MS{:},LW{1},1.5);
plot(oxd(phs(:,plg)>0,plg,Si),oxd(phs(:,plg)>0,plg,Ca),'s',CL{[1,4]},MS{:},LW{1},1.5);
plot(oxd(phs(:,mlt)>0,mlt,Si),oxd(phs(:,mlt)>0,mlt,Ca),'o',CL{[1,7]},MS{:},LW{1},1.5);
plot(cal.mnr_oxd(cal.mnr_for,Si),cal.mnr_oxd(cal.mnr_for,Ca),'d',CL{[1,2]},MS{1},12,LW{:});
plot(cal.mnr_oxd(cal.mnr_fay,Si),cal.mnr_oxd(cal.mnr_fay,Ca),'d',CL{[1,2]},MS{1},12,LW{:});
plot(cal.mnr_oxd(cal.mnr_px1,Si),cal.mnr_oxd(cal.mnr_px1,Ca),'v',CL{[1,3]},MS{1},12,LW{:});
plot(cal.mnr_oxd(cal.mnr_px2,Si),cal.mnr_oxd(cal.mnr_px2,Ca),'v',CL{[1,3]},MS{1},12,LW{:});
plot(cal.mnr_oxd(cal.mnr_px3,Si),cal.mnr_oxd(cal.mnr_px3,Ca),'v',CL{[1,3]},MS{1},12,LW{:});
plot(cal.mnr_oxd(cal.mnr_ant,Si),cal.mnr_oxd(cal.mnr_ant,Ca),'s',CL{[1,4]},MS{1},12,LW{:});
plot(cal.mnr_oxd(cal.mnr_alb,Si),cal.mnr_oxd(cal.mnr_alb,Ca),'s',CL{[1,4]},MS{1},12,LW{:});
plot(cal.oxd(cal.eut,Si),cal.oxd(cal.eut,Ca),'h',CL{[1,7]},MS{1},12,LW{:});
set(gca,TL{:}); xlabel('SiO$_2$ [wt\%]',TX{:}); ylabel('CaO [wt\%]',TX{:});

subplot(2,3,6);
plot(oxd(:,blk,Si),oxd(:,blk,Na),'k.',MS{1},12,LW{1},1.5); axis xy tight; box on; hold on
plot(oxd(phs(:,olv)>0,olv,Si),oxd(phs(:,olv)>0,olv,Na),'d',CL{[1,2]},MS{:},LW{1},1.5);
plot(oxd(phs(:,pxn)>0,pxn,Si),oxd(phs(:,pxn)>0,pxn,Na),'v',CL{[1,3]},MS{:},LW{1},1.5);
plot(oxd(phs(:,plg)>0,plg,Si),oxd(phs(:,plg)>0,plg,Na),'s',CL{[1,4]},MS{:},LW{1},1.5);
plot(oxd(phs(:,mlt)>0,mlt,Si),oxd(phs(:,mlt)>0,mlt,Na),'o',CL{[1,7]},MS{:},LW{1},1.5);
plot(cal.mnr_oxd(cal.mnr_for,Si),cal.mnr_oxd(cal.mnr_for,Na),'d',CL{[1,2]},MS{1},12,LW{:});
plot(cal.mnr_oxd(cal.mnr_fay,Si),cal.mnr_oxd(cal.mnr_fay,Na),'d',CL{[1,2]},MS{1},12,LW{:});
plot(cal.mnr_oxd(cal.mnr_px1,Si),cal.mnr_oxd(cal.mnr_px1,Na),'v',CL{[1,3]},MS{1},12,LW{:});
plot(cal.mnr_oxd(cal.mnr_px2,Si),cal.mnr_oxd(cal.mnr_px2,Na),'v',CL{[1,3]},MS{1},12,LW{:});
plot(cal.mnr_oxd(cal.mnr_px3,Si),cal.mnr_oxd(cal.mnr_px3,Na),'v',CL{[1,3]},MS{1},12,LW{:});
plot(cal.mnr_oxd(cal.mnr_ant,Si),cal.mnr_oxd(cal.mnr_ant,Na),'s',CL{[1,4]},MS{1},12,LW{:});
plot(cal.mnr_oxd(cal.mnr_alb,Si),cal.mnr_oxd(cal.mnr_alb,Na),'s',CL{[1,4]},MS{1},12,LW{:});
plot(cal.oxd(cal.eut,Si),cal.oxd(cal.eut,Na),'h',CL{[1,7]},MS{1},12,LW{:});
set(gca,TL{:}); xlabel('SiO$_2$ [wt\%]',TX{:}); ylabel('Na$_2$O [wt\%]',TX{:});
drawnow;

%% define and check functions for inversion

prsampfunc = @(Ni)    priorsamp(bnds.mat, Ni);
priorfunc  = @(model) prior(model, bnds.mat);
likefunc   = @(model) likefrommodel(model, cal0, c0, Tmp, Prs, stages, hasolv, haspxn, hasplg, hasspn, hasqtz, exp, sig, Tsol, Tliq, Psl, Tsl_sigma, wgts);
dhatfunc   = @(model) runmodel(model, cal0, c0, Tmp, Prs, stages, hasolv, haspxn, hasplg, hasspn, hasqtz, Psl);

%% run catmip

[mout, Pout, ~, ~, allmodels] = catmip(priorfunc,prsampfunc,likefunc,...
    'Niter',Niter,'Nsteps',Nstep,'Parallel',logical(pllopt),'Ncores',pllopt);

%% plot catmip results

% PlotTmperingSteps(allmodels,bnds.mat,vname);
mMAP = plotpdfs(mout, Pout, vname, bnds.mat);
% [mMAP] = plotcorner(mout, Pout, [], bnds.mat, 1, 1, vname);

% plot catmip histograms and predicted data from best-fit model

cal  = model2cal(cal0, mMAP)
mdl  = dhatfunc(mMAP);

% mdl = runmodel([], cal0, c0, Tmp, Prs, stages, hasolv, haspxn, hasplg, hasspn, hasqtz, Psl);;

figure(3); clf;
plot(exp.Tmp(hasolv),exp.phs(hasolv,olv),'d',CL{[1,8]},MS{:},LW{1},1.5); axis xy tight; box on; hold on
plot(exp.Tmp(haspxn),exp.phs(haspxn,pxn),'v',CL{[1,8]},MS{:},LW{1},1.5);
plot(exp.Tmp(hasplg),exp.phs(hasplg,plg),'s',CL{[1,8]},MS{:},LW{1},1.5);
plot(exp.Tmp(hasspn),exp.phs(hasspn,spn),'^',CL{[1,8]},MS{:},LW{1},1.5);
plot(exp.Tmp(hasqtz),exp.phs(hasqtz,qtz),'x',CL{[1,8]},MS{:},LW{1},1.5);
plot(exp.Tmp,exp.phs(:,mlt),'o',CL{[1,8]},MS{:},LW{1},1.5);
plot(exp.Tmp(hasolv),mdl.phs(hasolv,olv),'d',CL{[1,2]},MS{:},LW{1},1.5);
plot(exp.Tmp(haspxn),mdl.phs(haspxn,pxn),'v',CL{[1,3]},MS{:},LW{1},1.5);
plot(exp.Tmp(hasplg),mdl.phs(hasplg,plg),'s',CL{[1,4]},MS{:},LW{1},1.5);
plot(exp.Tmp(hasspn),mdl.phs(hasspn,spn),'^',CL{[1,5]},MS{:},LW{1},1.5);
plot(exp.Tmp(hasqtz),mdl.phs(hasqtz,qtz),'x',CL{[1,6]},MS{:},LW{1},1.5);
plot(exp.Tmp,mdl.phs(:,mlt),'o',CL{[1,7]},MS{:},LW{1},1.5);
set(gca,TL{:}); xlabel('$T [^\circ$C]',TX{:}); ylabel('$f_\mathrm{phs}$ [wt]',TX{:});

figure(4); clf;
subplot(2,3,1);
plot(exp.oxd(hasolv,olv,Si),exp.oxd(hasolv,olv,Ti),'d',CL{[1,8]},MS{:},LW{1},1.5); axis xy tight; box on; hold on
plot(exp.oxd(haspxn,pxn,Si),exp.oxd(haspxn,pxn,Ti),'v',CL{[1,8]},MS{:},LW{1},1.5);
plot(exp.oxd(hasplg,plg,Si),exp.oxd(hasplg,plg,Ti),'s',CL{[1,8]},MS{:},LW{1},1.5);
plot(exp.oxd(:,mlt,Si),exp.oxd(:,mlt,Ti),'o',CL{[1,8]},MS{:},LW{1},1.5);
plot(mdl.oxd(hasolv,olv,Si),mdl.oxd(hasolv,olv,Ti),'d',CL{[1,2]},MS{:},LW{1},1.5);
plot(mdl.oxd(haspxn,pxn,Si),mdl.oxd(haspxn,pxn,Ti),'v',CL{[1,3]},MS{:},LW{1},1.5);
plot(mdl.oxd(hasplg,plg,Si),mdl.oxd(hasplg,plg,Ti),'s',CL{[1,4]},MS{:},LW{1},1.5);
plot(mdl.oxd(:,mlt,Si),mdl.oxd(:,mlt,Ti),'o',CL{[1,7]},MS{:},LW{1},1.5);
plot(cal.mnr_oxd(cal.mnr_for,Si),cal.mnr_oxd(cal.mnr_for,Ti),'d',CL{[1,2]},MS{1},12,LW{:});
plot(cal.mnr_oxd(cal.mnr_fay,Si),cal.mnr_oxd(cal.mnr_fay,Ti),'d',CL{[1,2]},MS{1},12,LW{:});
plot(cal.mnr_oxd(cal.mnr_px1,Si),cal.mnr_oxd(cal.mnr_px1,Ti),'v',CL{[1,3]},MS{1},12,LW{:});
plot(cal.mnr_oxd(cal.mnr_px2,Si),cal.mnr_oxd(cal.mnr_px2,Ti),'v',CL{[1,3]},MS{1},12,LW{:});
plot(cal.mnr_oxd(cal.mnr_px3,Si),cal.mnr_oxd(cal.mnr_px3,Ti),'v',CL{[1,3]},MS{1},12,LW{:});
plot(cal.mnr_oxd(cal.mnr_ant,Si),cal.mnr_oxd(cal.mnr_ant,Ti),'s',CL{[1,4]},MS{1},12,LW{:});
plot(cal.mnr_oxd(cal.mnr_alb,Si),cal.mnr_oxd(cal.mnr_alb,Ti),'s',CL{[1,4]},MS{1},12,LW{:});
plot(cal.oxd(cal.eut,Si),cal.oxd(cal.eut,Ti),'h',CL{[1,7]},MS{1},12,LW{:});
set(gca,TL{:}); xlabel('SiO$_2$ [wt\%]',TX{:}); ylabel('TiO$_2$ [wt\%]',TX{:});

subplot(2,3,2);
plot(exp.oxd(hasolv,olv,Si),exp.oxd(hasolv,olv,Al),'d',CL{[1,8]},MS{:},LW{1},1.5); axis xy tight; box on; hold on
plot(exp.oxd(haspxn,pxn,Si),exp.oxd(haspxn,pxn,Al),'v',CL{[1,8]},MS{:},LW{1},1.5);
plot(exp.oxd(hasplg,plg,Si),exp.oxd(hasplg,plg,Al),'s',CL{[1,8]},MS{:},LW{1},1.5);
plot(exp.oxd(:,mlt,Si),exp.oxd(:,mlt,Al),'o',CL{[1,8]},MS{:},LW{1},1.5);
plot(mdl.oxd(hasolv,olv,Si),mdl.oxd(hasolv,olv,Al),'d',CL{[1,2]},MS{:},LW{1},1.5);
plot(mdl.oxd(haspxn,pxn,Si),mdl.oxd(haspxn,pxn,Al),'v',CL{[1,3]},MS{:},LW{1},1.5);
plot(mdl.oxd(hasplg,plg,Si),mdl.oxd(hasplg,plg,Al),'s',CL{[1,4]},MS{:},LW{1},1.5);
plot(mdl.oxd(:,mlt,Si),mdl.oxd(:,mlt,Al),'o',CL{[1,7]},MS{:},LW{1},1.5);
plot(cal.mnr_oxd(cal.mnr_for,Si),cal.mnr_oxd(cal.mnr_for,Al),'d',CL{[1,2]},MS{1},12,LW{:});
plot(cal.mnr_oxd(cal.mnr_fay,Si),cal.mnr_oxd(cal.mnr_fay,Al),'d',CL{[1,2]},MS{1},12,LW{:});
plot(cal.mnr_oxd(cal.mnr_px1,Si),cal.mnr_oxd(cal.mnr_px1,Al),'v',CL{[1,3]},MS{1},12,LW{:});
plot(cal.mnr_oxd(cal.mnr_px2,Si),cal.mnr_oxd(cal.mnr_px2,Al),'v',CL{[1,3]},MS{1},12,LW{:});
plot(cal.mnr_oxd(cal.mnr_px3,Si),cal.mnr_oxd(cal.mnr_px3,Al),'v',CL{[1,3]},MS{1},12,LW{:});
plot(cal.mnr_oxd(cal.mnr_ant,Si),cal.mnr_oxd(cal.mnr_ant,Al),'s',CL{[1,4]},MS{1},12,LW{:});
plot(cal.mnr_oxd(cal.mnr_alb,Si),cal.mnr_oxd(cal.mnr_alb,Al),'s',CL{[1,4]},MS{1},12,LW{:});
plot(cal.oxd(cal.eut,Si),cal.oxd(cal.eut,Al),'h',CL{[1,7]},MS{1},12,LW{:});
set(gca,TL{:}); xlabel('SiO$_2$ [wt\%]',TX{:}); ylabel('Al$_2$O$_3$ [wt\%]',TX{:});

subplot(2,3,3);
plot(exp.oxd(hasolv,olv,Si),exp.oxd(hasolv,olv,Fe),'d',CL{[1,8]},MS{:},LW{1},1.5); axis xy tight; box on; hold on
plot(exp.oxd(haspxn,pxn,Si),exp.oxd(haspxn,pxn,Fe),'v',CL{[1,8]},MS{:},LW{1},1.5);
plot(exp.oxd(hasplg,plg,Si),exp.oxd(hasplg,plg,Fe),'s',CL{[1,8]},MS{:},LW{1},1.5);
plot(exp.oxd(:,mlt,Si),exp.oxd(:,mlt,Fe),'o',CL{[1,8]},MS{:},LW{1},1.5);
plot(mdl.oxd(hasolv,olv,Si),mdl.oxd(hasolv,olv,Fe),'d',CL{[1,2]},MS{:},LW{1},1.5);
plot(mdl.oxd(haspxn,pxn,Si),mdl.oxd(haspxn,pxn,Fe),'v',CL{[1,3]},MS{:},LW{1},1.5);
plot(mdl.oxd(hasplg,plg,Si),mdl.oxd(hasplg,plg,Fe),'s',CL{[1,4]},MS{:},LW{1},1.5);
plot(mdl.oxd(:,mlt,Si),mdl.oxd(:,mlt,Fe),'o',CL{[1,7]},MS{:},LW{1},1.5);
plot(cal.mnr_oxd(cal.mnr_for,Si),cal.mnr_oxd(cal.mnr_for,Fe),'d',CL{[1,2]},MS{1},12,LW{:});
plot(cal.mnr_oxd(cal.mnr_fay,Si),cal.mnr_oxd(cal.mnr_fay,Fe),'d',CL{[1,2]},MS{1},12,LW{:});
plot(cal.mnr_oxd(cal.mnr_px1,Si),cal.mnr_oxd(cal.mnr_px1,Fe),'v',CL{[1,3]},MS{1},12,LW{:});
plot(cal.mnr_oxd(cal.mnr_px2,Si),cal.mnr_oxd(cal.mnr_px2,Fe),'v',CL{[1,3]},MS{1},12,LW{:});
plot(cal.mnr_oxd(cal.mnr_px3,Si),cal.mnr_oxd(cal.mnr_px3,Fe),'v',CL{[1,3]},MS{1},12,LW{:});
plot(cal.mnr_oxd(cal.mnr_ant,Si),cal.mnr_oxd(cal.mnr_ant,Fe),'s',CL{[1,4]},MS{1},12,LW{:});
plot(cal.mnr_oxd(cal.mnr_alb,Si),cal.mnr_oxd(cal.mnr_alb,Fe),'s',CL{[1,4]},MS{1},12,LW{:});
plot(cal.oxd(cal.eut,Si),cal.oxd(cal.eut,Fe),'h',CL{[1,7]},MS{1},12,LW{:});
set(gca,TL{:}); xlabel('SiO$_2$ [wt\%]',TX{:}); ylabel('FeO [wt\%]',TX{:});

subplot(2,3,4);
plot(exp.oxd(hasolv,olv,Si),exp.oxd(hasolv,olv,Mg),'d',CL{[1,8]},MS{:},LW{1},1.5); axis xy tight; box on; hold on
plot(exp.oxd(haspxn,pxn,Si),exp.oxd(haspxn,pxn,Mg),'v',CL{[1,8]},MS{:},LW{1},1.5);
plot(exp.oxd(hasplg,plg,Si),exp.oxd(hasplg,plg,Mg),'s',CL{[1,8]},MS{:},LW{1},1.5);
plot(exp.oxd(:,mlt,Si),exp.oxd(:,mlt,Mg),'o',CL{[1,8]},MS{:},LW{1},1.5);
plot(mdl.oxd(hasolv,olv,Si),mdl.oxd(hasolv,olv,Mg),'d',CL{[1,2]},MS{:},LW{1},1.5);
plot(mdl.oxd(haspxn,pxn,Si),mdl.oxd(haspxn,pxn,Mg),'v',CL{[1,3]},MS{:},LW{1},1.5);
plot(mdl.oxd(hasplg,plg,Si),mdl.oxd(hasplg,plg,Mg),'s',CL{[1,4]},MS{:},LW{1},1.5);
plot(mdl.oxd(:,mlt,Si),mdl.oxd(:,mlt,Mg),'o',CL{[1,7]},MS{:},LW{1},1.5);
plot(cal.mnr_oxd(cal.mnr_for,Si),cal.mnr_oxd(cal.mnr_for,Mg),'d',CL{[1,2]},MS{1},12,LW{:});
plot(cal.mnr_oxd(cal.mnr_fay,Si),cal.mnr_oxd(cal.mnr_fay,Mg),'d',CL{[1,2]},MS{1},12,LW{:});
plot(cal.mnr_oxd(cal.mnr_px1,Si),cal.mnr_oxd(cal.mnr_px1,Mg),'v',CL{[1,3]},MS{1},12,LW{:});
plot(cal.mnr_oxd(cal.mnr_px2,Si),cal.mnr_oxd(cal.mnr_px2,Mg),'v',CL{[1,3]},MS{1},12,LW{:});
plot(cal.mnr_oxd(cal.mnr_px3,Si),cal.mnr_oxd(cal.mnr_px3,Mg),'v',CL{[1,3]},MS{1},12,LW{:});
plot(cal.mnr_oxd(cal.mnr_ant,Si),cal.mnr_oxd(cal.mnr_ant,Mg),'s',CL{[1,4]},MS{1},12,LW{:});
plot(cal.mnr_oxd(cal.mnr_alb,Si),cal.mnr_oxd(cal.mnr_alb,Mg),'s',CL{[1,4]},MS{1},12,LW{:});
plot(cal.oxd(cal.eut,Si),cal.oxd(cal.eut,Mg),'h',CL{[1,7]},MS{1},12,LW{:});
set(gca,TL{:}); xlabel('SiO$_2$ [wt\%]',TX{:}); ylabel('MgO [wt\%]',TX{:});

subplot(2,3,5);
plot(exp.oxd(hasolv,olv,Si),exp.oxd(hasolv,olv,Ca),'d',CL{[1,8]},MS{:},LW{1},1.5); axis xy tight; box on; hold on
plot(exp.oxd(haspxn,pxn,Si),exp.oxd(haspxn,pxn,Ca),'v',CL{[1,8]},MS{:},LW{1},1.5);
plot(exp.oxd(hasplg,plg,Si),exp.oxd(hasplg,plg,Ca),'s',CL{[1,8]},MS{:},LW{1},1.5);
plot(exp.oxd(:,mlt,Si),exp.oxd(:,mlt,Ca),'o',CL{[1,8]},MS{:},LW{1},1.5);
plot(mdl.oxd(hasolv,olv,Si),mdl.oxd(hasolv,olv,Ca),'d',CL{[1,2]},MS{:},LW{1},1.5);
plot(mdl.oxd(haspxn,pxn,Si),mdl.oxd(haspxn,pxn,Ca),'v',CL{[1,3]},MS{:},LW{1},1.5);
plot(mdl.oxd(hasplg,plg,Si),mdl.oxd(hasplg,plg,Ca),'s',CL{[1,4]},MS{:},LW{1},1.5);
plot(mdl.oxd(:,mlt,Si),mdl.oxd(:,mlt,Ca),'o',CL{[1,7]},MS{:},LW{1},1.5);
plot(cal.mnr_oxd(cal.mnr_for,Si),cal.mnr_oxd(cal.mnr_for,Ca),'d',CL{[1,2]},MS{1},12,LW{:});
plot(cal.mnr_oxd(cal.mnr_fay,Si),cal.mnr_oxd(cal.mnr_fay,Ca),'d',CL{[1,2]},MS{1},12,LW{:});
plot(cal.mnr_oxd(cal.mnr_px1,Si),cal.mnr_oxd(cal.mnr_px1,Ca),'v',CL{[1,3]},MS{1},12,LW{:});
plot(cal.mnr_oxd(cal.mnr_px2,Si),cal.mnr_oxd(cal.mnr_px2,Ca),'v',CL{[1,3]},MS{1},12,LW{:});
plot(cal.mnr_oxd(cal.mnr_px3,Si),cal.mnr_oxd(cal.mnr_px3,Ca),'v',CL{[1,3]},MS{1},12,LW{:});
plot(cal.mnr_oxd(cal.mnr_ant,Si),cal.mnr_oxd(cal.mnr_ant,Ca),'s',CL{[1,4]},MS{1},12,LW{:});
plot(cal.mnr_oxd(cal.mnr_alb,Si),cal.mnr_oxd(cal.mnr_alb,Ca),'s',CL{[1,4]},MS{1},12,LW{:});
plot(cal.oxd(cal.eut,Si),cal.oxd(cal.eut,Ca),'h',CL{[1,7]},MS{1},12,LW{:});
set(gca,TL{:}); xlabel('SiO$_2$ [wt\%]',TX{:}); ylabel('CaO [wt\%]',TX{:});

subplot(2,3,6);
plot(exp.oxd(hasolv,olv,Si),exp.oxd(hasolv,olv,Na),'d',CL{[1,8]},MS{:},LW{1},1.5); axis xy tight; box on; hold on
plot(exp.oxd(haspxn,pxn,Si),exp.oxd(haspxn,pxn,Na),'v',CL{[1,8]},MS{:},LW{1},1.5);
plot(exp.oxd(hasplg,plg,Si),exp.oxd(hasplg,plg,Na),'s',CL{[1,8]},MS{:},LW{1},1.5);
plot(exp.oxd(:,mlt,Si),exp.oxd(:,mlt,Na),'o',CL{[1,8]},MS{:},LW{1},1.5);
plot(mdl.oxd(hasolv,olv,Si),mdl.oxd(hasolv,olv,Na),'d',CL{[1,2]},MS{:},LW{1},1.5);
plot(mdl.oxd(haspxn,pxn,Si),mdl.oxd(haspxn,pxn,Na),'v',CL{[1,3]},MS{:},LW{1},1.5);
plot(mdl.oxd(hasplg,plg,Si),mdl.oxd(hasplg,plg,Na),'s',CL{[1,4]},MS{:},LW{1},1.5);
plot(mdl.oxd(:,mlt,Si),mdl.oxd(:,mlt,Na),'o',CL{[1,7]},MS{:},LW{1},1.5);
plot(cal.mnr_oxd(cal.mnr_for,Si),cal.mnr_oxd(cal.mnr_for,Na),'d',CL{[1,2]},MS{1},12,LW{:});
plot(cal.mnr_oxd(cal.mnr_fay,Si),cal.mnr_oxd(cal.mnr_fay,Na),'d',CL{[1,2]},MS{1},12,LW{:});
plot(cal.mnr_oxd(cal.mnr_px1,Si),cal.mnr_oxd(cal.mnr_px1,Na),'v',CL{[1,3]},MS{1},12,LW{:});
plot(cal.mnr_oxd(cal.mnr_px2,Si),cal.mnr_oxd(cal.mnr_px2,Na),'v',CL{[1,3]},MS{1},12,LW{:});
plot(cal.mnr_oxd(cal.mnr_px3,Si),cal.mnr_oxd(cal.mnr_px3,Na),'v',CL{[1,3]},MS{1},12,LW{:});
plot(cal.mnr_oxd(cal.mnr_ant,Si),cal.mnr_oxd(cal.mnr_ant,Na),'s',CL{[1,4]},MS{1},12,LW{:});
plot(cal.mnr_oxd(cal.mnr_alb,Si),cal.mnr_oxd(cal.mnr_alb,Na),'s',CL{[1,4]},MS{1},12,LW{:});
plot(cal.oxd(cal.eut,Si),cal.oxd(cal.eut,Na),'h',CL{[1,7]},MS{1},12,LW{:});
set(gca,TL{:}); xlabel('SiO$_2$ [wt\%]',TX{:}); ylabel('Na$_2$O [wt\%]',TX{:});



%% plot best fit calibration

% plot fo-fay phase diagram
var.c = zeros(100,cal.nc);
var.c(:,cal.for) = linspace(0,1,100);
var.c(:,cal.fay) = 1-var.c(:,cal.for);
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
var.c(:,cal.ant) = linspace(0,1,100);
var.c(:,cal.eut) = 1-var.c(:,cal.ant);
var.T = 1300.*ones(100,1);
var.P =  0.1.*ones(100,1);

[~,cal]  =  meltmodel(var,cal,'T');

figure();
plot(var.c(:,cal.eut),cal.Tsol,'k-',var.c(:,cal.eut),cal.Tliq,'k-',LW{:}); axis xy tight; box on;
set(gca,TL{:},FS{[1,2]})
ylabel('$T \ [^\circ$C]',TX{:},FS{[1,3]})
xlabel('ab [wt]',TX{:},FS{[1,3]})

% plot partition coefficients (T-dependent)
var.c = ones(100,1)*c0(1,:);
var.T = linspace(1000,1800,100).';
var.P = linspace(1.5,1.5,100).';

[~,cal]  =  meltmodel(var,cal,'K');

figure(); clf;
semilogy(var.T,cal.K,LW{:}); axis xy tight; box on;
legend(cal.cmpStr{:},FS{[1,2]},TX{:})
set(gca,TL{:},FS{[1,2]})
xlabel('$T \ [^\circ$C]',TX{:},FS{[1,3]})
ylabel('log$_{10} \ K$',TX{:},FS{[1,3]})


% plot partition coefficients (P-dependent)
var.c = ones(100,1)*c0(1,:);
var.T = linspace(1300,1300,100).';
var.P = linspace(0,4,100).';

[~,cal]  =  meltmodel(var,cal,'K');

figure(); clf;
semilogx(cal.K,var.P,LW{:}); axis ij tight; box on;
legend(cal.cmpStr{:},FS{[1,2]},TX{:})
set(gca,TL{:},FS{[1,2]})
xlabel('log$_{10} \ K$',TX{:},FS{[1,3]})
ylabel('$P$ [GPa]',TX{:},FS{[1,3]})

% plot melting points, solidus, liquidus (P-dependent)
var.c = ones(100,1)*c0(1,:);
var.T = linspace(1300,1300,100).';
var.P = linspace(0,4,100).';

[~,cal]  =  meltmodel(var,cal,'T');

figure(); clf;
plot(cal.Tm,var.P,LW{:}); axis ij tight; box on; hold on
plot(cal.Tsol,var.P,'k--',LW{1},3); axis ij tight; box on;
plot(cal.Tliq,var.P,'k-.',LW{1},3); axis ij tight; box on;
plot(Tsol, Psl, 'k+', Tliq, Psl, 'kx')
legend(cal.cmpStr{:},'$T_\mathrm{sol}$','$T_\mathrm{liq}$',FS{[1,2]},TX{:})
set(gca,TL{:},FS{[1,2]})
xlabel('$T \ [^\circ$C]',TX{:},FS{[1,3]})
ylabel('$P$ [GPa]',TX{:},FS{[1,3]})

% plot melt fraction, phase compositions (T-dependent)
var.c = ones(100,1)*c0(1,:);
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
legend(cal.cmpStr{:},FS{[1,2]},TX{:})
set(gca,TL{:},FS{[1,2]})
xlabel('$T \ [^\circ$C]',TX{:},FS{[1,3]})
ylabel('$c^s$ [wt]',TX{:},FS{[1,3]})

figure(); clf;
plot(var.T,phs.cl,LW{:}); axis xy tight; box on;
legend(cal.cmpStr{:},FS{[1,2]},TX{:})
set(gca,TL{:},FS{[1,2]})
xlabel('$T \ [^\circ$C]',TX{:},FS{[1,3]})
ylabel('$c^\ell$ [wt]',TX{:},FS{[1,3]})

figure(); clf;
plot(var.T,phs.cl*cal.oxd,LW{:}); axis xy tight; box on;
legend('SiO$_2$','Al$_2$O$_3$','FeO','MgO','CaO','Na$_2$O',FS{[1,2]},TX{:})
set(gca,TL{:},FS{[1,2]})
xlabel('$T \ [^\circ$C]',TX{:},FS{[1,3]})
ylabel('melt composition [wt\%]',TX{:},FS{[1,3]})
