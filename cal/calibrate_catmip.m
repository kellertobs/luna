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
Nstep   = 200;            % number of steps in MCMC in catmip
pllopt  = 8;              % whether to run in parallel, number of workers

rng(10); % for reproducibility during testing, comment this for actual run

%% load calibration
% cal_catmip_4c;
cal_luna_4c;
cal0 = cal;

%% set up data for calibration

% load experimental data from Schmidt & Kraettli (2020), Table 3
Load_SKTable3;

% remove stages
phs(11,:) = [];  oxd(11,:,:) = [];
% phs( 7,:) = [];  oxd( 9,:,:) = [];
% phs( 7,:) = [];  oxd( 7,:,:) = [];
phs( 5,:) = [];  oxd( 5,:,:) = [];

nphs = 7; olv=1; opx=2; cpx=3; plg=4; ilm=5; qtz=6; mlt=7; blk=8;
ncmp = 11; Si=1; Ti=2; Al=3; Cr=4; Fe=5; Mn=6; Mg=7; Ca=8; Na=9; K=10; P=11;

oxd(:,:,P ) = [];  % ignore P2O5
oxd(:,:,K ) = [];  % ignore K2O
oxd(:,:,Mn) = [];  % ignore MnO
oxd(:,:,Cr) = [];  % ignore Cr2O3

phs = phs./(sum(phs,2)+TINY);      % normalise remaining phases to 1
oxd = oxd./(sum(oxd,3)+TINY)*100;  % normalise remaining oxides to 100%

nphs = 7; olv=1; opx=2; cpx=3; plg=4; ilm=5; qtz=6; mlt=7; blk=8; xtl=9;
ncmp = 7; Si=1; Ti=2; Al=3; Fe=4; Mg=5; Ca=6; Na=7;

% assign errors for phases and oxides - please review
factor           = 2; % to adjust sigmas, maybe too small
sig_phs         = factor*0.020.*ones(size(phs));
sig_oxd          = zeros(size(oxd));
sig_oxd(oxd>=20) = factor*max(0.02*oxd(oxd>=20), 1.00);
sig_oxd(oxd< 20) = factor*max(0.10*oxd(oxd< 20), 0.20);

% load data for solidus and liquidus
Psl = (0.1:0.2:4)'; % all across pressure space
[Tsol, Tliq] = solidusliquidus('johnson2021', Psl);
Tsl_sigma = 5;

% weights for datasets - [SK2020, Tsol/Tliq]
wgts = [1; 1];

%% prepare information for the parameter fitting
stages = 1:size(oxd,1);

for stg = stages
    
    % get fraction and oxide composition of xtal assemblage
    oxd(stg,xtl,:) = squeeze(phs(stg,1:nphs-1))*squeeze(oxd(stg,1:nphs-1,:))./sum(squeeze(phs(stg,1:nphs-1)));
    oxd(stg,blk,:) = squeeze(phs(stg,1:nphs-0))*squeeze(oxd(stg,1:nphs-0,:))./sum(squeeze(phs(stg,1:nphs-0)));


    % exclude results where experiments do not have one of the phases
    hasolv(stg) = (phs(stg,olv)>0);
    hasopx(stg) = (phs(stg,opx)>0);
    hascpx(stg) = (phs(stg,cpx)>0);
    hasplg(stg) = (phs(stg,plg)>0);
    hasilm(stg) = (phs(stg,ilm)>0);
    hasqtz(stg) = (phs(stg,qtz)>0);
    hasmlt(stg) = (phs(stg,mlt)>0);

    exp.Tmp(stg,:  ) = Tmp(stg);
    exp.phs(stg,:  ) = phs(stg,:);
    exp.oxd(stg,:,:) = oxd(stg,1:nphs,:);
    
    sig.phs(stg,:  ) = sig_phs(stg,:);
    sig.oxd(stg,:,:) = sig_oxd(stg,1:nphs,:);

    % fit component fractions to bulk compositions
    c0(stg,:)  =  fitc0(squeeze(oxd(stg,blk,:)).',cal.cmp_oxd);  % find end-member proportions of first step

end

cal.c0 = c0(1,:);

% plot Harker diagrams
figure(); clf;
plot(Tmp(hasolv),phs(hasolv,olv),'d',CL{[1,2]},MS{:},LW{1},1.5); axis xy tight; box on; hold on
plot(Tmp(hasopx),phs(hasopx,opx),'^',CL{[1,3]},MS{:},LW{1},1.5);
plot(Tmp(hascpx),phs(hascpx,cpx),'v',CL{[1,3]},MS{:},LW{1},1.5);
plot(Tmp(hasplg),phs(hasplg,plg),'s',CL{[1,4]},MS{:},LW{1},1.5);
plot(Tmp(hasilm),phs(hasilm,ilm),'+',CL{[1,5]},MS{:},LW{1},1.5);
plot(Tmp(hasqtz),phs(hasqtz,qtz),'x',CL{[1,6]},MS{:},LW{1},1.5);
plot(Tmp(hasmlt),phs(hasmlt,mlt),'o',CL{[1,7]},MS{:},LW{1},1.5);
set(gca,TL{:}); xlabel('$T [^\circ$C]',TX{:}); ylabel('$f_\mathrm{phs}$ [wt]',TX{:});

figure(); clf;
subplot(2,3,1);
plot(oxd(:,blk,Si),oxd(:,blk,Ti),'k.',MS{1},15,LW{1},1.5); axis xy tight; box on; hold on
plot(oxd(:,xtl,Si),oxd(:,xtl,Ti),'b.',MS{1},15,LW{1},1.5);
plot(oxd(:,mlt,Si),oxd(:,mlt,Ti),'r.',MS{1},15,LW{1},1.5);
plot(oxd(hasolv,olv,Si),oxd(hasolv,olv,Ti),'d',CL{[1,2]},MS{:},LW{1},1.5);
plot(oxd(hasopx,opx,Si),oxd(hasopx,opx,Ti),'^',CL{[1,3]},MS{:},LW{1},1.5); axis xy tight; box on; hold on
plot(oxd(hascpx,cpx,Si),oxd(hascpx,cpx,Ti),'v',CL{[1,3]},MS{:},LW{1},1.5);
plot(oxd(hasplg,plg,Si),oxd(hasplg,plg,Ti),'s',CL{[1,4]},MS{:},LW{1},1.5);
plot(oxd(hasilm,ilm,Si),oxd(hasilm,ilm,Ti),'+',CL{[1,5]},MS{:},LW{1},1.5);
plot(cal.mnr_oxd(cal.for,Si),cal.mnr_oxd(cal.for,Ti),'d',CL{[1,2]},MS{1},12,LW{:});
plot(cal.mnr_oxd(cal.fay,Si),cal.mnr_oxd(cal.fay,Ti),'d',CL{[1,2]},MS{1},12,LW{:});
plot(cal.mnr_oxd(cal.ens,Si),cal.mnr_oxd(cal.ens,Ti),'^',CL{[1,3]},MS{1},12,LW{:});
plot(cal.mnr_oxd(cal.hyp,Si),cal.mnr_oxd(cal.hyp,Ti),'^',CL{[1,3]},MS{1},12,LW{:});
plot(cal.mnr_oxd(cal.aug,Si),cal.mnr_oxd(cal.aug,Ti),'v',CL{[1,3]},MS{1},12,LW{:});
plot(cal.mnr_oxd(cal.pig,Si),cal.mnr_oxd(cal.pig,Ti),'v',CL{[1,3]},MS{1},12,LW{:});
plot(cal.mnr_oxd(cal.ant,Si),cal.mnr_oxd(cal.ant,Ti),'s',CL{[1,4]},MS{1},12,LW{:});
plot(cal.mnr_oxd(cal.alb,Si),cal.mnr_oxd(cal.alb,Ti),'s',CL{[1,4]},MS{1},12,LW{:});
plot(cal.mnr_oxd(cal.ilm,Si),cal.mnr_oxd(cal.ilm,Ti),'+',CL{[1,5]},MS{1},12,LW{:});
for ic=1:cal.ncmp
    plot(cal.cmp_oxd(ic,Si),cal.cmp_oxd(ic,Ti),'o',CL{1},[0,0,0]+ic/2/cal.ncmp,MS{1},12,LW{:});
end
set(gca,TL{:}); xlabel('SiO$_2$ [wt\%]',TX{:}); ylabel('TiO$_2$ [wt\%]',TX{:});

subplot(2,3,2);
plot(oxd(:,blk,Si),oxd(:,blk,Al),'k.',MS{1},15,LW{1},1.5); axis xy tight; box on; hold on
plot(oxd(:,xtl,Si),oxd(:,xtl,Al),'b.',MS{1},15,LW{1},1.5);
plot(oxd(:,mlt,Si),oxd(:,mlt,Al),'r.',MS{1},15,LW{1},1.5);
plot(oxd(hasolv,olv,Si),oxd(hasolv,olv,Al),'d',CL{[1,2]},MS{:},LW{1},1.5);
plot(oxd(hasopx,opx,Si),oxd(hasopx,opx,Al),'^',CL{[1,3]},MS{:},LW{1},1.5); axis xy tight; box on; hold on
plot(oxd(hascpx,cpx,Si),oxd(hascpx,cpx,Al),'v',CL{[1,3]},MS{:},LW{1},1.5);
plot(oxd(hasplg,plg,Si),oxd(hasplg,plg,Al),'s',CL{[1,4]},MS{:},LW{1},1.5);
plot(cal.mnr_oxd(cal.for,Si),cal.mnr_oxd(cal.for,Al),'d',CL{[1,2]},MS{1},12,LW{:});
plot(cal.mnr_oxd(cal.fay,Si),cal.mnr_oxd(cal.fay,Al),'d',CL{[1,2]},MS{1},12,LW{:});
plot(cal.mnr_oxd(cal.ens,Si),cal.mnr_oxd(cal.ens,Al),'^',CL{[1,3]},MS{1},12,LW{:});
plot(cal.mnr_oxd(cal.hyp,Si),cal.mnr_oxd(cal.hyp,Al),'^',CL{[1,3]},MS{1},12,LW{:});
plot(cal.mnr_oxd(cal.aug,Si),cal.mnr_oxd(cal.aug,Al),'v',CL{[1,3]},MS{1},12,LW{:});
plot(cal.mnr_oxd(cal.pig,Si),cal.mnr_oxd(cal.pig,Al),'v',CL{[1,3]},MS{1},12,LW{:});
plot(cal.mnr_oxd(cal.ant,Si),cal.mnr_oxd(cal.ant,Al),'s',CL{[1,4]},MS{1},12,LW{:});
plot(cal.mnr_oxd(cal.alb,Si),cal.mnr_oxd(cal.alb,Al),'s',CL{[1,4]},MS{1},12,LW{:});
for ic=1:cal.ncmp
    plot(cal.cmp_oxd(ic,Si),cal.cmp_oxd(ic,Al),'o',CL{1},[0,0,0]+ic/2/cal.ncmp,MS{1},12,LW{:});
end
set(gca,TL{:}); xlabel('SiO$_2$ [wt\%]',TX{:}); ylabel('Al$_2$O$_3$ [wt\%]',TX{:});

subplot(2,3,3);
plot(oxd(:,blk,Si),oxd(:,blk,Fe),'k.',MS{1},15,LW{1},1.5); axis xy tight; box on; hold on
plot(oxd(:,xtl,Si),oxd(:,xtl,Fe),'b.',MS{1},15,LW{1},1.5);
plot(oxd(:,mlt,Si),oxd(:,mlt,Fe),'r.',MS{1},15,LW{1},1.5);
plot(oxd(hasolv,olv,Si),oxd(hasolv,olv,Fe),'d',CL{[1,2]},MS{:},LW{1},1.5);
plot(oxd(hasopx,opx,Si),oxd(hasopx,opx,Fe),'^',CL{[1,3]},MS{:},LW{1},1.5); axis xy tight; box on; hold on
plot(oxd(hascpx,cpx,Si),oxd(hascpx,cpx,Fe),'v',CL{[1,3]},MS{:},LW{1},1.5);
plot(oxd(hasplg,plg,Si),oxd(hasplg,plg,Fe),'s',CL{[1,4]},MS{:},LW{1},1.5);
plot(cal.mnr_oxd(cal.for,Si),cal.mnr_oxd(cal.for,Fe),'d',CL{[1,2]},MS{1},12,LW{:});
plot(cal.mnr_oxd(cal.fay,Si),cal.mnr_oxd(cal.fay,Fe),'d',CL{[1,2]},MS{1},12,LW{:});
plot(cal.mnr_oxd(cal.ens,Si),cal.mnr_oxd(cal.ens,Fe),'^',CL{[1,3]},MS{1},12,LW{:});
plot(cal.mnr_oxd(cal.hyp,Si),cal.mnr_oxd(cal.hyp,Fe),'^',CL{[1,3]},MS{1},12,LW{:});
plot(cal.mnr_oxd(cal.aug,Si),cal.mnr_oxd(cal.aug,Fe),'v',CL{[1,3]},MS{1},12,LW{:});
plot(cal.mnr_oxd(cal.pig,Si),cal.mnr_oxd(cal.pig,Fe),'v',CL{[1,3]},MS{1},12,LW{:});
plot(cal.mnr_oxd(cal.ant,Si),cal.mnr_oxd(cal.ant,Fe),'s',CL{[1,4]},MS{1},12,LW{:});
plot(cal.mnr_oxd(cal.alb,Si),cal.mnr_oxd(cal.alb,Fe),'s',CL{[1,4]},MS{1},12,LW{:});
for ic=1:cal.ncmp
    plot(cal.cmp_oxd(ic,Si),cal.cmp_oxd(ic,Fe),'o',CL{1},[0,0,0]+ic/2/cal.ncmp,MS{1},12,LW{:});
end
set(gca,TL{:}); xlabel('SiO$_2$ [wt\%]',TX{:}); ylabel('FeO [wt\%]',TX{:});

subplot(2,3,4);
plot(oxd(:,blk,Si),oxd(:,blk,Mg),'k.',MS{1},15,LW{1},1.5); axis xy tight; box on; hold on
plot(oxd(:,xtl,Si),oxd(:,xtl,Mg),'b.',MS{1},15,LW{1},1.5);
plot(oxd(:,mlt,Si),oxd(:,mlt,Mg),'r.',MS{1},15,LW{1},1.5);
plot(oxd(hasolv,olv,Si),oxd(hasolv,olv,Mg),'d',CL{[1,2]},MS{:},LW{1},1.5);
plot(oxd(hasopx,opx,Si),oxd(hasopx,opx,Mg),'^',CL{[1,3]},MS{:},LW{1},1.5); axis xy tight; box on; hold on
plot(oxd(hascpx,cpx,Si),oxd(hascpx,cpx,Mg),'v',CL{[1,3]},MS{:},LW{1},1.5);
plot(oxd(hasplg,plg,Si),oxd(hasplg,plg,Mg),'s',CL{[1,4]},MS{:},LW{1},1.5);
plot(cal.mnr_oxd(cal.for,Si),cal.mnr_oxd(cal.for,Mg),'d',CL{[1,2]},MS{1},12,LW{:});
plot(cal.mnr_oxd(cal.fay,Si),cal.mnr_oxd(cal.fay,Mg),'d',CL{[1,2]},MS{1},12,LW{:});
plot(cal.mnr_oxd(cal.ens,Si),cal.mnr_oxd(cal.ens,Mg),'^',CL{[1,3]},MS{1},12,LW{:});
plot(cal.mnr_oxd(cal.hyp,Si),cal.mnr_oxd(cal.hyp,Mg),'^',CL{[1,3]},MS{1},12,LW{:});
plot(cal.mnr_oxd(cal.aug,Si),cal.mnr_oxd(cal.aug,Mg),'v',CL{[1,3]},MS{1},12,LW{:});
plot(cal.mnr_oxd(cal.pig,Si),cal.mnr_oxd(cal.pig,Mg),'v',CL{[1,3]},MS{1},12,LW{:});
plot(cal.mnr_oxd(cal.ant,Si),cal.mnr_oxd(cal.ant,Mg),'s',CL{[1,4]},MS{1},12,LW{:});
plot(cal.mnr_oxd(cal.alb,Si),cal.mnr_oxd(cal.alb,Mg),'s',CL{[1,4]},MS{1},12,LW{:});
for ic=1:cal.ncmp
    plot(cal.cmp_oxd(ic,Si),cal.cmp_oxd(ic,Mg),'o',CL{1},[0,0,0]+ic/2/cal.ncmp,MS{1},12,LW{:});
end
set(gca,TL{:}); xlabel('SiO$_2$ [wt\%]',TX{:}); ylabel('MgO [wt\%]',TX{:});

subplot(2,3,5);
plot(oxd(:,blk,Si),oxd(:,blk,Ca),'k.',MS{1},15,LW{1},1.5); axis xy tight; box on; hold on
plot(oxd(:,xtl,Si),oxd(:,xtl,Ca),'b.',MS{1},15,LW{1},1.5);
plot(oxd(:,mlt,Si),oxd(:,mlt,Ca),'r.',MS{1},15,LW{1},1.5);
plot(oxd(hasolv,olv,Si),oxd(hasolv,olv,Ca),'d',CL{[1,2]},MS{:},LW{1},1.5);
plot(oxd(hasopx,opx,Si),oxd(hasopx,opx,Ca),'^',CL{[1,3]},MS{:},LW{1},1.5); axis xy tight; box on; hold on
plot(oxd(hascpx,cpx,Si),oxd(hascpx,cpx,Ca),'v',CL{[1,3]},MS{:},LW{1},1.5);
plot(oxd(hasplg,plg,Si),oxd(hasplg,plg,Ca),'s',CL{[1,4]},MS{:},LW{1},1.5);
plot(cal.mnr_oxd(cal.for,Si),cal.mnr_oxd(cal.for,Ca),'d',CL{[1,2]},MS{1},12,LW{:});
plot(cal.mnr_oxd(cal.fay,Si),cal.mnr_oxd(cal.fay,Ca),'d',CL{[1,2]},MS{1},12,LW{:});
plot(cal.mnr_oxd(cal.ens,Si),cal.mnr_oxd(cal.ens,Ca),'^',CL{[1,3]},MS{1},12,LW{:});
plot(cal.mnr_oxd(cal.hyp,Si),cal.mnr_oxd(cal.hyp,Ca),'^',CL{[1,3]},MS{1},12,LW{:});
plot(cal.mnr_oxd(cal.aug,Si),cal.mnr_oxd(cal.aug,Ca),'v',CL{[1,3]},MS{1},12,LW{:});
plot(cal.mnr_oxd(cal.pig,Si),cal.mnr_oxd(cal.pig,Ca),'v',CL{[1,3]},MS{1},12,LW{:});
plot(cal.mnr_oxd(cal.ant,Si),cal.mnr_oxd(cal.ant,Ca),'s',CL{[1,4]},MS{1},12,LW{:});
plot(cal.mnr_oxd(cal.alb,Si),cal.mnr_oxd(cal.alb,Ca),'s',CL{[1,4]},MS{1},12,LW{:});
for ic=1:cal.ncmp
    plot(cal.cmp_oxd(ic,Si),cal.cmp_oxd(ic,Ca),'o',CL{1},[0,0,0]+ic/2/cal.ncmp,MS{1},12,LW{:});
end
set(gca,TL{:}); xlabel('SiO$_2$ [wt\%]',TX{:}); ylabel('CaO [wt\%]',TX{:});

subplot(2,3,6);
plot(oxd(:,blk,Si),oxd(:,blk,Na),'k.',MS{1},15,LW{1},1.5); axis xy tight; box on; hold on
plot(oxd(:,xtl,Si),oxd(:,xtl,Na),'b.',MS{1},15,LW{1},1.5);
plot(oxd(:,mlt,Si),oxd(:,mlt,Na),'r.',MS{1},15,LW{1},1.5);
plot(oxd(hasolv,olv,Si),oxd(hasolv,olv,Na),'d',CL{[1,2]},MS{:},LW{1},1.5);
plot(oxd(hasopx,opx,Si),oxd(hasopx,opx,Na),'^',CL{[1,3]},MS{:},LW{1},1.5); axis xy tight; box on; hold on
plot(oxd(hascpx,cpx,Si),oxd(hascpx,cpx,Na),'v',CL{[1,3]},MS{:},LW{1},1.5);
plot(oxd(hasplg,plg,Si),oxd(hasplg,plg,Na),'s',CL{[1,4]},MS{:},LW{1},1.5);
plot(cal.mnr_oxd(cal.for,Si),cal.mnr_oxd(cal.for,Na),'d',CL{[1,2]},MS{1},12,LW{:});
plot(cal.mnr_oxd(cal.fay,Si),cal.mnr_oxd(cal.fay,Na),'d',CL{[1,2]},MS{1},12,LW{:});
plot(cal.mnr_oxd(cal.ens,Si),cal.mnr_oxd(cal.ens,Na),'^',CL{[1,3]},MS{1},12,LW{:});
plot(cal.mnr_oxd(cal.hyp,Si),cal.mnr_oxd(cal.hyp,Na),'^',CL{[1,3]},MS{1},12,LW{:});
plot(cal.mnr_oxd(cal.aug,Si),cal.mnr_oxd(cal.aug,Na),'v',CL{[1,3]},MS{1},12,LW{:});
plot(cal.mnr_oxd(cal.pig,Si),cal.mnr_oxd(cal.pig,Na),'v',CL{[1,3]},MS{1},12,LW{:});
plot(cal.mnr_oxd(cal.ant,Si),cal.mnr_oxd(cal.ant,Na),'s',CL{[1,4]},MS{1},12,LW{:});
plot(cal.mnr_oxd(cal.alb,Si),cal.mnr_oxd(cal.alb,Na),'s',CL{[1,4]},MS{1},12,LW{:});
for ic=1:cal.ncmp
    plot(cal.cmp_oxd(ic,Si),cal.cmp_oxd(ic,Na),'o',CL{1},[0,0,0]+ic/2/cal.ncmp,MS{1},12,LW{:});
end
set(gca,TL{:}); xlabel('SiO$_2$ [wt\%]',TX{:}); ylabel('Na$_2$O [wt\%]',TX{:});
drawnow;

%% define and check functions for inversion

prsampfunc = @(Ni)    priorsamp(bnds.mat, Ni);
priorfunc  = @(model) prior(model, bnds.mat);
likefunc   = @(model) likefrommodel(model, cal0, c0, Tmp, Prs, stages, exp, sig, Tsol, Tliq, Psl, Tsl_sigma, wgts);
dhatfunc   = @(model) runmodel(model, cal0, c0, Tmp, Prs, stages, Psl);

% %% run catmip
% 
% [mout, Pout, ~, ~, allmodels] = catmip(priorfunc,prsampfunc,likefunc,...
%     'Niter',Niter,'Nsteps',Nstep,'Parallel',logical(pllopt),'Ncores',pllopt);
% 
% %% plot catmip results
% 
% % PlotTmperingSteps(allmodels,bnds.mat,vname);
% mMAP = plotpdfs(mout, Pout, vname, bnds.mat);
% % [mMAP] = plotcorner(mout, Pout, [], bnds.mat, 1, 1, vname);
% 
% % plot catmip histograms and predicted data from best-fit model
% 
% cal  = model2cal(cal0, mMAP)
% mdl  = dhatfunc(mMAP);

cal_luna_4c; cal0 = cal;
mdl = runmodel([], cal0, c0, Tmp, Prs, stages, Psl);
L   = likelihood (mdl, exp, sig, Tsol, Tliq, Tsl_sigma, wgts)

figure(3); clf;
plot(exp.Tmp(hasolv),exp.phs(hasolv,olv),'d',CL{[1,8]},MS{:},LW{1},1.5); axis xy tight; box on; hold on
plot(exp.Tmp(hasopx),exp.phs(hasopx,opx),'^',CL{[1,8]},MS{:},LW{1},1.5);
plot(exp.Tmp(hascpx),exp.phs(hascpx,cpx),'v',CL{[1,8]},MS{:},LW{1},1.5);
plot(exp.Tmp(hasplg),exp.phs(hasplg,plg),'s',CL{[1,8]},MS{:},LW{1},1.5);
plot(exp.Tmp(hasilm),exp.phs(hasilm,ilm),'+',CL{[1,8]},MS{:},LW{1},1.5);
plot(exp.Tmp(hasqtz),exp.phs(hasqtz,qtz),'x',CL{[1,8]},MS{:},LW{1},1.5);
plot(exp.Tmp(hasmlt),exp.phs(hasmlt,mlt),'o',CL{[1,8]},MS{:},LW{1},1.5);
plot(exp.Tmp(hasolv),mdl.phs(hasolv,olv),'d',CL{[1,2]},MS{:},LW{1},1.5);
plot(exp.Tmp(hasopx),mdl.phs(hasopx,opx),'^',CL{[1,3]},MS{:},LW{1},1.5);
plot(exp.Tmp(hascpx),mdl.phs(hascpx,cpx),'v',CL{[1,3]},MS{:},LW{1},1.5);
plot(exp.Tmp(hasplg),mdl.phs(hasplg,plg),'s',CL{[1,4]},MS{:},LW{1},1.5);
plot(exp.Tmp(hasilm),mdl.phs(hasilm,ilm),'+',CL{[1,5]},MS{:},LW{1},1.5);
plot(exp.Tmp(hasqtz),mdl.phs(hasqtz,qtz),'x',CL{[1,6]},MS{:},LW{1},1.5);
plot(exp.Tmp(hasmlt),mdl.phs(hasmlt,mlt),'o',CL{[1,7]},MS{:},LW{1},1.5);
set(gca,TL{:}); xlabel('$T [^\circ$C]',TX{:}); ylabel('$f_\mathrm{phs}$ [wt]',TX{:});

figure(4); clf;
subplot(2,3,1);
plot(exp.oxd(hasolv,olv,Si),exp.oxd(hasolv,olv,Ti),'d',CL{[1,8]},MS{:},LW{1},1.5); axis xy tight; box on; hold on
plot(exp.oxd(hasopx,opx,Si),exp.oxd(hasopx,opx,Ti),'^',CL{[1,8]},MS{:},LW{1},1.5);
plot(exp.oxd(hascpx,cpx,Si),exp.oxd(hascpx,cpx,Ti),'v',CL{[1,8]},MS{:},LW{1},1.5);
plot(exp.oxd(hasplg,plg,Si),exp.oxd(hasplg,plg,Ti),'s',CL{[1,8]},MS{:},LW{1},1.5);
plot(exp.oxd(hasmlt,mlt,Si),exp.oxd(hasmlt,mlt,Ti),'.',CL{[1,8]},MS{1},15,LW{1},1.5);
plot(mdl.oxd(hasolv,olv,Si),mdl.oxd(hasolv,olv,Ti),'d',CL{[1,2]},MS{:},LW{1},1.5);
plot(mdl.oxd(hasopx,opx,Si),mdl.oxd(hasopx,opx,Ti),'^',CL{[1,3]},MS{:},LW{1},1.5);
plot(mdl.oxd(hascpx,cpx,Si),mdl.oxd(hascpx,cpx,Ti),'v',CL{[1,3]},MS{:},LW{1},1.5);
plot(mdl.oxd(hasplg,plg,Si),mdl.oxd(hasplg,plg,Ti),'s',CL{[1,4]},MS{:},LW{1},1.5);
plot(mdl.oxd(hasmlt,mlt,Si),mdl.oxd(hasmlt,mlt,Ti),'.',CL{[1,7]},MS{1},15,LW{1},1.5);
plot(cal.mnr_oxd(cal.for,Si),cal.mnr_oxd(cal.for,Ti),'d',CL{[1,2]},MS{1},12,LW{:});
plot(cal.mnr_oxd(cal.fay,Si),cal.mnr_oxd(cal.fay,Ti),'d',CL{[1,2]},MS{1},12,LW{:});
plot(cal.mnr_oxd(cal.ens,Si),cal.mnr_oxd(cal.ens,Ti),'^',CL{[1,3]},MS{1},12,LW{:});
plot(cal.mnr_oxd(cal.hyp,Si),cal.mnr_oxd(cal.hyp,Ti),'^',CL{[1,3]},MS{1},12,LW{:});
plot(cal.mnr_oxd(cal.aug,Si),cal.mnr_oxd(cal.aug,Ti),'v',CL{[1,3]},MS{1},12,LW{:});
plot(cal.mnr_oxd(cal.pig,Si),cal.mnr_oxd(cal.pig,Ti),'v',CL{[1,3]},MS{1},12,LW{:});
plot(cal.mnr_oxd(cal.ant,Si),cal.mnr_oxd(cal.ant,Ti),'s',CL{[1,4]},MS{1},12,LW{:});
plot(cal.mnr_oxd(cal.alb,Si),cal.mnr_oxd(cal.alb,Ti),'s',CL{[1,4]},MS{1},12,LW{:});
set(gca,TL{:}); xlabel('SiO$_2$ [wt\%]',TX{:}); ylabel('TiO$_2$ [wt\%]',TX{:});

subplot(2,3,2);
plot(exp.oxd(hasolv,olv,Si),exp.oxd(hasolv,olv,Al),'d',CL{[1,8]},MS{:},LW{1},1.5); axis xy tight; box on; hold on
plot(exp.oxd(hasopx,opx,Si),exp.oxd(hasopx,opx,Al),'^',CL{[1,8]},MS{:},LW{1},1.5);
plot(exp.oxd(hascpx,cpx,Si),exp.oxd(hascpx,cpx,Al),'v',CL{[1,8]},MS{:},LW{1},1.5);
plot(exp.oxd(hasplg,plg,Si),exp.oxd(hasplg,plg,Al),'s',CL{[1,8]},MS{:},LW{1},1.5);
plot(exp.oxd(hasmlt,mlt,Si),exp.oxd(hasmlt,mlt,Al),'.',CL{[1,8]},MS{1},15,LW{1},1.5);
plot(mdl.oxd(hasolv,olv,Si),mdl.oxd(hasolv,olv,Al),'d',CL{[1,2]},MS{:},LW{1},1.5);
plot(mdl.oxd(hasopx,opx,Si),mdl.oxd(hasopx,opx,Al),'^',CL{[1,3]},MS{:},LW{1},1.5);
plot(mdl.oxd(hascpx,cpx,Si),mdl.oxd(hascpx,cpx,Al),'v',CL{[1,3]},MS{:},LW{1},1.5);
plot(mdl.oxd(hasplg,plg,Si),mdl.oxd(hasplg,plg,Al),'s',CL{[1,4]},MS{:},LW{1},1.5);
plot(mdl.oxd(hasmlt,mlt,Si),mdl.oxd(hasmlt,mlt,Al),'.',CL{[1,7]},MS{1},15,LW{1},1.5);
plot(cal.mnr_oxd(cal.for,Si),cal.mnr_oxd(cal.for,Al),'d',CL{[1,2]},MS{1},12,LW{:});
plot(cal.mnr_oxd(cal.fay,Si),cal.mnr_oxd(cal.fay,Al),'d',CL{[1,2]},MS{1},12,LW{:});
plot(cal.mnr_oxd(cal.ens,Si),cal.mnr_oxd(cal.ens,Al),'^',CL{[1,3]},MS{1},12,LW{:});
plot(cal.mnr_oxd(cal.hyp,Si),cal.mnr_oxd(cal.hyp,Al),'^',CL{[1,3]},MS{1},12,LW{:});
plot(cal.mnr_oxd(cal.aug,Si),cal.mnr_oxd(cal.aug,Al),'v',CL{[1,3]},MS{1},12,LW{:});
plot(cal.mnr_oxd(cal.pig,Si),cal.mnr_oxd(cal.pig,Al),'v',CL{[1,3]},MS{1},12,LW{:});
plot(cal.mnr_oxd(cal.ant,Si),cal.mnr_oxd(cal.ant,Al),'s',CL{[1,4]},MS{1},12,LW{:});
plot(cal.mnr_oxd(cal.alb,Si),cal.mnr_oxd(cal.alb,Al),'s',CL{[1,4]},MS{1},12,LW{:});
set(gca,TL{:}); xlabel('SiO$_2$ [wt\%]',TX{:}); ylabel('Al$_2$O$_3$ [wt\%]',TX{:});

subplot(2,3,3);
plot(exp.oxd(hasolv,olv,Si),exp.oxd(hasolv,olv,Fe),'d',CL{[1,8]},MS{:},LW{1},1.5); axis xy tight; box on; hold on
plot(exp.oxd(hasopx,opx,Si),exp.oxd(hasopx,opx,Fe),'^',CL{[1,8]},MS{:},LW{1},1.5);
plot(exp.oxd(hascpx,cpx,Si),exp.oxd(hascpx,cpx,Fe),'v',CL{[1,8]},MS{:},LW{1},1.5);
plot(exp.oxd(hasplg,plg,Si),exp.oxd(hasplg,plg,Fe),'s',CL{[1,8]},MS{:},LW{1},1.5);
plot(exp.oxd(hasmlt,mlt,Si),exp.oxd(hasmlt,mlt,Fe),'.',CL{[1,8]},MS{1},15,LW{1},1.5);
plot(mdl.oxd(hasolv,olv,Si),mdl.oxd(hasolv,olv,Fe),'d',CL{[1,2]},MS{:},LW{1},1.5);
plot(mdl.oxd(hasopx,opx,Si),mdl.oxd(hasopx,opx,Fe),'^',CL{[1,3]},MS{:},LW{1},1.5);
plot(mdl.oxd(hascpx,cpx,Si),mdl.oxd(hascpx,cpx,Fe),'v',CL{[1,3]},MS{:},LW{1},1.5);
plot(mdl.oxd(hasplg,plg,Si),mdl.oxd(hasplg,plg,Fe),'s',CL{[1,4]},MS{:},LW{1},1.5);
plot(mdl.oxd(hasmlt,mlt,Si),mdl.oxd(hasmlt,mlt,Fe),'.',CL{[1,7]},MS{1},15,LW{1},1.5);
plot(cal.mnr_oxd(cal.for,Si),cal.mnr_oxd(cal.for,Fe),'d',CL{[1,2]},MS{1},12,LW{:});
plot(cal.mnr_oxd(cal.fay,Si),cal.mnr_oxd(cal.fay,Fe),'d',CL{[1,2]},MS{1},12,LW{:});
plot(cal.mnr_oxd(cal.ens,Si),cal.mnr_oxd(cal.ens,Fe),'^',CL{[1,3]},MS{1},12,LW{:});
plot(cal.mnr_oxd(cal.hyp,Si),cal.mnr_oxd(cal.hyp,Fe),'^',CL{[1,3]},MS{1},12,LW{:});
plot(cal.mnr_oxd(cal.aug,Si),cal.mnr_oxd(cal.aug,Fe),'v',CL{[1,3]},MS{1},12,LW{:});
plot(cal.mnr_oxd(cal.pig,Si),cal.mnr_oxd(cal.pig,Fe),'v',CL{[1,3]},MS{1},12,LW{:});
plot(cal.mnr_oxd(cal.ant,Si),cal.mnr_oxd(cal.ant,Fe),'s',CL{[1,4]},MS{1},12,LW{:});
plot(cal.mnr_oxd(cal.alb,Si),cal.mnr_oxd(cal.alb,Fe),'s',CL{[1,4]},MS{1},12,LW{:});
set(gca,TL{:}); xlabel('SiO$_2$ [wt\%]',TX{:}); ylabel('FeO [wt\%]',TX{:});

subplot(2,3,4);
plot(exp.oxd(hasolv,olv,Si),exp.oxd(hasolv,olv,Mg),'d',CL{[1,8]},MS{:},LW{1},1.5); axis xy tight; box on; hold on
plot(exp.oxd(hasopx,opx,Si),exp.oxd(hasopx,opx,Mg),'^',CL{[1,8]},MS{:},LW{1},1.5);
plot(exp.oxd(hascpx,cpx,Si),exp.oxd(hascpx,cpx,Mg),'v',CL{[1,8]},MS{:},LW{1},1.5);
plot(exp.oxd(hasplg,plg,Si),exp.oxd(hasplg,plg,Mg),'s',CL{[1,8]},MS{:},LW{1},1.5);
plot(exp.oxd(hasmlt,mlt,Si),exp.oxd(hasmlt,mlt,Mg),'.',CL{[1,8]},MS{1},15,LW{1},1.5);
plot(mdl.oxd(hasolv,olv,Si),mdl.oxd(hasolv,olv,Mg),'d',CL{[1,2]},MS{:},LW{1},1.5);
plot(mdl.oxd(hasopx,opx,Si),mdl.oxd(hasopx,opx,Mg),'^',CL{[1,3]},MS{:},LW{1},1.5);
plot(mdl.oxd(hascpx,cpx,Si),mdl.oxd(hascpx,cpx,Mg),'v',CL{[1,3]},MS{:},LW{1},1.5);
plot(mdl.oxd(hasplg,plg,Si),mdl.oxd(hasplg,plg,Mg),'s',CL{[1,4]},MS{:},LW{1},1.5);
plot(mdl.oxd(hasmlt,mlt,Si),mdl.oxd(hasmlt,mlt,Mg),'.',CL{[1,7]},MS{1},15,LW{1},1.5);
plot(cal.mnr_oxd(cal.for,Si),cal.mnr_oxd(cal.for,Mg),'d',CL{[1,2]},MS{1},12,LW{:});
plot(cal.mnr_oxd(cal.fay,Si),cal.mnr_oxd(cal.fay,Mg),'d',CL{[1,2]},MS{1},12,LW{:});
plot(cal.mnr_oxd(cal.ens,Si),cal.mnr_oxd(cal.ens,Mg),'^',CL{[1,3]},MS{1},12,LW{:});
plot(cal.mnr_oxd(cal.hyp,Si),cal.mnr_oxd(cal.hyp,Mg),'^',CL{[1,3]},MS{1},12,LW{:});
plot(cal.mnr_oxd(cal.aug,Si),cal.mnr_oxd(cal.aug,Mg),'v',CL{[1,3]},MS{1},12,LW{:});
plot(cal.mnr_oxd(cal.pig,Si),cal.mnr_oxd(cal.pig,Mg),'v',CL{[1,3]},MS{1},12,LW{:});
plot(cal.mnr_oxd(cal.ant,Si),cal.mnr_oxd(cal.ant,Mg),'s',CL{[1,4]},MS{1},12,LW{:});
plot(cal.mnr_oxd(cal.alb,Si),cal.mnr_oxd(cal.alb,Mg),'s',CL{[1,4]},MS{1},12,LW{:});
set(gca,TL{:}); xlabel('SiO$_2$ [wt\%]',TX{:}); ylabel('MgO [wt\%]',TX{:});

subplot(2,3,5);
plot(exp.oxd(hasolv,olv,Si),exp.oxd(hasolv,olv,Ca),'d',CL{[1,8]},MS{:},LW{1},1.5); axis xy tight; box on; hold on
plot(exp.oxd(hasopx,opx,Si),exp.oxd(hasopx,opx,Ca),'^',CL{[1,8]},MS{:},LW{1},1.5);
plot(exp.oxd(hascpx,cpx,Si),exp.oxd(hascpx,cpx,Ca),'v',CL{[1,8]},MS{:},LW{1},1.5);
plot(exp.oxd(hasplg,plg,Si),exp.oxd(hasplg,plg,Ca),'s',CL{[1,8]},MS{:},LW{1},1.5);
plot(exp.oxd(hasmlt,mlt,Si),exp.oxd(hasmlt,mlt,Ca),'.',CL{[1,8]},MS{1},15,LW{1},1.5);
plot(mdl.oxd(hasolv,olv,Si),mdl.oxd(hasolv,olv,Ca),'d',CL{[1,2]},MS{:},LW{1},1.5);
plot(mdl.oxd(hasopx,opx,Si),mdl.oxd(hasopx,opx,Ca),'^',CL{[1,3]},MS{:},LW{1},1.5);
plot(mdl.oxd(hascpx,cpx,Si),mdl.oxd(hascpx,cpx,Ca),'v',CL{[1,3]},MS{:},LW{1},1.5);
plot(mdl.oxd(hasplg,plg,Si),mdl.oxd(hasplg,plg,Ca),'s',CL{[1,4]},MS{:},LW{1},1.5);
plot(mdl.oxd(hasmlt,mlt,Si),mdl.oxd(hasmlt,mlt,Ca),'.',CL{[1,7]},MS{1},15,LW{1},1.5);
plot(cal.mnr_oxd(cal.for,Si),cal.mnr_oxd(cal.for,Ca),'d',CL{[1,2]},MS{1},12,LW{:});
plot(cal.mnr_oxd(cal.fay,Si),cal.mnr_oxd(cal.fay,Ca),'d',CL{[1,2]},MS{1},12,LW{:});
plot(cal.mnr_oxd(cal.ens,Si),cal.mnr_oxd(cal.ens,Ca),'^',CL{[1,3]},MS{1},12,LW{:});
plot(cal.mnr_oxd(cal.hyp,Si),cal.mnr_oxd(cal.hyp,Ca),'^',CL{[1,3]},MS{1},12,LW{:});
plot(cal.mnr_oxd(cal.aug,Si),cal.mnr_oxd(cal.aug,Ca),'v',CL{[1,3]},MS{1},12,LW{:});
plot(cal.mnr_oxd(cal.pig,Si),cal.mnr_oxd(cal.pig,Ca),'v',CL{[1,3]},MS{1},12,LW{:});
plot(cal.mnr_oxd(cal.ant,Si),cal.mnr_oxd(cal.ant,Ca),'s',CL{[1,4]},MS{1},12,LW{:});
plot(cal.mnr_oxd(cal.alb,Si),cal.mnr_oxd(cal.alb,Ca),'s',CL{[1,4]},MS{1},12,LW{:});
set(gca,TL{:}); xlabel('SiO$_2$ [wt\%]',TX{:}); ylabel('CaO [wt\%]',TX{:});

subplot(2,3,6);
plot(exp.oxd(hasolv,olv,Si),exp.oxd(hasolv,olv,Na),'d',CL{[1,8]},MS{:},LW{1},1.5); axis xy tight; box on; hold on
plot(exp.oxd(hasopx,opx,Si),exp.oxd(hasopx,opx,Na),'^',CL{[1,8]},MS{:},LW{1},1.5);
plot(exp.oxd(hascpx,cpx,Si),exp.oxd(hascpx,cpx,Na),'v',CL{[1,8]},MS{:},LW{1},1.5);
plot(exp.oxd(hasplg,plg,Si),exp.oxd(hasplg,plg,Na),'s',CL{[1,8]},MS{:},LW{1},1.5);
plot(exp.oxd(hasmlt,mlt,Si),exp.oxd(hasmlt,mlt,Na),'.',CL{[1,8]},MS{1},15,LW{1},1.5);
plot(mdl.oxd(hasolv,olv,Si),mdl.oxd(hasolv,olv,Na),'d',CL{[1,2]},MS{:},LW{1},1.5);
plot(mdl.oxd(hasopx,opx,Si),mdl.oxd(hasopx,opx,Na),'^',CL{[1,3]},MS{:},LW{1},1.5);
plot(mdl.oxd(hascpx,cpx,Si),mdl.oxd(hascpx,cpx,Na),'v',CL{[1,3]},MS{:},LW{1},1.5);
plot(mdl.oxd(hasplg,plg,Si),mdl.oxd(hasplg,plg,Na),'s',CL{[1,4]},MS{:},LW{1},1.5);
plot(mdl.oxd(hasmlt,mlt,Si),mdl.oxd(hasmlt,mlt,Na),'.',CL{[1,7]},MS{1},15,LW{1},1.5);
plot(cal.mnr_oxd(cal.for,Si),cal.mnr_oxd(cal.for,Na),'d',CL{[1,2]},MS{1},12,LW{:});
plot(cal.mnr_oxd(cal.fay,Si),cal.mnr_oxd(cal.fay,Na),'d',CL{[1,2]},MS{1},12,LW{:});
plot(cal.mnr_oxd(cal.ens,Si),cal.mnr_oxd(cal.ens,Na),'^',CL{[1,3]},MS{1},12,LW{:});
plot(cal.mnr_oxd(cal.hyp,Si),cal.mnr_oxd(cal.hyp,Na),'^',CL{[1,3]},MS{1},12,LW{:});
plot(cal.mnr_oxd(cal.aug,Si),cal.mnr_oxd(cal.aug,Na),'v',CL{[1,3]},MS{1},12,LW{:});
plot(cal.mnr_oxd(cal.pig,Si),cal.mnr_oxd(cal.pig,Na),'v',CL{[1,3]},MS{1},12,LW{:});
plot(cal.mnr_oxd(cal.ant,Si),cal.mnr_oxd(cal.ant,Na),'s',CL{[1,4]},MS{1},12,LW{:});
plot(cal.mnr_oxd(cal.alb,Si),cal.mnr_oxd(cal.alb,Na),'s',CL{[1,4]},MS{1},12,LW{:});
set(gca,TL{:}); xlabel('SiO$_2$ [wt\%]',TX{:}); ylabel('Na$_2$O [wt\%]',TX{:});



%% plot best fit calibration
stg = 1;  % select which stage to plot

% plot fo-fay phase diagram
var.c = zeros(100,cal.ncmp);
var.c(:,cal.dun) = linspace(0,1,100);
var.c(:,cal.pxn) = 1-var.c(:,cal.dun);
var.T = 1200.*ones(100,1);
var.P =    1.*ones(100,1);

[~,cal]  =  meltmodel(var,cal,'T');

figure();
plot(var.c(:,cal.pxn),cal.Tsol,'k-',var.c(:,cal.pxn),cal.Tliq,'k-',LW{:}); axis xy tight; box on;
set(gca,TL{:},FS{[1,2]})
ylabel('$T \ [^\circ$C]',TX{:},FS{[1,3]})
xlabel('dun-pxn [wt]',TX{:},FS{[1,3]})
cal = rmfield(cal,{'Tsol','Tliq'});

% plot opx-cpx phase diagram
var.c = zeros(100,cal.ncmp);
var.c(:,cal.pxn) = linspace(0,1,100);
var.c(:,cal.bas) = 1-var.c(:,cal.pxn);
var.T = 1200.*ones(100,1);
var.P =    1.*ones(100,1);

[~,cal]  =  meltmodel(var,cal,'T');

figure();
plot(var.c(:,cal.bas),cal.Tsol,'k-',var.c(:,cal.bas),cal.Tliq,'k-',LW{:}); axis xy tight; box on;
set(gca,TL{:},FS{[1,2]})
ylabel('$T \ [^\circ$C]',TX{:},FS{[1,3]})
xlabel('pxn-bas [wt]',TX{:},FS{[1,3]})
cal = rmfield(cal,{'Tsol','Tliq'});

% plot an-ab phase diagram
var.c = zeros(100,cal.ncmp);
var.c(:,cal.bas) = linspace(0,1,100);
var.c(:,cal.eut) = 1-var.c(:,cal.bas);
var.T = 1200.*ones(100,1);
var.P =    1.*ones(100,1);

[~,cal]  =  meltmodel(var,cal,'T');

figure();
plot(var.c(:,cal.eut),cal.Tsol,'k-',var.c(:,cal.eut),cal.Tliq,'k-',LW{:}); axis xy tight; box on;
set(gca,TL{:},FS{[1,2]})
ylabel('$T \ [^\circ$C]',TX{:},FS{[1,3]})
xlabel('bas-eut [wt]',TX{:},FS{[1,3]})
cal = rmfield(cal,{'Tsol','Tliq'});


% plot partition coefficients (T-dependent)
var.c = ones(100,1)*c0(stg,:);
var.T = linspace(1000,1800,100).';
var.P = linspace(1,1,100).';

[~,cal]  =  meltmodel(var,cal,'K');

figure(); clf;
semilogy(var.T,cal.K,LW{:}); axis xy tight; box on;
legend(cal.cmpStr{:},FS{[1,2]},TX{:})
set(gca,TL{:},FS{[1,2]})
xlabel('$T \ [^\circ$C]',TX{:},FS{[1,3]})
ylabel('log$_{10} \ K$',TX{:},FS{[1,3]})


% plot partition coefficients (P-dependent)
var.c = ones(100,1)*c0(stg,:);
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
cal = rmfield(cal,{'Tsol','Tliq'});


% plot melt fraction, phase compositions (T-dependent)
var.c = ones(100,1)*c0(stg,:);
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
plot(var.T,phs.cl*cal.cmp_oxd,LW{:}); axis xy tight; box on;
legend(cal.oxdStr{:},FS{[1,2]},TX{:})
set(gca,TL{:},FS{[1,2]})
xlabel('$T \ [^\circ$C]',TX{:},FS{[1,3]})
ylabel('melt composition [wt\%]',TX{:},FS{[1,3]})
