function [mdl] = runmodel (model, cal0, c0, Tmp, Prs, stages, Psl)

% run melt model and cast outputs into the same form as the data

nphs = 7; olv=1; opx=2; cpx=3; plg=4; ilm=5; qtz=6; mlt=7; blk=8; xtl=9;
TINY = 1e-16;

% cal = model2cal(cal0, model);
cal = cal0;

% extract experimental data and compute model outcome for each stage
for stg = stages

    % get computed phase proportions, compositions
    var.c    =  c0(stg,:);
    var.T    =  Tmp(stg);
    var.P    =  Prs(stg);
    [phs,~]  =  meltmodel(var,cal,'E');
    
    % record whether the newton solver in meltmodel was successful
    mdl.flag = phs.flag;
    if mdl.flag==0; return; end % newton solver not successful
    
    % record model results for all stages
    mdl.f (stg  ) = phs.f;
    mdl.cs(stg,:) = phs.cs;
    mdl.cl(stg,:) = phs.cl;
    mdl.c0(stg,:) = var.c;

    % bring model result into form of experimental data
    mdl.phs(stg,olv) = sum(phs.cs.'.*cal.cmp_mnr(:,[cal.for,cal.fay]),'all') .* (1-phs.f);
    mdl.phs(stg,opx) = sum(phs.cs.'.*cal.cmp_mnr(:,[cal.ens,cal.hyp]),'all') .* (1-phs.f);
    mdl.phs(stg,cpx) = sum(phs.cs.'.*cal.cmp_mnr(:,[cal.aug,cal.pig]),'all') .* (1-phs.f);
    mdl.phs(stg,plg) = sum(phs.cs.'.*cal.cmp_mnr(:,[cal.ant,cal.alb]),'all') .* (1-phs.f);
    mdl.phs(stg,ilm) = sum(phs.cs.'.*cal.cmp_mnr(:,[cal.ilm        ]),'all') .* (1-phs.f);
    mdl.phs(stg,qtz) = sum(phs.cs.'.*cal.cmp_mnr(:,[cal.qtz        ]),'all') .* (1-phs.f);
    mdl.phs(stg,mlt) =  phs.f;

    oxd = zeros(nphs,cal.noxd);
    wgt = zeros(nphs,cal.noxd);
    for ic = 1:cal.ncmp
        oxd(olv,:) = oxd(olv,:) + phs.cs(ic).'.*cal.cmp_mnr(ic,[cal.for,cal.fay])*cal.mnr_oxd([cal.for,cal.fay],:);
        wgt(olv,:) = wgt(olv,:) + sum(phs.cs(ic).'.*cal.cmp_mnr(ic,[cal.for,cal.fay]));
        oxd(opx,:) = oxd(opx,:) + phs.cs(ic).'.*cal.cmp_mnr(ic,[cal.ens,cal.hyp])*cal.mnr_oxd([cal.ens,cal.hyp],:);
        wgt(opx,:) = wgt(opx,:) + sum(phs.cs(ic).'.*cal.cmp_mnr(ic,[cal.ens,cal.hyp]));
        oxd(cpx,:) = oxd(cpx,:) + phs.cs(ic).'.*cal.cmp_mnr(ic,[cal.aug,cal.pig])*cal.mnr_oxd([cal.aug,cal.pig],:);
        wgt(cpx,:) = wgt(cpx,:) + sum(phs.cs(ic).'.*cal.cmp_mnr(ic,[cal.aug,cal.pig]));
        oxd(plg,:) = oxd(plg,:) + phs.cs(ic).'.*cal.cmp_mnr(ic,[cal.ant,cal.alb])*cal.mnr_oxd([cal.ant,cal.alb],:);
        wgt(plg,:) = wgt(plg,:) + sum(phs.cs(ic).'.*cal.cmp_mnr(ic,[cal.ant,cal.alb]));
        oxd(ilm,:) = oxd(ilm,:) + phs.cs(ic).'.*cal.cmp_mnr(ic,[cal.ilm        ])*cal.mnr_oxd([cal.ilm        ],:);
        wgt(ilm,:) = wgt(ilm,:) + sum(phs.cs(ic).'.*cal.cmp_mnr(ic,[cal.ilm        ]));
        oxd(qtz,:) = oxd(qtz,:) + phs.cs(ic).'.*cal.cmp_mnr(ic,[cal.qtz        ])*cal.mnr_oxd([cal.qtz        ],:);
        wgt(qtz,:) = wgt(qtz,:) + sum(phs.cs(ic).'.*cal.cmp_mnr(ic,[cal.qtz        ]));
        oxd(mlt,:) = oxd(mlt,:) + phs.cl(ic).'.*cal.cmp_mnr(ic,:)*cal.mnr_oxd(:,:);
        wgt(mlt,:) = wgt(mlt,:) + sum(phs.cl(ic).'.*cal.cmp_mnr(ic,:));
    end
    mdl.oxd(stg,:,:) = oxd./(wgt+TINY);
end

% calculate solidus and liquidus from the bulk composition
var.P    = Psl;
var.c    = mdl.c0(1,:).*ones(size(Psl));
[~,cal]  = meltmodel(var,cal,'T');
mdl.Tsol = cal.Tsol; 
mdl.Tliq = cal.Tliq;
mdl.Tm   = cal.Tm;

end