function [mdl] = runmodel (model, cal0, c0, Tmp, Prs, stages, hasolv, haspxn, hasplg, hasspn, hasqtz, Psl)

% run melt model and cast outputs into the same form as the data

olv=1; pxn=2; plg=3; spn=4; qtz=5; mlt=6; blk=7;
TINY = 1e-16;

cal = model2cal(cal0, model);

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
    mdl.phs(stg,olv) = (                                            phs.cs(cal.for) + phs.cs(cal.fay)) .* (1-phs.f) .* hasolv(stg);
    mdl.phs(stg,pxn) = (phs.cs(cal.eut).*cal.eut_mnr(cal.mnr_px3) + phs.cs(cal.opx) + phs.cs(cal.cpx)) .* (1-phs.f) .* haspxn(stg);
    mdl.phs(stg,plg) = (phs.cs(cal.eut).*cal.eut_mnr(cal.mnr_alb) + phs.cs(cal.ant)                  ) .* (1-phs.f) .* hasplg(stg);
    mdl.phs(stg,spn) = (phs.cs(cal.eut).*cal.eut_mnr(cal.mnr_spn)                                    ) .* (1-phs.f) .* hasspn(stg);
    mdl.phs(stg,qtz) = (phs.cs(cal.eut).*cal.eut_mnr(cal.mnr_qtz)                                    ) .* (1-phs.f) .* hasqtz(stg);
    mdl.phs(stg,mlt) =  phs.f;

    mdl.oxd(stg,olv,:) = (                                                                       phs.cs(cal.for)*cal.oxd(cal.for,:) + phs.cs(cal.fay)*cal.oxd(cal.fay,:)) ./ (                                          phs.cs(cal.for)+phs.cs(cal.fay)+TINY) .* hasolv(stg);
    mdl.oxd(stg,pxn,:) = (phs.cs(cal.eut).*cal.eut_mnr(cal.mnr_px3)*cal.mnr_oxd(cal.mnr_px3,:) + phs.cs(cal.opx)*cal.oxd(cal.opx,:) + phs.cs(cal.cpx)*cal.oxd(cal.cpx,:)) ./ (phs.cs(cal.eut).*cal.eut_mnr(cal.mnr_px3)+phs.cs(cal.opx)+phs.cs(cal.cpx)+TINY) .* haspxn(stg);
    mdl.oxd(stg,plg,:) = (phs.cs(cal.eut).*cal.eut_mnr(cal.mnr_alb)*cal.mnr_oxd(cal.mnr_alb,:) + phs.cs(cal.ant)*cal.oxd(cal.ant,:)                                     ) ./ (phs.cs(cal.eut).*cal.eut_mnr(cal.mnr_alb)+phs.cs(cal.ant)                +TINY) .* hasplg(stg);
    mdl.oxd(stg,spn,:) = (phs.cs(cal.eut).*cal.eut_mnr(cal.mnr_spn)*cal.mnr_oxd(cal.mnr_spn,:)                                                                          ) ./ (phs.cs(cal.eut).*cal.eut_mnr(cal.mnr_spn)                                +TINY) .* hasspn(stg);
    mdl.oxd(stg,qtz,:) = (phs.cs(cal.eut).*cal.eut_mnr(cal.mnr_qtz)*cal.mnr_oxd(cal.mnr_qtz,:)                                                                          ) ./ (phs.cs(cal.eut).*cal.eut_mnr(cal.mnr_qtz)                                +TINY) .* hasqtz(stg);
    mdl.oxd(stg,mlt,:) =  phs.cl*cal.oxd;
    
end

% calculate solidus and liquidus from the bulk composition
var.P    = Psl;
var.c    = mdl.c0(1,:).*ones(size(Psl));
[~,cal]  = meltmodel(var,cal,'T');
mdl.Tsol = cal.Tsol; 
mdl.Tliq = cal.Tliq;
mdl.Tm   = cal.Tm;

end