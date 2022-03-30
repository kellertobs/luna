function [mdl] = runmodel (model, cal0, oxds, Temp, Pres, stages, hasolv, haspxn, hasplg)

% run melt model and cast outputs into the same form as the data

olv=1; pxn=2; plg=3; mlt=4; blk=5;
TINY = 1e-16;

cal = model2cal(cal0, model);

% extract experimental data and compute model outcome for each stage
for stg = stages

    % get computed phase proportions, compositions
    c0       =  max(0,min(1,cal.oxds.'\squeeze(oxds(stg,blk,:))));  % find end-member proportions of first step
    var.c    =  c0.'./sum(c0);
    var.T    =  Temp(stg);
    var.P    =  Pres(stg);
    [phs,~]  =  meltmodel(var,cal,'E');
    
    % record whether the newton solver in meltmodel was successful
    mdl.flag = phs.flag;
    if mdl.flag==0; return; end % newton solver not successful
    
    % record model results for all stages
    mdl.f(stg   ) = phs.f;
    mdl.cs(stg,:) = phs.cs;
    mdl.cl(stg,:) = phs.cl;
    mdl.c0(stg,:) = var.c;

    % bring model result into form of experimental data
    mdl.fphs(stg,olv) = (phs.cs(cal.fo ) + phs.cs(cal.fay)) .* (1-phs.f) .* hasolv(stg);
    mdl.fphs(stg,pxn) = (phs.cs(cal.opx) + phs.cs(cal.cpx)) .* (1-phs.f) .* haspxn(stg);
    mdl.fphs(stg,plg) = (phs.cs(cal.an ) + phs.cs(cal.ab )) .* (1-phs.f) .* hasplg(stg);
    mdl.fphs(stg,mlt) = phs.f;

    mdl.oxds(stg,olv,:) = (phs.cs(cal.fo )*cal.oxds(cal.fo ,:) + phs.cs(cal.fay)*cal.oxds(cal.fay,:)) ./ (phs.cs(cal.fo )+phs.cs(cal.fay)+TINY) .* hasolv(stg);
    mdl.oxds(stg,pxn,:) = (phs.cs(cal.opx)*cal.oxds(cal.opx,:) + phs.cs(cal.cpx)*cal.oxds(cal.cpx,:)) ./ (phs.cs(cal.opx)+phs.cs(cal.cpx)+TINY) .* haspxn(stg);
    mdl.oxds(stg,plg,:) = (phs.cs(cal.an )*cal.oxds(cal.an ,:) + phs.cs(cal.ab )*cal.oxds(cal.ab ,:)) ./ (phs.cs(cal.an )+phs.cs(cal.ab )+TINY) .* hasplg(stg);
    mdl.oxds(stg,mlt,:) = phs.cl*cal.oxds;
    
end

end