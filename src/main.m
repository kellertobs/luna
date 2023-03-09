% initialise model run
init;
    
% physical time stepping loop
while time <= tend && step <= Nt
    
    fprintf(1,'*****  step %d;  dt = %4.4e;  time = %4.4e [yr]\n\n',step,dt./yr,time./yr);
    tic;
    
    if     strcmp(TINT,'be1im') || step==1 % first step / 1st-order backward-Euler implicit scheme
        a1 = 1; a2 = 1; a3 = 0;
        b1 = 1; b2 = 0; b3 = 0;
    elseif strcmp(TINT,'bd2im') || step==2 % second step / 2nd-order 3-point backward-difference implicit scheme
        a1 = 3/2; a2 = 4/2; a3 = -1/2;
        b1 = 1;   b2 =  0;  b3 = 0;
    elseif strcmp(TINT,'cn2si')            % other steps / 2nd-order Crank-Nicolson semi-implicit scheme
        a1 = 1;   a2 = 1;   a3 = 0;
        b1 = 1/2; b2 = 1/2; b3 = 0;
    elseif strcmp(TINT,'bd2si')            % other steps / 2nd-order 3-point backward-difference semi-implicit scheme
        a1 = 3/2; a2 = 4/2; a3 = -1/2;
        b1 = 3/4; b2 = 2/4; b3 = -1/4;
    end

    % store previous solution
    Soo = So; So = S;
    Coo = Co; Co = C;
    Xoo = Xo; Xo = X;
    TEoo = TEo; TEo = TE;
    IRoo = IRo; IRo = IR;
    dSdtoo = dSdto; dSdto = dSdt;
    dCdtoo = dCdto; dCdto = dCdt;
    dXdtoo = dXdto; dXdto = dXdt;
    dTEdtoo = dTEdto; dTEdto = dTEdt;
    dIRdtoo = dIRdto; dIRdto = dIRdt;
    rhooo = rhoo; rhoo = rho;
    Div_rhoVoo = Div_rhoVo; Div_rhoVo = Div_rhoV;
    Div_Vo = Div_V;
    dto    = dt;
    
    % reset residuals and iteration count
    resnorm  = 1;
    resnorm0 = resnorm;
    iter     = 1;
    
    % non-linear iteration loop
    while resnorm/resnorm0 >= rtol && resnorm >= atol && iter <= maxit || iter < 2
        
        % solve thermo-chemical equations
        thermochem;
        
        % update non-linear parameters and auxiliary variables
        update;

        % solve fluid-mechanics equations
        fluidmech;
        
        % update geochemical evolution
        geochem;
        
        % report convergence
        report;

        iter = iter+1;
    end
        
    % record model history
    history;
    
    % print diagnostics
    fprintf(1,'\n         time to solution = %4.4f sec\n\n',toc);
    
    fprintf(1,'         min T   =  %4.1f;    mean T   = %4.1f;    max T   = %4.1f;   [degC]\n' ,min(T(:)-273.15),mean(T(:)-273.15),max(T(:)-273.15));
    fprintf(1,'         min c   =  %1.4f;    mean c   = %1.4f;    max c   = %1.4f;   [wt]\n\n'   ,min(c(:)  ),mean(c(:)  ),max(c(:)  ));
    
    fprintf(1,'         min x   =  %1.4f;    mean x   = %1.4f;    max x   = %1.4f;   [wt]\n'   ,min(x(:)  ),mean(x(:)  ),max(x(:)  ));
    fprintf(1,'         min m   =  %1.4f;    mean m   = %1.4f;    max m   = %1.4f;   [wt]\n\n' ,min(m(:)  ),mean(m(:)  ),max(m(:)  ));

    fprintf(1,'         min rho =  %4.1f;    mean rho = %4.1f;    max rho = %4.1f;   [kg/m3]\n'  ,min(rho(:)),mean(rho(:))   ,max(rho(:)));
    fprintf(1,'         min eta =  %1.2e;  mean eta = %1.2e;  max eta = %1.2e; [Pas]\n\n',min(eta(:)),geomean(eta(:)),max(eta(:)));

    fprintf(1,'         min U   = %1.4f;    mean U   = %1.4f;    max U   = %1.4f;   [m/s]\n'  ,min(U(:)  ),mean(U(:)  ),max(U(:)  ));
    fprintf(1,'         min W   = %1.4f;    mean W   = %1.4f;    max W   = %1.4f;   [m/s]\n'  ,min(-W(:) ),mean(-W(:) ),max(-W(:) ));
    fprintf(1,'         min P   = %2.4f;    mean P   = %2.4f;    max P   = %2.4f;  [kPa]\n\n',min(P(:)./1e3),mean(P(:)./1e3),max(P(:)./1e3));

    % plot results
    if ~mod(step,nop); output; end
    
    % increment time/step
    time = time+dt;
    step = step+1;
    
end

diary off
