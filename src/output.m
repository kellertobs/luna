
if plot_op
    % prepare for plotting
    TX = {'Interpreter','Latex'}; FS = {'FontSize',12};
    TL = {'TickLabelInterpreter','Latex'}; TS = {'FontSize',10};
    UN = {'Units','Centimeters'};
    CL = {'Color',[0.0 0.0 0.0],[0.80 0.15 0.10],[0.10 0.15 0.65],[0.45 0.60 0.95]};
    LW = {'LineWidth',2};
    
    if Nx <= 10 && Nz <= 10  % create 0D plots
            
        fh1 = figure(1); clf;
        subplot(4,1,1)
        plot(hist.time/hr,hist.mu (:,2)*100.*(hist.mu (:,2)>1e-9),CL{[1,3]},LW{:}); axis xy tight; box on; hold on;
        plot(hist.time/hr,hist.chi(:,2)*100.*(hist.chi(:,2)>1e-9),CL{[1,4]},LW{:});
        title(['$\mu$, $\chi$ [vol\%]'],TX{:},FS{:}); set(gca,TL{:},TS{:});
        subplot(4,1,2)
        for i =1:cal.nc; plot(hist.time/hr,hist.c (i,:,2)*100,LW{:});  axis xy tight; box on; hold on; end
        title('$\bar{c}$ [wt\%]',TX{:},FS{:}); set(gca,TL{:},TS{:});
        subplot(4,1,3)
        for i =1:cal.nc; plot(hist.time/hr,hist.cx(i,:,2)*100,LW{:});  axis xy tight; box on; hold on; end
        title('$c^x$ [wt\%]',TX{:},FS{:}); set(gca,TL{:},TS{:});
        legend(cal.CompStr,TX{:},FS{:},'location','best');
        subplot(4,1,4)
        for i =1:cal.nc; plot(hist.time/hr,hist.cm(i,:,2)*100,LW{:});  axis xy tight; box on; hold on; end
        title('$c^m$ [wt\%]',TX{:},FS{:}); set(gca,TL{:},TS{:});
        xlabel('Time [hr]',TX{:},FS{:});
        
        fh2 = figure(2); clf;
        subplot(4,1,1)
        plot(hist.time/hr,hist.T(:,2)-273.15,CL{[1,2]},LW{:}); axis xy tight; box on;
        title('$T [^\circ$C]',TX{:},FS{:}); set(gca,TL{:},TS{:});
        subplot(4,1,2)
        plot(hist.time/hr,hist.rhom(:,2),'-',CL{[1,3]},LW{:}); axis xy tight; box on; hold on;
        plot(hist.time/hr,hist.rhox(:,2),'-',CL{[1,4]},LW{:});
        plot(hist.time/hr,hist.rho (:,2),'-',CL{[1,2]},LW{:});
        title('$\bar{\rho}$ [kg/m$^3$]',TX{:},FS{:}); set(gca,TL{:},TS{:});
        subplot(4,1,3)
        semilogy(hist.time/hr,hist.etam(:,2),'r-',LW{:}); axis xy tight; box on; hold on;
        semilogy(hist.time/hr,hist.eta (:,2),'k-',LW{:});
        title('$\bar{\eta}$ [Pas]',TX{:},FS{:}); set(gca,TL{:},TS{:});
        subplot(4,1,4)
        plot(hist.time/hr,hist.Gx(:,2)./hist.rho(:,2)*hr*100.*(hist.chi(:,2)>1e-9),CL{[1,4]},LW{:}); axis xy tight; box on; hold on;
        title('$\Gamma_x/\bar{\rho}$ [wt\%/hr]',TX{:},FS{:}); set(gca,TL{:},TS{:});
        xlabel('Time [hr]',TX{:},FS{:});
        
        fh3 = figure(3); clf;
        subplot(3,1,1)
        for i =1:cal.nc; plot(hist.time/hr,hist.oxd (i,:,2),LW{:});  axis xy tight; box on; hold on; end
        title('$\bar{c}_\mathrm{oxd}$ [wt\%]',TX{:},FS{:}); set(gca,TL{:},TS{:});
        subplot(3,1,2)
        for i =1:cal.nc; plot(hist.time/hr,hist.oxdx(i,:,2),LW{:});  axis xy tight; box on; hold on; end
        title('$c^x_\mathrm{oxd}$ [wt\%]',TX{:},FS{:}); set(gca,TL{:},TS{:});
        legend(cal.oxdStr,TX{:},FS{:},'location','best');
        subplot(3,1,3)
        for i =1:cal.nc; plot(hist.time/hr,hist.oxdm(i,:,2),LW{:});  axis xy tight; box on; hold on; end
        title('$c^m_\mathrm{oxd}$ [wt\%]',TX{:},FS{:}); set(gca,TL{:},TS{:});
        xlabel('Time [hr]',TX{:},FS{:});
        
    elseif Nx <= 10  % create 1D plots
        
        fh1 = figure(1); clf;
        subplot(1,5,1)
        plot(T(2:end-1,2:end-1)-273.15,Z(2:end-1).'./1e3,CL{[1,2]},LW{:}); axis ij tight; box on;
        title('$T [^\circ$C]',TX{:},FS{:}); ylabel('Depth [km]',TX{:},FS{:}); set(gca,TL{:},TS{:});
        subplot(1,5,2)
        plot(mu (2:end-1,2:end-1)*100.*(mu (2:end-1,2:end-1)>1e-9),Z(2:end-1).'./1e3,CL{[1,3]},LW{:}); axis ij tight; box on; hold on;
        plot(chi(2:end-1,2:end-1)*100.*(chi(2:end-1,2:end-1)>1e-9),Z(2:end-1).'./1e3,CL{[1,4]},LW{:});
        title(['$\mu$, $\chi$ [vol\%]'],TX{:},FS{:}); set(gca,TL{:},TS{:});
        subplot(1,5,3)
        for i =1:cal.nc; plot(c(i,2:end-1,2:end-1)*100,Z(2:end-1).'./1e3,LW{:}); axis ij tight; box on; hold on; end
        title('$\bar{c}$ [wt\%]',TX{:},FS{:}); set(gca,TL{:},TS{:});
        subplot(1,5,4)
        for i =1:cal.nc; plot(cx(i,2:end-1,2:end-1)*100,Z(2:end-1).'./1e3,LW{:}); axis ij tight; box on; hold on; end
        title('$c^x$ [wt\%]',TX{:},FS{:}); set(gca,TL{:},TS{:});
        legend(cal.CompStr,TX{:},FS{:},'location','best');
        subplot(1,5,5)
        for i =1:cal.nc; plot(cm(i,2:end-1,2:end-1)*100,Z(2:end-1).'./1e3,LW{:}); axis ij tight; box on; hold on; end
        title('$c^m$ [wt\%]',TX{:},FS{:}); set(gca,TL{:},TS{:});

        fh2 = figure(2); clf;
        subplot(1,5,1)
        plot(rhox(2:end-1,2:end-1),Z(2:end-1).'./1e3,CL{[1,4]},LW{:}); axis ij tight; box on; hold on;
        plot(rhom(2:end-1,2:end-1),Z(2:end-1).'./1e3,CL{[1,3]},LW{:});
        plot(rho (2:end-1,2:end-1),Z(2:end-1).'./1e3,CL{[1,2]},LW{:});
        title('$\bar{\rho}$ [kg/m$^3$]',TX{:},FS{:}); ylabel('Depth [km]',TX{:},FS{:}); set(gca,TL{:},TS{:});
        subplot(1,5,2)
        semilogx(min(eta(2:end-1,2:end-1),etam(2:end-1,2:end-1)),Z(2:end-1).'./1e3,CL{[1,3]},LW{:}); axis ij tight; box on; hold on;
        semilogx(eta(2:end-1,2:end-1),Z(2:end-1).'./1e3,CL{[1,2]},LW{:});
        title('$\bar{\eta}$ [Pas]',TX{:},FS{:}); set(gca,TL{:},TS{:});
        subplot(1,5,3)
        plot(Gx(2:end-1,2:end-1)./rho(2:end-1,2:end-1)*hr*100.*(chi(2:end-1,2:end-1)>1e-9),Z(2:end-1).'./1e3,CL{[1,4]},LW{:}); axis ij tight; box on; hold on;
        title('$\Gamma_x/\bar{\rho}$ [wt/hr]',TX{:},FS{:}); set(gca,TL{:},TS{:});
        subplot(1,5,4)
        plot(-W(:,2:end-1)*hr,Zfc.',CL{[1,2]},LW{:}); axis ij tight; box on; hold on;
        plot(-(chi(1:end-1,2:end-1)+chi(2:end,2:end-1))/2.*wx(:,2:end-1)*hr,Zfc.',CL{[1,4]},LW{:});
        plot(-(mu (1:end-1,2:end-1)+mu (2:end,2:end-1))/2.*wm(:,2:end-1)*hr,Zfc.',CL{[1,3]},LW{:});
        title('$W$, $w_\Delta^x$ [m/hr]',TX{:},FS{:}); ylabel('Depth [km]',TX{:},FS{:}); set(gca,TL{:},TS{:});
        subplot(1,5,5)
        plot(P(2:end-1,2:end-1)/1e3,Z(2:end-1).'./1e3,CL{[1,2]},LW{:}); axis ij tight; box on;
        title('$P$ [kPa]',TX{:},FS{:}); set(gca,TL{:},TS{:});
        
        fh3 = figure(3); clf;
        subplot(1,3,1)
        for i =1:cal.nc; plot(oxd (i,2:end-1,2:end-1),Z(2:end-1).'./1e3,LW{:}); axis ij tight; box on; hold on; end
        title('$\bar{c}_\mathrm{oxd}$ [wt\%]',TX{:},FS{:}); set(gca,TL{:},TS{:});
        subplot(1,3,2)
        for i =1:cal.nc; plot(oxdx(i,2:end-1,2:end-1),Z(2:end-1).'./1e3,LW{:}); axis ij tight; box on; hold on; end
        title('$c^x_\mathrm{oxd}$ [wt\%]',TX{:},FS{:}); set(gca,TL{:},TS{:});
        legend(cal.oxdStr,TX{:},FS{:},'location','best');
        subplot(1,3,3)
        for i =1:cal.nc; plot(oxdm(i,2:end-1,2:end-1),Z(2:end-1).'./1e3,LW{:}); axis ij tight; box on; hold on; end
        title('$c^m_\mathrm{oxd}$ [wt\%]',TX{:},FS{:}); set(gca,TL{:},TS{:});
        
        fh4 = figure(4); clf;
        subplot(1,5,1)
        semilogx(max(1e-3,min(1e3,itx(2:end-1,2:end-1)))*100,Z(2:end-1).'./1e3,CL{[1,4]},LW{:}); axis ij tight; box on; hold on;
        semilogx(max(1e-3,min(1e3,itm(2:end-1,2:end-1)))*100,Z(2:end-1).'./1e3,CL{[1,3]},LW{:});
        semilogx(max(1e-3,min(1e3,it (2:end-1,2:end-1)))*100,Z(2:end-1).'./1e3,CL{[1,2]},LW{:});
        title('incomp. trace',TX{:},FS{:}); ylabel('Depth [km]',TX{:},FS{:}); set(gca,TL{:},TS{:});
        subplot(1,5,2)
        semilogx(max(1e-3,min(1e3,ctx(2:end-1,2:end-1)))*100,Z(2:end-1).'./1e3,CL{[1,4]},LW{:}); axis ij tight; box on; hold on;
        semilogx(max(1e-3,min(1e3,ctm(2:end-1,2:end-1)))*100,Z(2:end-1).'./1e3,CL{[1,3]},LW{:});
        semilogx(max(1e-3,min(1e3,ct (2:end-1,2:end-1)))*100,Z(2:end-1).'./1e3,CL{[1,2]},LW{:});
        title('comp. trace',TX{:},FS{:}); set(gca,TL{:},TS{:});
        subplot(1,5,3)
        plot(si(2:end-1,2:end-1),Z(2:end-1).'./1e3,CL{[1,2]},LW{:}); axis ij tight; box on;
        title('stable isotope',TX{:},FS{:}); set(gca,TL{:},TS{:});
        subplot(1,5,4)
        semilogx(max(1e-3,min(1e3,ripx(2:end-1,2:end-1)))*100,Z(2:end-1).'./1e3,CL{[1,4]},LW{:}); axis ij tight; box on; hold on;
        semilogx(max(1e-3,min(1e3,ripm(2:end-1,2:end-1)))*100,Z(2:end-1).'./1e3,CL{[1,3]},LW{:});
        semilogx(max(1e-3,min(1e3,rip (2:end-1,2:end-1)))*100,Z(2:end-1).'./1e3,CL{[1,2]},LW{:});
        title(['radiogenic parent'],TX{:},FS{:}); set(gca,TL{:},TS{:});
        subplot(1,5,5)
        semilogx(max(1e-3,min(1e3,ridx(2:end-1,2:end-1)))*100,Z(2:end-1).'./1e3,CL{[1,4]},LW{:}); axis ij tight; box on; hold on;
        semilogx(max(1e-3,min(1e3,ridm(2:end-1,2:end-1)))*100,Z(2:end-1).'./1e3,CL{[1,3]},LW{:});
        semilogx(max(1e-3,min(1e3,rid (2:end-1,2:end-1)))*100,Z(2:end-1).'./1e3,CL{[1,2]},LW{:});
        title(['radiogenic daughter'],TX{:},FS{:}); set(gca,TL{:},TS{:});
        
    else % create 2D plots
        
    % set axis and border dimensions
    axh = 6.00; axw = axh*L/D;
    ahs = 0.30; avs = 0.0;
    axb = 0.40; axt = 0.6;
    axl = 1.40; axr = 0.4;
    
    % initialize figures and axes
    fh1 = figure(1); clf; colormap(ocean);
    fh = axb + 2*axh + 1*avs + axt;
    fw = axl + 2*axw + 1*ahs + axr;
    set(fh1,UN{:},'Position',[1 1 fw fh]);
    set(fh1,'PaperUnits','Centimeters','PaperPosition',[0 0 fw fh],'PaperSize',[fw fh]);
    set(fh1,'Color','w','InvertHardcopy','off');
    set(fh1,'Resize','off');
    ax(11) = axes(UN{:},'position',[axl+0*axw+0*ahs axb+1*axh+1*avs axw axh]);
    ax(12) = axes(UN{:},'position',[axl+1*axw+1*ahs axb+1*axh+1*avs axw axh]);
    ax(13) = axes(UN{:},'position',[axl+0*axw+0*ahs axb+0*axh+0*avs axw axh]);
    ax(14) = axes(UN{:},'position',[axl+1*axw+1*ahs axb+0*axh+0*avs axw axh]);

    fh2 = figure(2); clf; colormap(ocean);
    fh = axb + 2*axh + 1*avs + axt;
    fw = axl + 3*axw + 2*ahs + axr;
    set(fh2,UN{:},'Position',[3 3 fw fh]);
    set(fh2,'PaperUnits','Centimeters','PaperPosition',[0 0 fw fh],'PaperSize',[fw fh]);
    set(fh2,'Color','w','InvertHardcopy','off');
    set(fh2,'Resize','off');
    ax(21) = axes(UN{:},'position',[axl+0*axw+0*ahs axb+1*axh+1*avs axw axh]);
    ax(22) = axes(UN{:},'position',[axl+1*axw+1*ahs axb+1*axh+1*avs axw axh]);
    ax(23) = axes(UN{:},'position',[axl+2*axw+2*ahs axb+1*axh+1*avs axw axh]);
    ax(24) = axes(UN{:},'position',[axl+0*axw+0*ahs axb+0*axh+0*avs axw axh]);
    ax(25) = axes(UN{:},'position',[axl+1*axw+1*ahs axb+0*axh+0*avs axw axh]);
    ax(26) = axes(UN{:},'position',[axl+2*axw+2*ahs axb+0*axh+0*avs axw axh]);
    
    fh3 = figure(3); clf; colormap(ocean);
    fh = axb + 2*axh + 1*avs + axt;
    fw = axl + 3*axw + 2*ahs + axr;
    set(fh3,UN{:},'Position',[5 5 fw fh]);
    set(fh3,'PaperUnits','Centimeters','PaperPosition',[0 0 fw fh],'PaperSize',[fw fh]);
    set(fh3,'Color','w','InvertHardcopy','off');
    set(fh3,'Resize','off');
    ax(31) = axes(UN{:},'position',[axl+0*axw+0*ahs axb+1*axh+1*avs axw axh]);
    ax(32) = axes(UN{:},'position',[axl+1*axw+1*ahs axb+1*axh+1*avs axw axh]);
    ax(33) = axes(UN{:},'position',[axl+2*axw+2*ahs axb+1*axh+1*avs axw axh]);
    ax(34) = axes(UN{:},'position',[axl+0*axw+0*ahs axb+0*axh+0*avs axw axh]);
    ax(35) = axes(UN{:},'position',[axl+1*axw+1*ahs axb+0*axh+0*avs axw axh]);
    ax(36) = axes(UN{:},'position',[axl+2*axw+2*ahs axb+0*axh+0*avs axw axh]);
    
    fh4 = figure(4); clf; colormap(ocean);
    fh = axb + 2*axh + 1*avs + axt;
    fw = axl + 3*axw + 2*ahs + axr;
    set(fh4,UN{:},'Position',[7 7 fw fh]);
    set(fh4,'PaperUnits','Centimeters','PaperPosition',[0 0 fw fh],'PaperSize',[fw fh]);
    set(fh4,'Color','w','InvertHardcopy','off');
    set(fh4,'Resize','off');
    ax(41) = axes(UN{:},'position',[axl+0*axw+0*ahs axb+1*axh+1*avs axw axh]);
    ax(42) = axes(UN{:},'position',[axl+1*axw+1*ahs axb+1*axh+1*avs axw axh]);
    ax(43) = axes(UN{:},'position',[axl+2*axw+2*ahs axb+1*axh+1*avs axw axh]);
    ax(44) = axes(UN{:},'position',[axl+0*axw+0*ahs axb+0*axh+0*avs axw axh]);
    ax(45) = axes(UN{:},'position',[axl+1*axw+1*ahs axb+0*axh+0*avs axw axh]);
    ax(46) = axes(UN{:},'position',[axl+2*axw+2*ahs axb+0*axh+0*avs axw axh]);
    
    fh5 = figure(5); clf; colormap(ocean);
    fh = axb + 2*axh + 1*avs + axt;
    fw = axl + 3*axw + 2*ahs + axr;
    set(fh5,UN{:},'Position',[9 9 fw fh]);
    set(fh5,'PaperUnits','Centimeters','PaperPosition',[0 0 fw fh],'PaperSize',[fw fh]);
    set(fh5,'Color','w','InvertHardcopy','off');
    set(fh5,'Resize','off');
    ax(51) = axes(UN{:},'position',[axl+0*axw+0*ahs axb+1*axh+1*avs axw axh]);
    ax(52) = axes(UN{:},'position',[axl+1*axw+1*ahs axb+1*axh+1*avs axw axh]);
    ax(53) = axes(UN{:},'position',[axl+2*axw+2*ahs axb+1*axh+1*avs axw axh]);
    ax(54) = axes(UN{:},'position',[axl+0*axw+0*ahs axb+0*axh+0*avs axw axh]);
    ax(55) = axes(UN{:},'position',[axl+1*axw+1*ahs axb+0*axh+0*avs axw axh]);
    ax(56) = axes(UN{:},'position',[axl+2*axw+2*ahs axb+0*axh+0*avs axw axh]);
    
    if plot_cv
        fh6 = figure(6); clf; colormap(ocean);
        fh = axb + 1*axh + 0*avs + axt;
        fw = axl + 3*axw + 2*ahs + axr;
        set(fh6,UN{:},'Position',[11 11 fw fh]);
        set(fh6,'PaperUnits','Centimeters','PaperPosition',[0 0 fw fh],'PaperSize',[fw fh]);
        set(fh6,'Color','w','InvertHardcopy','off');
        set(fh6,'Resize','off');
        ax(61) = axes(UN{:},'position',[axl+0*axw+0*ahs axb+0*axh+0*avs axw axh]);
        ax(62) = axes(UN{:},'position',[axl+1*axw+1*ahs axb+0*axh+0*avs axw axh]);
        ax(63) = axes(UN{:},'position',[axl+2*axw+2*ahs axb+0*axh+0*avs axw axh]);
    end
    
    % plot velocity-pressure solution in Fig. 1
    figure(1);
    sgtitle(['time = ',num2str(time/hr,3),' [hr]'],TX{:},FS{:},'Color','k');
    axes(ax(11));
    imagesc(X(2:end-1)./1e3,Z(2:end-1)./1e3,-W(:      ,2:end-1).*hr); axis ij equal tight; box on; cb = colorbar;
    set(cb,TL{:},TS{:}); set(gca,TL{:},TS{:}); title(['$W$ [m/hr]'],TX{:},FS{:}); ylabel('Depth [km]',TX{:},FS{:}); ylabel('Depth [km]',TX{:},FS{:});  set(gca,'XTickLabel',[]);
    axes(ax(12));
    imagesc(X(2:end-1)./1e3,Z(2:end-1)./1e3, U(2:end-1,:      ).*hr); axis ij equal tight; box on; cb = colorbar;
    set(cb,TL{:},TS{:}); set(gca,TL{:},TS{:}); title(['$U$ [m/hr]'],TX{:},FS{:}); set(gca,'XTickLabel',[],'YTickLabel',[]);
    axes(ax(13));
    imagesc(X(2:end-1)./1e3,Z(2:end-1)./1e3, P(2:end-1,2:end-1)./1e3); axis ij equal tight; box on; cb = colorbar;
    set(cb,TL{:},TS{:}); set(gca,TL{:},TS{:}); title(['$P$ [kPa]'],TX{:},FS{:}); ylabel('Depth [km]',TX{:},FS{:}); xlabel('Width [km]',TX{:},FS{:});
    axes(ax(14));
    imagesc(X(2:end-1)./1e3,Z(2:end-1)./1e3,T(2:end-1,2:end-1)-273.15); axis ij equal tight; box on; cb = colorbar;
    set(cb,TL{:},TS{:}); set(gca,TL{:},TS{:}); title(['$T [^\circ$C]'],TX{:},FS{:}); xlabel('Width [km]',TX{:},FS{:}); set(gca,'YTickLabel',[]);
    
    % plot pseudo-component composition in Fig. 2
    figure(2);
    sgtitle(['time = ',num2str(time/hr,3),' [hr]'],TX{:},FS{:},'Color','k');
    axes(ax(21));
    imagesc(X(2:end-1)./1e3,Z(2:end-1)./1e3,squeeze(c(1,2:end-1,2:end-1))*100); axis ij equal tight; box on; cb = colorbar;
    set(cb,TL{:},TS{:}); set(gca,TL{:},TS{:}); title([cal.CompStr{1},' [wt\%]'],TX{:},FS{:}); set(gca,'XTickLabel',[]); ylabel('Depth [km]',TX{:},FS{:});
    axes(ax(22));
    imagesc(X(2:end-1)./1e3,Z(2:end-1)./1e3,squeeze(c(2,2:end-1,2:end-1))*100); axis ij equal tight; box on; cb = colorbar;
    set(cb,TL{:},TS{:}); set(gca,TL{:},TS{:}); title([cal.CompStr{2},' [wt\%]'],TX{:},FS{:}); set(gca,'XTickLabel',[],'YTickLabel',[]);
    axes(ax(23));
    imagesc(X(2:end-1)./1e3,Z(2:end-1)./1e3,squeeze(c(3,2:end-1,2:end-1))*100); axis ij equal tight; box on; cb = colorbar;
    set(cb,TL{:},TS{:}); set(gca,TL{:},TS{:}); title([cal.CompStr{3},' [wt\%]'],TX{:},FS{:}); set(gca,'XTickLabel',[],'YTickLabel',[]);
    axes(ax(24));
    imagesc(X(2:end-1)./1e3,Z(2:end-1)./1e3,squeeze(c(4,2:end-1,2:end-1))*100); axis ij equal tight; box on; cb = colorbar;
    set(cb,TL{:},TS{:}); set(gca,TL{:},TS{:}); title([cal.CompStr{4},' [wt\%]'],TX{:},FS{:}); ylabel('Depth [km]',TX{:},FS{:});
    axes(ax(25));
    imagesc(X(2:end-1)./1e3,Z(2:end-1)./1e3,squeeze(c(5,2:end-1,2:end-1))*100); axis ij equal tight; box on; cb = colorbar;
    set(cb,TL{:},TS{:}); set(gca,TL{:},TS{:}); title([cal.CompStr{5},' [wt\%]'],TX{:},FS{:}); set(gca,'YTickLabel',[]); xlabel('Width [km]',TX{:},FS{:});
    axes(ax(26));
    imagesc(X(2:end-1)./1e3,Z(2:end-1)./1e3,squeeze(c(6,2:end-1,2:end-1))*100); axis ij equal tight; box on; cb = colorbar;
    set(cb,TL{:},TS{:}); set(gca,TL{:},TS{:}); title([cal.CompStr{6},' [wt\%]'],TX{:},FS{:}); set(gca,'YTickLabel',[]);

    % plot major oxide composition in Fig. 2
    figure(3);
    sgtitle(['time = ',num2str(time/hr,3),' [hr]'],TX{:},FS{:},'Color','k');
    axes(ax(31));
    imagesc(X(2:end-1)./1e3,Z(2:end-1)./1e3,squeeze(oxd(1,2:end-1,2:end-1))); axis ij equal tight; box on; cb = colorbar;
    set(cb,TL{:},TS{:}); set(gca,TL{:},TS{:}); title([cal.oxdStr{1},' [wt\%]'],TX{:},FS{:}); set(gca,'XTickLabel',[]); ylabel('Depth [km]',TX{:},FS{:});
    axes(ax(32));
    imagesc(X(2:end-1)./1e3,Z(2:end-1)./1e3,squeeze(oxd(2,2:end-1,2:end-1))); axis ij equal tight; box on; cb = colorbar;
    set(cb,TL{:},TS{:}); set(gca,TL{:},TS{:}); title([cal.oxdStr{2},' [wt\%]'],TX{:},FS{:}); set(gca,'XTickLabel',[],'YTickLabel',[]);
    axes(ax(33));
    imagesc(X(2:end-1)./1e3,Z(2:end-1)./1e3,squeeze(oxd(3,2:end-1,2:end-1))); axis ij equal tight; box on; cb = colorbar;
    set(cb,TL{:},TS{:}); set(gca,TL{:},TS{:}); title([cal.oxdStr{3},' [wt\%]'],TX{:},FS{:}); set(gca,'XTickLabel',[],'YTickLabel',[]);
    axes(ax(34));
    imagesc(X(2:end-1)./1e3,Z(2:end-1)./1e3,squeeze(oxd(4,2:end-1,2:end-1))); axis ij equal tight; box on; cb = colorbar;
    set(cb,TL{:},TS{:}); set(gca,TL{:},TS{:}); title([cal.oxdStr{4},' [wt\%]'],TX{:},FS{:}); ylabel('Depth [km]',TX{:},FS{:});
    axes(ax(35));
    imagesc(X(2:end-1)./1e3,Z(2:end-1)./1e3,squeeze(oxd(5,2:end-1,2:end-1))); axis ij equal tight; box on; cb = colorbar;
    set(cb,TL{:},TS{:}); set(gca,TL{:},TS{:}); title([cal.oxdStr{5},' [wt\%]'],TX{:},FS{:}); set(gca,'YTickLabel',[]); xlabel('Width [km]',TX{:},FS{:});
    axes(ax(36));
    imagesc(X(2:end-1)./1e3,Z(2:end-1)./1e3,squeeze(oxd(6,2:end-1,2:end-1))); axis ij equal tight; box on; cb = colorbar;
    set(cb,TL{:},TS{:}); set(gca,TL{:},TS{:}); title([cal.oxdStr{6},' [wt\%]'],TX{:},FS{:}); set(gca,'YTickLabel',[]);
    
    % plot phase fractions and reaction rates in Fig. 3
    figure(4);
    sgtitle(['time = ',num2str(time/hr,3),' [hr]'],TX{:},FS{:},'Color','k');
    axes(ax(41));
    imagesc(X(2:end-1)./1e3,Z(2:end-1)./1e3,chi(2:end-1,2:end-1).*100.*(chi(2:end-1,2:end-1)>1e-9) ); axis ij equal tight; box on; cb = colorbar;
    set(cb,TL{:},TS{:}); set(gca,TL{:},TS{:}); title(['$\chi$ [vol\%]'],TX{:},FS{:}); set(gca,'XTickLabel',[]); ylabel('Depth [km]',TX{:},FS{:});
    axes(ax(42));
    imagesc(X(2:end-1)./1e3,Z(2:end-1)./1e3,Gx(2:end-1,2:end-1)./rho(2:end-1,2:end-1)*hr*100.*(chi(2:end-1,2:end-1)>1e-9)); axis ij equal tight; box on; cb = colorbar;
    set(cb,TL{:},TS{:}); set(gca,TL{:},TS{:}); title(['$\Gamma_x/\bar{\rho}$ [wt\%/hr]'],TX{:},FS{:}); set(gca,'XTickLabel',[],'YTickLabel',[]);
    axes(ax(43));
    imagesc(X(2:end-1)./1e3,Z(2:end-1)./1e3,      rho(2:end-1,2:end-1) ); axis ij equal tight; box on; cb = colorbar;
    set(cb,TL{:},TS{:}); set(gca,TL{:},TS{:}); title(['$\bar{\rho}$ [kg/m$^3$]'],TX{:},FS{:}); set(gca,'XTickLabel',[],'YTickLabel',[]);
    axes(ax(44));
    imagesc(X(2:end-1)./1e3,Z(2:end-1)./1e3,-(chi(1:end-1,2:end-1)+chi(2:end,2:end-1))/2.*wx(:,2:end-1).*hr); axis ij equal tight; box on; cb = colorbar;
    set(cb,TL{:},TS{:}); set(gca,TL{:},TS{:}); title(['$w_\Delta^x$ [m/hr]'],TX{:},FS{:}); ylabel('Depth [km]',TX{:},FS{:}); ylabel('Depth [km]',TX{:},FS{:});
    axes(ax(45));
    imagesc(X(2:end-1)./1e3,Z(2:end-1)./1e3,-(mu (1:end-1,2:end-1)+mu (2:end,2:end-1))/2.*wm(:,2:end-1).*hr); axis ij equal tight; box on; cb = colorbar;
    set(cb,TL{:},TS{:}); set(gca,TL{:},TS{:}); title(['$w_\Delta^m$ [m/hr]'],TX{:},FS{:}); set(gca,'YTickLabel',[]);
    axes(ax(46));
    imagesc(X(2:end-1)./1e3,Z(2:end-1)./1e3,log10(eta(2:end-1,2:end-1))); axis ij equal tight; box on; cb = colorbar;
    set(cb,TL{:},TS{:}); set(gca,TL{:},TS{:}); title(['$\bar{\eta}$ [log$_{10}$ Pas]'],TX{:},FS{:}); set(gca,'YTickLabel',[]);
    
    % plot geochemical variables in Fig. 5
    figure(5);
    sgtitle(['time = ',num2str(time/hr,3),' [hr]'],TX{:},FS{:},'Color','k');
    axes(ax(51));
    imagesc(X(2:end-1)./1e3,Z(2:end-1)./1e3,it(2:end-1,2:end-1)); axis ij equal tight; box on; cb = colorbar;
    set(cb,TL{:},TS{:}); set(gca,TL{:},TS{:}); title(['incomp. trace'],TX{:},FS{:}); set(gca,'XTickLabel',[]); ylabel('Depth [km]',TX{:},FS{:});
    axes(ax(52));
    imagesc(X(2:end-1)./1e3,Z(2:end-1)./1e3,ct(2:end-1,2:end-1)); axis ij equal tight; box on; cb = colorbar;
    set(cb,TL{:},TS{:}); set(gca,TL{:},TS{:}); title(['comp. trace'],TX{:},FS{:}); set(gca,'XTickLabel',[],'YTickLabel',[]);
    axes(ax(53));
    imagesc(X(2:end-1)./1e3,Z(2:end-1)./1e3,si(2:end-1,2:end-1)); axis ij equal tight; box on; cb = colorbar;
    set(cb,TL{:},TS{:}); set(gca,TL{:},TS{:}); title(['stable isotope'],TX{:},FS{:}); set(gca,'XTickLabel',[],'YTickLabel',[]);
    axes(ax(54));
    imagesc(X(2:end-1)./1e3,Z(2:end-1)./1e3,rip(2:end-1,2:end-1)); axis ij equal tight; box on; cb = colorbar;
    set(cb,TL{:},TS{:}); set(gca,TL{:},TS{:}); title(['radiogen. parent'],TX{:},FS{:}); ylabel('Depth [km]',TX{:},FS{:});
    axes(ax(55));
    imagesc(X(2:end-1)./1e3,Z(2:end-1)./1e3,rid(2:end-1,2:end-1)); axis ij equal tight; box on; cb = colorbar;
    set(cb,TL{:},TS{:}); set(gca,TL{:},TS{:}); title(['radiogen. daughter'],TX{:},FS{:}); set(gca,'YTickLabel',[]); xlabel('Width [km]',TX{:},FS{:});
    axes(ax(56));
    imagesc(X(2:end-1)./1e3,Z(2:end-1)./1e3,(dcy_rip(2:end-1,2:end-1)-dcy_rid(2:end-1,2:end-1))./(dcy_rip(2:end-1,2:end-1)+dcy_rid(2:end-1,2:end-1))); axis ij equal tight; box on; cb = colorbar;
    set(cb,TL{:},TS{:}); set(gca,TL{:},TS{:}); title(['radiogen. disequilibrium'],TX{:},FS{:}); set(gca,'YTickLabel',[]);
    
    if plot_cv && iter > 0
        % plot residual fields in Fig. 4
        figure(6);
        axes(ax(61));
        imagesc(X(2:end-1)./1e3,Z(2:end-1)./1e3,-res_W(:      ,2:end-1)); axis ij equal tight; box on; cb = colorbar;
        set(cb,TL{:},TS{:}); set(gca,TL{:},TS{:}); title(['res. $W$'],TX{:},FS{:}); ylabel('Depth [km]',TX{:},FS{:});
        axes(ax(62));
        imagesc(X(2:end-1)./1e3,Z(2:end-1)./1e3, res_U(2:end-1,:      )); axis ij equal tight; box on; cb = colorbar;
        set(cb,TL{:},TS{:}); set(gca,TL{:},TS{:}); title(['res. $U$'],TX{:},FS{:}); set(gca,'YTickLabel',[]); xlabel('Width [km]',TX{:},FS{:});
        axes(ax(63));
        imagesc(X(2:end-1)./1e3,Z(2:end-1)./1e3, res_P(2:end-1,2:end-1)); axis ij equal tight; box on; cb = colorbar;
        set(cb,TL{:},TS{:}); set(gca,TL{:},TS{:}); title(['res. $P$'],TX{:},FS{:}); set(gca,'YTickLabel',[]); 
    end
    
    end
    
       
    % plot model history
    if plot_cv
        figure(8); clf;
        subplot(3,1,1);
        plot(hist.time/hr,hist.DM./hist.sumM,'k-',LW{:}); hold on; axis tight; box on;
        ylabel('consv. $M$',TX{:},FS{:}); set(gca,TL{:},TS{:},'XTickLabel',[]);
        subplot(3,1,2);
        plot(hist.time/hr,hist.DH./hist.sumH,'k-',LW{:}); hold on; axis tight; box on;
        ylabel('consv. $H$',TX{:},FS{:}); set(gca,TL{:},TS{:},'XTickLabel',[]);
        subplot(3,1,3);
        plot(hist.time/hr,hist.DC./hist.sumC,'-',LW{:}); hold on; axis tight; box on;
        ylabel('consv. $C$',TX{:},FS{:}); set(gca,TL{:},TS{:},'XTickLabel',[]);
        xlabel('Time [hr]',TX{:},FS{:});
        
        figure(9); clf;
        subplot(3,1,1);
        plot(hist.time/hr,hist.EM,'k-',LW{:}); hold on; axis tight; box on;
        ylabel('error $M$',TX{:},FS{:}); set(gca,TL{:},TS{:},'XTickLabel',[]);
        subplot(3,1,2);
        plot(hist.time/hr,hist.EH,'k-',LW{:}); hold on; axis tight; box on;
        ylabel('error $H$',TX{:},FS{:}); set(gca,TL{:},TS{:},'XTickLabel',[]);
        subplot(3,1,3);
        plot(hist.time/hr,hist.EC,'-',LW{:}); hold on; axis tight; box on;
        ylabel('error $C$',TX{:},FS{:}); set(gca,TL{:},TS{:},'XTickLabel',[]);
        xlabel('Time [hr]',TX{:},FS{:});
    end
    
    drawnow

end

% save output to file
if save_op
    if plot_op
        if Nx <= 10 && Nz <= 10  % print 0D plots
            name = [opdir,'/',runID,'/',runID,'_tch_',num2str(floor(step/nop))];
            print(fh1,name,'-dpng','-r300','-opengl');
            name = [opdir,'/',runID,'/',runID,'_aux_',num2str(floor(step/nop))];
            print(fh2,name,'-dpng','-r300','-opengl');
            name = [opdir,'/',runID,'/',runID,'_oxd_',num2str(floor(step/nop))];
            print(fh3,name,'-dpng','-r300','-opengl');
        elseif Nx <= 10  % create 1D plots
            name = [opdir,'/',runID,'/',runID,'_sol_',num2str(floor(step/nop))];
            print(fh1,name,'-dpng','-r300','-opengl');
            name = [opdir,'/',runID,'/',runID,'_aux_',num2str(floor(step/nop))];
            print(fh2,name,'-dpng','-r300','-opengl');
            name = [opdir,'/',runID,'/',runID,'_oxd_',num2str(floor(step/nop))];
            print(fh3,name,'-dpng','-r300','-opengl');
            name = [opdir,'/',runID,'/',runID,'_gch_',num2str(floor(step/nop))];
            print(fh4,name,'-dpng','-r300','-opengl');
        else
            name = [opdir,'/',runID,'/',runID,'_sol_',num2str(floor(step/nop))];
            print(fh1,name,'-dpng','-r300','-opengl');
            name = [opdir,'/',runID,'/',runID,'_cmp_',num2str(floor(step/nop))];
            print(fh2,name,'-dpng','-r300','-opengl');
            name = [opdir,'/',runID,'/',runID,'_oxd_',num2str(floor(step/nop))];
            print(fh3,name,'-dpng','-r300','-opengl');
            name = [opdir,'/',runID,'/',runID,'_aux_',num2str(floor(step/nop))];
            print(fh4,name,'-dpng','-r300','-opengl');
            name = [opdir,'/',runID,'/',runID,'_gch',num2str(floor(step/nop))];
            print(fh5,name,'-dpng','-r300','-opengl');
        end
    end
    
    name = [opdir,'/',runID,'/',runID,'_',num2str(floor(step/nop))];
    save(name,'U','W','P','Pt','x','m','chi','mu','H','C','T','c','cm','cx','IT','CT','SI','RIP','RID','it','ct','si','rip','rid','dHdt','dCdt','dITdt','dCTdt','dSIdt','dxdt','Gx','rho','eta','exx','ezz','exz','txx','tzz','txz','eII','tII','dt','time','step','hist','VolSrc','wx');
    name = [opdir,'/',runID,'/',runID,'_cont'];
    save(name,'U','W','P','Pt','x','m','chi','mu','H','C','T','c','cm','cx','IT','CT','SI','RIP','RID','it','ct','si','rip','rid','dHdt','dCdt','dITdt','dCTdt','dSIdt','dxdt','Gx','rho','eta','exx','ezz','exz','txx','tzz','txz','eII','tII','dt','time','step','hist','VolSrc','wx');
    
    if step == 0
        logfile = [opdir,'/',runID,'/',runID,'.log'];
        if exist(logfile,'file'); delete(logfile); end
        diary(logfile)
    end
end

clear fh1 fh2 fh3 fh4 fh5 fh6 fh7
    