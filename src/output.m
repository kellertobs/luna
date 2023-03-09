
% prepare for plotting
TX = {'Interpreter','Latex'}; FS = {'FontSize',12};
TL = {'TickLabelInterpreter','Latex'}; TS = {'FontSize',10};
UN = {'Units','Centimeters'};
CL = {'Color',[0.0 0.0 0.0],[0.80 0.15 0.10],[0.10 0.15 0.65],[0.45 0.60 0.95]};
LW = {'LineWidth',2};
if plot_op
    VIS = {'Visible','on'};
else
    VIS = {'Visible','off'};
end

if Nx <= 10 && Nz <= 10  % create 0D plots

    if ~exist('fh1','var'); fh1 = figure(VIS{:});
    else; set(0, 'CurrentFigure', fh1); clf;
    end
    subplot(4,1,1)
    plot(hist.time/yr,hist.mu (:,2)*100.*(hist.mu (:,2)>1e-9),CL{[1,3]},LW{:}); axis xy tight; box on; hold on;
    plot(hist.time/yr,hist.chi(:,2)*100.*(hist.chi(:,2)>1e-9),CL{[1,4]},LW{:});
    title(['$\mu$, $\chi$ [vol\%]'],TX{:},FS{:}); set(gca,TL{:},TS{:});
    subplot(4,1,2)
    plot(hist.time/yr,squeeze(hist.c(:,2,:)),LW{:});  axis xy tight; box on;
    title('$\bar{c}$ [wt\%]',TX{:},FS{:}); set(gca,TL{:},TS{:});
    subplot(4,1,3)
    plot(hist.time/yr,squeeze(hist.cx(:,2,:)),LW{:});  axis xy tight; box on;
    title('$c^x$ [wt\%]',TX{:},FS{:}); set(gca,TL{:},TS{:});
    legend(cal.cmpStr,TX{:},FS{:},'location','best');
    subplot(4,1,4)
    plot(hist.time/yr,squeeze(hist.cm(:,2,:)),LW{:});  axis xy tight; box on;
    title('$c^m$ [wt\%]',TX{:},FS{:}); set(gca,TL{:},TS{:});
    xlabel('Time [yr]',TX{:},FS{:});

    if ~exist('fh2','var'); fh2 = figure(VIS{:});
    else; set(0, 'CurrentFigure', fh2); clf;
    end
    subplot(4,1,1)
    plot(hist.time/yr,hist.T(:,2)-273.15,CL{[1,2]},LW{:}); axis xy tight; box on;
    title('$T [^\circ$C]',TX{:},FS{:}); set(gca,TL{:},TS{:});
    subplot(4,1,2)
    plot(hist.time/yr,hist.rhom(:,2),'-',CL{[1,3]},LW{:}); axis xy tight; box on; hold on;
    plot(hist.time/yr,hist.rhox(:,2),'-',CL{[1,4]},LW{:});
    plot(hist.time/yr,hist.rho (:,2),'-',CL{[1,2]},LW{:});
    title('$\bar{\rho}$ [kg/m$^3$]',TX{:},FS{:}); set(gca,TL{:},TS{:});
    subplot(4,1,3)
    semilogy(hist.time/yr,hist.etam(:,2),'r-',LW{:}); axis xy tight; box on; hold on;
    semilogy(hist.time/yr,hist.eta (:,2),'k-',LW{:});
    title('$\bar{\eta}$ [Pas]',TX{:},FS{:}); set(gca,TL{:},TS{:});
    subplot(4,1,4)
    plot(hist.time/yr,hist.Gx(:,2)./hist.rho(:,2)*yr*100.*(hist.chi(:,2)>1e-9),CL{[1,4]},LW{:}); axis xy tight; box on; hold on;
    title('$\Gamma_x/\bar{\rho}$ [wt\%/yr]',TX{:},FS{:}); set(gca,TL{:},TS{:});
    xlabel('Time [yr]',TX{:},FS{:});

    if ~exist('fh3','var'); fh3 = figure(VIS{:});
    else; set(0, 'CurrentFigure', fh3); clf;
    end
    subplot(3,1,1)
    plot(hist.time/yr,squeeze(hist.c_oxd(:,2,:)),LW{:});  axis xy tight; box on;
    title('$\bar{c}_\mathrm{oxd}$ [wt\%]',TX{:},FS{:}); set(gca,TL{:},TS{:});
    subplot(3,1,2)
    plot(hist.time/yr,squeeze(hist.cx_oxd(:,2,:)),LW{:});  axis xy tight; box on;
    title('$c^x_\mathrm{oxd}$ [wt\%]',TX{:},FS{:}); set(gca,TL{:},TS{:});
    legend(cal.oxdStr,TX{:},FS{:},'location','best');
    subplot(3,1,3)
    plot(hist.time/yr,squeeze(hist.cm_oxd(:,2,:)),LW{:});  axis xy tight; box on;
    title('$c^m_\mathrm{oxd}$ [wt\%]',TX{:},FS{:}); set(gca,TL{:},TS{:});
    xlabel('Time [yr]',TX{:},FS{:});

elseif Nx <= 10  % create 1D plots

    if ~exist('fh1','var'); fh1 = figure(VIS{:});
    else; set(0, 'CurrentFigure', fh1); clf;
    end
    sgtitle(['time = ',num2str(time/yr,3),' [yr]'],TX{:},FS{:},'Color','k');
    subplot(1,5,1)
    plot(reshape(cal.Tliq,Nz-2,Nx-2),Zc(2:end-1).'./1e3,CL{[1,3]},LW{:}); axis ij tight; box on; hold on;
    plot(reshape(cal.Tsol,Nz-2,Nx-2),Zc(2:end-1).'./1e3,CL{[1,4]},LW{:});
    plot(T(2:end-1,2:end-1)-273.15,Zc(2:end-1).'./1e3,CL{[1,2]},LW{:}); 
    title('$T [^\circ$C]',TX{:},FS{:}); ylabel('Depth [km]',TX{:},FS{:}); set(gca,TL{:},TS{:});
    subplot(1,5,2)
    plot(mu (2:end-1,2:end-1)*100.*(mu (2:end-1,2:end-1)>1e-9),Zc(2:end-1).'./1e3,CL{[1,3]},LW{:}); axis ij tight; box on; hold on;
    plot(chi(2:end-1,2:end-1)*100.*(chi(2:end-1,2:end-1)>1e-9),Zc(2:end-1).'./1e3,CL{[1,4]},LW{:});
    title(['$\mu$, $\chi$ [vol\%]'],TX{:},FS{:}); set(gca,TL{:},TS{:}); set(gca,'YTickLabel',[]);
    subplot(1,5,3)
    plot(squeeze(c(2:end-1,2:end-1,:))*100,Zc(2:end-1).'./1e3,LW{:}); axis ij tight; box on; hold on;
    title('$\bar{c}$ [wt\%]',TX{:},FS{:}); set(gca,TL{:},TS{:}); set(gca,'YTickLabel',[]);
    subplot(1,5,4)
    plot(squeeze(cx(2:end-1,2:end-1,:))*100,Zc(2:end-1).'./1e3,LW{:}); axis ij tight; box on; hold on;
    title('$c^x$ [wt\%]',TX{:},FS{:}); set(gca,TL{:},TS{:}); set(gca,'YTickLabel',[]);
    subplot(1,5,5)
    plot(squeeze(cm(2:end-1,2:end-1,:))*100,Zc(2:end-1).'./1e3,LW{:}); axis ij tight; box on; hold on;
    title('$c^m$ [wt\%]',TX{:},FS{:}); set(gca,TL{:},TS{:}); set(gca,'YTickLabel',[]);
    legend(cal.cmpStr,TX{:},FS{:},'location','east');

    if ~exist('fh2','var'); fh2 = figure(VIS{:});
    else; set(0, 'CurrentFigure', fh2); clf;
    end
    sgtitle(['time = ',num2str(time/yr,3),' [yr]'],TX{:},FS{:},'Color','k');
    subplot(1,5,1)
    plot(rhox(2:end-1,2:end-1),Zc(2:end-1).'./1e3,CL{[1,4]},LW{:}); axis ij tight; box on; hold on;
    plot(rhom(2:end-1,2:end-1),Zc(2:end-1).'./1e3,CL{[1,3]},LW{:});
    plot(rho (2:end-1,2:end-1),Zc(2:end-1).'./1e3,CL{[1,2]},LW{:});
    title('$\bar{\rho}$ [kg/m$^3$]',TX{:},FS{:}); ylabel('Depth [km]',TX{:},FS{:}); set(gca,TL{:},TS{:});
    subplot(1,5,2)
    semilogx(min(eta(2:end-1,2:end-1),etam(2:end-1,2:end-1)),Zc(2:end-1).'./1e3,CL{[1,3]},LW{:}); axis ij tight; box on; hold on;
    semilogx(eta(2:end-1,2:end-1),Zc(2:end-1).'./1e3,CL{[1,2]},LW{:});
    title('$\bar{\eta}$ [Pas]',TX{:},FS{:}); set(gca,TL{:},TS{:}); set(gca,'YTickLabel',[]);
    subplot(1,5,3)
    plot(Gx./rho(2:end-1,2:end-1)*yr*100.*(chi(2:end-1,2:end-1)>1e-9),Zc(2:end-1).'./1e3,CL{[1,4]},LW{:}); axis ij tight; box on; hold on;
    title('$\Gamma_x/\bar{\rho}$ [wt/yr]',TX{:},FS{:}); set(gca,TL{:},TS{:}); set(gca,'YTickLabel',[]);
    subplot(1,5,4)
    plot(-( mu(1:end-1,2:end-1)+ mu(2:end,2:end-1))/2.*wm(:,2:end-1)*hr,Zf.'./1e3,CL{[1,3]},LW{:}); axis ij tight; box on; hold on;
    plot(-(chi(1:end-1,2:end-1)+chi(2:end,2:end-1))/2.*wx(:,2:end-1)*hr,Zf.'./1e3,CL{[1,4]},LW{:});
    plot(-                                              W(:,2:end-1)*hr,Zf.'./1e3,CL{[1,2]},LW{:});
    title('$W$, $w_\Delta^x$ [m/hr]',TX{:},FS{:}); set(gca,TL{:},TS{:}); set(gca,'YTickLabel',[]);
    subplot(1,5,5)
    plot(P(2:end-1,2:end-1),Zc(2:end-1).'./1e3,CL{[1,2]},LW{:}); axis ij tight; box on;
    title('$P$ [Pa]',TX{:},FS{:}); set(gca,TL{:},TS{:}); set(gca,'YTickLabel',[]);

    if ~exist('fh3','var'); fh3 = figure(VIS{:});
    else; set(0, 'CurrentFigure', fh3); clf;
    end
    sgtitle(['time = ',num2str(time/yr,3),' [yr]'],TX{:},FS{:},'Color','k');
    subplot(1,3,1)
    plot(squeeze(c_oxd(2:end-1,2:end-1,:)),Zc(2:end-1).'./1e3,LW{:}); axis ij tight; box on; hold on;
    title('$\bar{c}_\mathrm{oxd}$ [wt\%]',TX{:},FS{:}); set(gca,TL{:},TS{:});
    subplot(1,3,2)
    plot(squeeze(cx_oxd(2:end-1,2:end-1,:)),Zc(2:end-1).'./1e3,LW{:}); axis ij tight; box on; hold on;
    title('$c^x_\mathrm{oxd}$ [wt\%]',TX{:},FS{:}); set(gca,TL{:},TS{:}); set(gca,'YTickLabel',[]);
    subplot(1,3,3)
    plot(squeeze(cm_oxd(2:end-1,2:end-1,:)),Zc(2:end-1).'./1e3,LW{:}); axis ij tight; box on; hold on;
    title('$c^m_\mathrm{oxd}$ [wt\%]',TX{:},FS{:}); set(gca,TL{:},TS{:}); set(gca,'YTickLabel',[]);
    legend(cal.oxdStr,TX{:},FS{:},'location','east');

else % create 2D plots

    % set axis and border dimensions
    axh = 6.00; axw = axh*L/D;
    ahs = 0.30; avs = 0.0;
    axb = 0.40; axt = 0.6;
    axl = 1.40; axr = 0.4;

    % initialize figures and axes
    if ~exist('fh1','var'); fh1 = figure(VIS{:});
        colormap(ocean);
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
    end

    if ~exist('fh2','var'); fh2 = figure(VIS{:});
        colormap(ocean);
        fh = axb + 2*axh + 1*avs + axt;
        fw = axl + 2*axw + 1*ahs + axr;
        set(fh2,UN{:},'Position',[2 2 fw fh]);
        set(fh2,'PaperUnits','Centimeters','PaperPosition',[0 0 fw fh],'PaperSize',[fw fh]);
        set(fh2,'Color','w','InvertHardcopy','off');
        set(fh2,'Resize','off');
        ax(21) = axes(UN{:},'position',[axl+0*axw+0*ahs axb+1*axh+1*avs axw axh]);
        ax(22) = axes(UN{:},'position',[axl+1*axw+1*ahs axb+1*axh+1*avs axw axh]);
        ax(23) = axes(UN{:},'position',[axl+0*axw+0*ahs axb+0*axh+0*avs axw axh]);
        ax(24) = axes(UN{:},'position',[axl+1*axw+1*ahs axb+0*axh+0*avs axw axh]);
    end

    if ~exist('fh3','var'); fh3 = figure(VIS{:});
        colormap(ocean);
        fh = axb + 2*axh + 1*avs + axt;
        fw = axl + 3*axw + 2*ahs + axr;
        set(fh3,UN{:},'Position',[3 3 fw fh]);
        set(fh3,'PaperUnits','Centimeters','PaperPosition',[0 0 fw fh],'PaperSize',[fw fh]);
        set(fh3,'Color','w','InvertHardcopy','off');
        set(fh3,'Resize','off');
        ax(31) = axes(UN{:},'position',[axl+0*axw+0*ahs axb+1*axh+1*avs axw axh]);
        ax(32) = axes(UN{:},'position',[axl+1*axw+1*ahs axb+1*axh+1*avs axw axh]);
        ax(33) = axes(UN{:},'position',[axl+2*axw+2*ahs axb+1*axh+1*avs axw axh]);
        ax(34) = axes(UN{:},'position',[axl+0*axw+0*ahs axb+0*axh+0*avs axw axh]);
        ax(35) = axes(UN{:},'position',[axl+1*axw+1*ahs axb+0*axh+0*avs axw axh]);
        ax(36) = axes(UN{:},'position',[axl+2*axw+2*ahs axb+0*axh+0*avs axw axh]);
    end

    if ~exist('fh4','var'); fh4 = figure(VIS{:});
        colormap(ocean);
        fh = axb + 2*axh + 1*avs + axt;
        fw = axl + 3*axw + 2*ahs + axr;
        set(fh4,UN{:},'Position',[4 4 fw fh]);
        set(fh4,'PaperUnits','Centimeters','PaperPosition',[0 0 fw fh],'PaperSize',[fw fh]);
        set(fh4,'Color','w','InvertHardcopy','off');
        set(fh4,'Resize','off');
        ax(41) = axes(UN{:},'position',[axl+0*axw+0*ahs axb+1*axh+1*avs axw axh]);
        ax(42) = axes(UN{:},'position',[axl+1*axw+1*ahs axb+1*axh+1*avs axw axh]);
        ax(43) = axes(UN{:},'position',[axl+2*axw+2*ahs axb+1*axh+1*avs axw axh]);
        ax(44) = axes(UN{:},'position',[axl+0*axw+0*ahs axb+0*axh+0*avs axw axh]);
        ax(45) = axes(UN{:},'position',[axl+1*axw+1*ahs axb+0*axh+0*avs axw axh]);
        ax(46) = axes(UN{:},'position',[axl+2*axw+2*ahs axb+0*axh+0*avs axw axh]);
    end

    if ~exist('fh5','var'); fh5 = figure(VIS{:});
        colormap(ocean);
        fh = axb + 2*axh + 1*avs + axt;
        fw = axl + 3*axw + 2*ahs + axr;
        set(fh5,UN{:},'Position',[5 5 fw fh]);
        set(fh5,'PaperUnits','Centimeters','PaperPosition',[0 0 fw fh],'PaperSize',[fw fh]);
        set(fh5,'Color','w','InvertHardcopy','off');
        set(fh5,'Resize','off');
        ax(51) = axes(UN{:},'position',[axl+0*axw+0*ahs axb+1*axh+1*avs axw axh]);
        ax(52) = axes(UN{:},'position',[axl+1*axw+1*ahs axb+1*axh+1*avs axw axh]);
        ax(53) = axes(UN{:},'position',[axl+2*axw+2*ahs axb+1*axh+1*avs axw axh]);
        ax(54) = axes(UN{:},'position',[axl+0*axw+0*ahs axb+0*axh+0*avs axw axh]);
        ax(55) = axes(UN{:},'position',[axl+1*axw+1*ahs axb+0*axh+0*avs axw axh]);
        ax(56) = axes(UN{:},'position',[axl+2*axw+2*ahs axb+0*axh+0*avs axw axh]);
    end

    if ~exist('fh6','var'); fh6 = figure(VIS{:});
        colormap(ocean);
        fh = axb + 3*axh + 2*avs + axt;
        fw = axl + 3*axw + 2*ahs + axr;
        set(fh6,UN{:},'Position',[6 6 fw fh]);
        set(fh6,'PaperUnits','Centimeters','PaperPosition',[0 0 fw fh],'PaperSize',[fw fh]);
        set(fh6,'Color','w','InvertHardcopy','off');
        set(fh6,'Resize','off');
        ax(61) = axes(UN{:},'position',[axl+0*axw+0*ahs axb+2*axh+2*avs axw axh]);
        ax(62) = axes(UN{:},'position',[axl+1*axw+1*ahs axb+2*axh+2*avs axw axh]);
        ax(63) = axes(UN{:},'position',[axl+2*axw+2*ahs axb+2*axh+2*avs axw axh]);
        ax(64) = axes(UN{:},'position',[axl+0*axw+0*ahs axb+1*axh+1*avs axw axh]);
        ax(65) = axes(UN{:},'position',[axl+1*axw+1*ahs axb+1*axh+1*avs axw axh]);
        ax(66) = axes(UN{:},'position',[axl+2*axw+2*ahs axb+1*axh+1*avs axw axh]);
        ax(67) = axes(UN{:},'position',[axl+0*axw+0*ahs axb+0*axh+0*avs axw axh]);
        ax(68) = axes(UN{:},'position',[axl+1*axw+1*ahs axb+0*axh+0*avs axw axh]);
        ax(69) = axes(UN{:},'position',[axl+2*axw+2*ahs axb+0*axh+0*avs axw axh]);
    end

    if plot_cv
        if ~exist('fh7','var'); fh7 = figure(VIS{:});
            colormap(ocean);
            fh = axb + 1*axh + 0*avs + axt;
            fw = axl + 3*axw + 2*ahs + axr;
            set(fh7,UN{:},'Position',[7 7 fw fh]);
            set(fh7,'PaperUnits','Centimeters','PaperPosition',[0 0 fw fh],'PaperSize',[fw fh]);
            set(fh7,'Color','w','InvertHardcopy','off');
            set(fh7,'Resize','off');
            ax(71) = axes(UN{:},'position',[axl+0*axw+0*ahs axb+0*axh+0*avs axw axh]);
            ax(72) = axes(UN{:},'position',[axl+1*axw+1*ahs axb+0*axh+0*avs axw axh]);
            ax(73) = axes(UN{:},'position',[axl+2*axw+2*ahs axb+0*axh+0*avs axw axh]);
        end
    end

    % plot velocity-pressure solution in Fig. 1
    set(0,'CurrentFigure',fh1)
    set(fh1,'CurrentAxes',ax(11));
    imagesc(Xc(2:end-1)./1e3,Zc(2:end-1)./1e3,-W(:      ,2:end-1).*hr); axis ij equal tight; box on; cb = colorbar;
    set(cb,TL{:},TS{:}); set(gca,TL{:},TS{:}); title(['$W$ [m/hr]'],TX{:},FS{:}); ylabel('Depth [km]',TX{:},FS{:}); ylabel('Depth [km]',TX{:},FS{:});  set(gca,'XTickLabel',[]);
    set(fh1,'CurrentAxes',ax(12));
    imagesc(Xc(2:end-1)./1e3,Zc(2:end-1)./1e3, U(2:end-1,:      ).*hr); axis ij equal tight; box on; cb = colorbar;
    set(cb,TL{:},TS{:}); set(gca,TL{:},TS{:}); title(['$U$ [m/hr]'],TX{:},FS{:}); set(gca,'XTickLabel',[],'YTickLabel',[]);
    set(fh1,'CurrentAxes',ax(13));
    imagesc(Xc(2:end-1)./1e3,Zc(2:end-1)./1e3, P(2:end-1,2:end-1)); axis ij equal tight; box on; cb = colorbar;
    set(cb,TL{:},TS{:}); set(gca,TL{:},TS{:}); title(['$P$ [Pa]'],TX{:},FS{:}); ylabel('Depth [km]',TX{:},FS{:}); xlabel('Width [km]',TX{:},FS{:});
    set(fh1,'CurrentAxes',ax(14));
    imagesc(Xc(2:end-1)./1e3,Zc(2:end-1)./1e3,T(2:end-1,2:end-1)-273.15); axis ij equal tight; box on; cb = colorbar;
    set(cb,TL{:},TS{:}); set(gca,TL{:},TS{:}); title(['$T [^\circ$C]'],TX{:},FS{:}); xlabel('Width [km]',TX{:},FS{:}); set(gca,'YTickLabel',[]);
    sgtitle(['time = ',num2str(time/yr,3),' [yr]'],TX{:},FS{:},'Color','k');

    % plot pseudo-component composition in Fig. 2
    set(0,'CurrentFigure',fh2)
    set(fh2,'CurrentAxes',ax(21));
    imagesc(Xc(2:end-1)./1e3,Zc(2:end-1)./1e3,c(2:end-1,2:end-1,1)); axis ij equal tight; box on; cb = colorbar;
    set(cb,TL{:},TS{:}); set(gca,TL{:},TS{:}); title([cal.cmpStr{1},' [wt\%]'],TX{:},FS{:}); set(gca,'XTickLabel',[]); ylabel('Depth [km]',TX{:},FS{:});
    set(fh2,'CurrentAxes',ax(22));
    imagesc(Xc(2:end-1)./1e3,Zc(2:end-1)./1e3,c(2:end-1,2:end-1,2)); axis ij equal tight; box on; cb = colorbar;
    set(cb,TL{:},TS{:}); set(gca,TL{:},TS{:}); title([cal.cmpStr{2},' [wt\%]'],TX{:},FS{:}); set(gca,'XTickLabel',[],'YTickLabel',[]);
    set(fh2,'CurrentAxes',ax(23));
    imagesc(Xc(2:end-1)./1e3,Zc(2:end-1)./1e3,c(2:end-1,2:end-1,3)); axis ij equal tight; box on; cb = colorbar;
    set(cb,TL{:},TS{:}); set(gca,TL{:},TS{:}); title([cal.cmpStr{3},' [wt\%]'],TX{:},FS{:}); ylabel('Depth [km]',TX{:},FS{:});
    set(fh2,'CurrentAxes',ax(24));
    imagesc(Xc(2:end-1)./1e3,Zc(2:end-1)./1e3,c(2:end-1,2:end-1,4)); axis ij equal tight; box on; cb = colorbar;
    set(cb,TL{:},TS{:}); set(gca,TL{:},TS{:}); title([cal.cmpStr{4},' [wt\%]'],TX{:},FS{:}); set(gca,'YTickLabel',[]); xlabel('Width [km]',TX{:},FS{:});
    sgtitle(['time = ',num2str(time/yr,3),' [yr]'],TX{:},FS{:},'Color','k');

    % plot major oxide composition in Fig. 3
    set(0,'CurrentFigure',fh3)
    set(fh3,'CurrentAxes',ax(31));
    imagesc(Xc(2:end-1)./1e3,Zc(2:end-1)./1e3,c_oxd(2:end-1,2:end-1,1)); axis ij equal tight; box on; cb = colorbar;
    set(cb,TL{:},TS{:}); set(gca,TL{:},TS{:}); title([cal.oxdStr{1},' [wt\%]'],TX{:},FS{:}); set(gca,'XTickLabel',[]); ylabel('Depth [km]',TX{:},FS{:});
    set(fh3,'CurrentAxes',ax(32));
    imagesc(Xc(2:end-1)./1e3,Zc(2:end-1)./1e3,c_oxd(2:end-1,2:end-1,2)); axis ij equal tight; box on; cb = colorbar;
    set(cb,TL{:},TS{:}); set(gca,TL{:},TS{:}); title([cal.oxdStr{2},' [wt\%]'],TX{:},FS{:}); set(gca,'XTickLabel',[],'YTickLabel',[]);
    set(fh3,'CurrentAxes',ax(33));
    imagesc(Xc(2:end-1)./1e3,Zc(2:end-1)./1e3,c_oxd(2:end-1,2:end-1,3)); axis ij equal tight; box on; cb = colorbar;
    set(cb,TL{:},TS{:}); set(gca,TL{:},TS{:}); title([cal.oxdStr{3},' [wt\%]'],TX{:},FS{:}); set(gca,'XTickLabel',[],'YTickLabel',[]);
    set(fh3,'CurrentAxes',ax(34));
    imagesc(Xc(2:end-1)./1e3,Zc(2:end-1)./1e3,c_oxd(2:end-1,2:end-1,4)); axis ij equal tight; box on; cb = colorbar;
    set(cb,TL{:},TS{:}); set(gca,TL{:},TS{:}); title([cal.oxdStr{4},' [wt\%]'],TX{:},FS{:}); ylabel('Depth [km]',TX{:},FS{:});
    set(fh3,'CurrentAxes',ax(35));
    imagesc(Xc(2:end-1)./1e3,Zc(2:end-1)./1e3,c_oxd(2:end-1,2:end-1,5)); axis ij equal tight; box on; cb = colorbar;
    set(cb,TL{:},TS{:}); set(gca,TL{:},TS{:}); title([cal.oxdStr{5},' [wt\%]'],TX{:},FS{:}); set(gca,'YTickLabel',[]); xlabel('Width [km]',TX{:},FS{:});
    set(fh3,'CurrentAxes',ax(36));
    imagesc(Xc(2:end-1)./1e3,Zc(2:end-1)./1e3,c_oxd(2:end-1,2:end-1,6)+c_oxd(2:end-1,2:end-1,7)); axis ij equal tight; box on; cb = colorbar;
    set(cb,TL{:},TS{:}); set(gca,TL{:},TS{:}); title([cal.oxdStr{6},'+',cal.oxdStr{7},' [wt\%]'],TX{:},FS{:}); set(gca,'YTickLabel',[]);
    sgtitle(['time = ',num2str(time/yr,3),' [yr]'],TX{:},FS{:},'Color','k');

    % plot phase fractions and reaction rates in Fig. 4
    set(0,'CurrentFigure',fh4)
    set(fh4,'CurrentAxes',ax(41));
    imagesc(Xc(2:end-1)./1e3,Zc(2:end-1)./1e3,chi(2:end-1,2:end-1).*100.*(chi(2:end-1,2:end-1)>1e-9) ); axis ij equal tight; box on; cb = colorbar;
    set(cb,TL{:},TS{:}); set(gca,TL{:},TS{:}); title(['$\chi$ [vol\%]'],TX{:},FS{:}); set(gca,'XTickLabel',[]); ylabel('Depth [km]',TX{:},FS{:});
    set(fh4,'CurrentAxes',ax(42));
    imagesc(Xc(2:end-1)./1e3,Zc(2:end-1)./1e3,Gx./rho(2:end-1,2:end-1)*yr*100.*(chi(2:end-1,2:end-1)>1e-9)); axis ij equal tight; box on; cb = colorbar;
    set(cb,TL{:},TS{:}); set(gca,TL{:},TS{:}); title(['$\Gamma_x/\bar{\rho}$ [wt\%/yr]'],TX{:},FS{:}); set(gca,'XTickLabel',[],'YTickLabel',[]);
    set(fh4,'CurrentAxes',ax(43));
    imagesc(Xc(2:end-1)./1e3,Zc(2:end-1)./1e3,-(chi(1:end-1,2:end-1)+chi(2:end,2:end-1))/2.*wx(:,2:end-1).*hr); axis ij equal tight; box on; cb = colorbar;
    set(cb,TL{:},TS{:}); set(gca,TL{:},TS{:}); title(['$w_\Delta^x$ [m/hr]'],TX{:},FS{:}); set(gca,'XTickLabel',[],'YTickLabel',[]);
    set(fh4,'CurrentAxes',ax(44));
    imagesc(Xc(2:end-1)./1e3,Zc(2:end-1)./1e3,rho(2:end-1,2:end-1)-mean(rho(2:end-1,2:end-1),2)); axis ij equal tight; box on; cb = colorbar;
    set(cb,TL{:},TS{:}); set(gca,TL{:},TS{:}); title(['$\bar{\rho}-\langle\bar{\rho}\rangle_h$ [kg/m$^3$]'],TX{:},FS{:}); ylabel('Depth [km]',TX{:},FS{:});
    set(fh4,'CurrentAxes',ax(45));
    imagesc(Xc(2:end-1)./1e3,Zc(2:end-1)./1e3,log10(etam(2:end-1,2:end-1))); axis ij equal tight; box on; cb = colorbar;
    set(cb,TL{:},TS{:}); set(gca,TL{:},TS{:}); title(['$\eta^m$ [log$_{10}$ Pas]'],TX{:},FS{:}); set(gca,'YTickLabel',[]); xlabel('Width [km]',TX{:},FS{:});
    set(fh4,'CurrentAxes',ax(46));
    imagesc(Xc(2:end-1)./1e3,Zc(2:end-1)./1e3,log10(eta(2:end-1,2:end-1))); axis ij equal tight; box on; cb = colorbar;
    set(cb,TL{:},TS{:}); set(gca,TL{:},TS{:}); title(['$\bar{\eta}$ [log$_{10}$ Pas]'],TX{:},FS{:}); set(gca,'YTickLabel',[]);
    sgtitle(['time = ',num2str(time/yr,3),' [yr]'],TX{:},FS{:},'Color','k');
    
    % plot geochemical variables in Fig. 5
    set(0,'CurrentFigure',fh5)
    set(fh5,'CurrentAxes',ax(51));
    imagesc(Xc(2:end-1)./1e3,Zc(2:end-1)./1e3,te(2:end-1,2:end-1,1)); axis ij equal tight; box on; cb = colorbar;
    set(cb,TL{:},TS{:}); set(gca,TL{:},TS{:}); title(['trace element 1'],TX{:},FS{:}); set(gca,'XTickLabel',[]); ylabel('Depth [m]',TX{:},FS{:});
    set(fh5,'CurrentAxes',ax(52));
    imagesc(Xc(2:end-1)./1e3,Zc(2:end-1)./1e3,te(2:end-1,2:end-1,2)); axis ij equal tight; box on; cb = colorbar;
    set(cb,TL{:},TS{:}); set(gca,TL{:},TS{:}); title(['trace element 2'],TX{:},FS{:}); set(gca,'XTickLabel',[],'YTickLabel',[]);
    set(fh5,'CurrentAxes',ax(53));
    imagesc(Xc(2:end-1)./1e3,Zc(2:end-1)./1e3,ir(2:end-1,2:end-1,1)); axis ij equal tight; box on; cb = colorbar;
    set(cb,TL{:},TS{:}); set(gca,TL{:},TS{:}); title(['isotope ratio 1'],TX{:},FS{:}); set(gca,'XTickLabel',[],'YTickLabel',[]);
    set(fh5,'CurrentAxes',ax(54));
    imagesc(Xc(2:end-1)./1e3,Zc(2:end-1)./1e3,te(2:end-1,2:end-1,3)); axis ij equal tight; box on; cb = colorbar;
    set(cb,TL{:},TS{:}); set(gca,TL{:},TS{:}); title(['trace element 3'],TX{:},FS{:}); ylabel('Depth [m]',TX{:},FS{:});
    set(fh5,'CurrentAxes',ax(55));
    imagesc(Xc(2:end-1)./1e3,Zc(2:end-1)./1e3,te(2:end-1,2:end-1,4)); axis ij equal tight; box on; cb = colorbar;
    set(cb,TL{:},TS{:}); set(gca,TL{:},TS{:}); title(['trace element 4'],TX{:},FS{:}); set(gca,'YTickLabel',[]); xlabel('Width [m]',TX{:},FS{:});
    set(fh5,'CurrentAxes',ax(56));
    imagesc(Xc(2:end-1)./1e3,Zc(2:end-1)./1e3,ir(2:end-1,2:end-1,2)); axis ij equal tight; box on; cb = colorbar;
    set(cb,TL{:},TS{:}); set(gca,TL{:},TS{:}); title(['isotope ratio 2'],TX{:},FS{:}); set(gca,'YTickLabel',[]);
    sgtitle(['time = ',num2str(time/yr,3),' [yr]'],TX{:},FS{:},'Color','k');

    % plot major oxide composition in Fig. 6
    set(0,'CurrentFigure',fh6)
    set(fh6,'CurrentAxes',ax(61));
    imagesc(Xc(2:end-1)./1e3,Zc(2:end-1)./1e3,x(2:end-1,2:end-1).*cx_mnr(2:end-1,2:end-1,1).*100); axis ij equal tight; box on; cb = colorbar;
    set(cb,TL{:},TS{:}); set(gca,TL{:},TS{:}); title([cal.mnrStr{1},' [wt\%]'],TX{:},FS{:}); set(gca,'XTickLabel',[]); ylabel('Depth [km]',TX{:},FS{:});
    set(fh6,'CurrentAxes',ax(62));
    imagesc(Xc(2:end-1)./1e3,Zc(2:end-1)./1e3,x(2:end-1,2:end-1).*cx_mnr(2:end-1,2:end-1,2).*100); axis ij equal tight; box on; cb = colorbar;
    set(cb,TL{:},TS{:}); set(gca,TL{:},TS{:}); title([cal.mnrStr{2},' [wt\%]'],TX{:},FS{:}); set(gca,'XTickLabel',[],'YTickLabel',[]);
    set(fh6,'CurrentAxes',ax(63));
    imagesc(Xc(2:end-1)./1e3,Zc(2:end-1)./1e3,x(2:end-1,2:end-1).*cx_mnr(2:end-1,2:end-1,3).*100); axis ij equal tight; box on; cb = colorbar;
    set(cb,TL{:},TS{:}); set(gca,TL{:},TS{:}); title([cal.mnrStr{3},' [wt\%]'],TX{:},FS{:}); set(gca,'XTickLabel',[],'YTickLabel',[]);
    set(fh6,'CurrentAxes',ax(64));
    imagesc(Xc(2:end-1)./1e3,Zc(2:end-1)./1e3,x(2:end-1,2:end-1).*cx_mnr(2:end-1,2:end-1,4).*100); axis ij equal tight; box on; cb = colorbar;
    set(cb,TL{:},TS{:}); set(gca,TL{:},TS{:}); title([cal.mnrStr{4},' [wt\%]'],TX{:},FS{:}); set(gca,'XTickLabel',[]); ylabel('Depth [km]',TX{:},FS{:});
    set(fh6,'CurrentAxes',ax(65));
    imagesc(Xc(2:end-1)./1e3,Zc(2:end-1)./1e3,x(2:end-1,2:end-1).*cx_mnr(2:end-1,2:end-1,5).*100); axis ij equal tight; box on; cb = colorbar;
    set(cb,TL{:},TS{:}); set(gca,TL{:},TS{:}); title([cal.mnrStr{5},' [wt\%]'],TX{:},FS{:}); set(gca,'XTickLabel',[],'YTickLabel',[]);
    set(fh6,'CurrentAxes',ax(66));
    imagesc(Xc(2:end-1)./1e3,Zc(2:end-1)./1e3,x(2:end-1,2:end-1).*cx_mnr(2:end-1,2:end-1,6).*100); axis ij equal tight; box on; cb = colorbar;
    set(cb,TL{:},TS{:}); set(gca,TL{:},TS{:}); title([cal.mnrStr{6},' [wt\%]'],TX{:},FS{:}); set(gca,'XTickLabel',[],'YTickLabel',[]);
    set(fh6,'CurrentAxes',ax(67));
    imagesc(Xc(2:end-1)./1e3,Zc(2:end-1)./1e3,x(2:end-1,2:end-1).*cx_mnr(2:end-1,2:end-1,7).*100); axis ij equal tight; box on; cb = colorbar;
    set(cb,TL{:},TS{:}); set(gca,TL{:},TS{:}); title([cal.mnrStr{7},' [wt\%]'],TX{:},FS{:}); ylabel('Depth [km]',TX{:},FS{:});
    set(fh6,'CurrentAxes',ax(68));
    imagesc(Xc(2:end-1)./1e3,Zc(2:end-1)./1e3,x(2:end-1,2:end-1).*cx_mnr(2:end-1,2:end-1,8).*100); axis ij equal tight; box on; cb = colorbar;
    set(cb,TL{:},TS{:}); set(gca,TL{:},TS{:}); title([cal.mnrStr{8},' [wt\%]'],TX{:},FS{:}); set(gca,'YTickLabel',[]); xlabel('Width [km]',TX{:},FS{:});
    set(fh6,'CurrentAxes',ax(69));
    imagesc(Xc(2:end-1)./1e3,Zc(2:end-1)./1e3,x(2:end-1,2:end-1).*cx_mnr(2:end-1,2:end-1,9).*100); axis ij equal tight; box on; cb = colorbar;
    set(cb,TL{:},TS{:}); set(gca,TL{:},TS{:}); title([cal.mnrStr{9},' [wt\%]'],TX{:},FS{:}); set(gca,'YTickLabel',[]);
    sgtitle(['time = ',num2str(time/yr,3),' [yr]'],TX{:},FS{:},'Color','k');

    if plot_cv && iter > 0
        % plot residual fields in Fig. 7
        set(0,'CurrentFigure',fh7)
        set(fh7,'CurrentAxes',ax(61));
        imagesc(Xc(2:end-1)./1e3,Zc(2:end-1)./1e3,-res_W(:      ,2:end-1)); axis ij equal tight; box on; cb = colorbar;
        set(cb,TL{:},TS{:}); set(gca,TL{:},TS{:}); title(['res. $W$'],TX{:},FS{:}); ylabel('Depth [km]',TX{:},FS{:});
        set(fh7,'CurrentAxes',ax(62));
        imagesc(Xc(2:end-1)./1e3,Zc(2:end-1)./1e3, res_U(2:end-1,:      )); axis ij equal tight; box on; cb = colorbar;
        set(cb,TL{:},TS{:}); set(gca,TL{:},TS{:}); title(['res. $U$'],TX{:},FS{:}); set(gca,'YTickLabel',[]); xlabel('Width [km]',TX{:},FS{:});
        set(fh7,'CurrentAxes',ax(63));
        imagesc(Xc(2:end-1)./1e3,Zc(2:end-1)./1e3, res_P(2:end-1,2:end-1)); axis ij equal tight; box on; cb = colorbar;
        set(cb,TL{:},TS{:}); set(gca,TL{:},TS{:}); title(['res. $P$'],TX{:},FS{:}); set(gca,'YTickLabel',[]);
        sgtitle(['time = ',num2str(time/yr,3),' [yr]'],TX{:},FS{:},'Color','k');
    end

end

% plot model history
if plot_cv
    if ~exist('fh7','var'); fh8 = figure(VIS{:});
    else; set(0, 'CurrentFigure', fh8); clf;
    end
    subplot(3,1,1);
    plot(hist.time/yr,hist.DM./hist.sumM,'k-',LW{:}); hold on; axis tight; box on;
    ylabel('consv. $M$',TX{:},FS{:}); set(gca,TL{:},TS{:},'XTickLabel',[]);
    subplot(3,1,2);
    plot(hist.time/yr,hist.DS./hist.sumS,'k-',LW{:}); hold on; axis tight; box on;
    ylabel('consv. $S$',TX{:},FS{:}); set(gca,TL{:},TS{:},'XTickLabel',[]);
    subplot(3,1,3);
    plot(hist.time/yr,hist.DC./hist.sumC,'-',LW{:}); hold on; axis tight; box on;
    ylabel('consv. $C$',TX{:},FS{:}); set(gca,TL{:},TS{:});
    xlabel('Time [yr]',TX{:},FS{:});

    if ~exist('fh8','var'); fh9 = figure(VIS{:});
    else; set(0, 'CurrentFigure', fh9); clf;
    end
    subplot(3,1,1);
    plot(hist.time/yr,hist.EM,'k-',LW{:}); hold on; axis tight; box on;
    ylabel('error $M$',TX{:},FS{:}); set(gca,TL{:},TS{:},'XTickLabel',[]);
    subplot(3,1,2);
    plot(hist.time/yr,hist.ES,'k-',LW{:}); hold on; axis tight; box on;
    ylabel('error $S$',TX{:},FS{:}); set(gca,TL{:},TS{:},'XTickLabel',[]);
    subplot(3,1,3);
    plot(hist.time/yr,hist.EC,'-',LW{:}); hold on; axis tight; box on;
    ylabel('error $C$',TX{:},FS{:}); set(gca,TL{:},TS{:});
    xlabel('Time [yr]',TX{:},FS{:});
end

drawnow


% save output to file
if save_op && ~restart
    if Nx <= 10 && Nz <= 10  % print 0D plots
        name = [opdir,'/',runID,'/',runID,'_tch_',num2str(floor(step/nop))];
        print(fh1,name,'-dpng','-r300','-image');
        name = [opdir,'/',runID,'/',runID,'_aux_',num2str(floor(step/nop))];
        print(fh2,name,'-dpng','-r300','-image');
        name = [opdir,'/',runID,'/',runID,'_oxd_',num2str(floor(step/nop))];
        print(fh3,name,'-dpng','-r300','-image');
    elseif Nx <= 10  % create 1D plots
        name = [opdir,'/',runID,'/',runID,'_sol_',num2str(floor(step/nop))];
        print(fh1,name,'-dpng','-r300','-image');
        name = [opdir,'/',runID,'/',runID,'_aux_',num2str(floor(step/nop))];
        print(fh2,name,'-dpng','-r300','-image');
        name = [opdir,'/',runID,'/',runID,'_oxd_',num2str(floor(step/nop))];
        print(fh3,name,'-dpng','-r300','-image');
    else
        name = [opdir,'/',runID,'/',runID,'_sol_',num2str(floor(step/nop))];
        print(fh1,name,'-dpng','-r300','-image');
        name = [opdir,'/',runID,'/',runID,'_cmp_',num2str(floor(step/nop))];
        print(fh2,name,'-dpng','-r300','-image');
        name = [opdir,'/',runID,'/',runID,'_oxd_',num2str(floor(step/nop))];
        print(fh3,name,'-dpng','-r300','-image');
        name = [opdir,'/',runID,'/',runID,'_aux_',num2str(floor(step/nop))];
        print(fh4,name,'-dpng','-r300','-image');
        name = [opdir,'/',runID,'/',runID,'_gch',num2str(floor(step/nop))];
        print(fh5,name,'-dpng','-r300','-image');
        name = [opdir,'/',runID,'/',runID,'_mnr',num2str(floor(step/nop))];
        print(fh6,name,'-dpng','-r300','-image');
    end

    name = [opdir,'/',runID,'/',runID,'_',num2str(floor(step/nop))];
    save(name,'U','W','P','Pt','x','m','chi','mu','X','S','C','T','c','cm','cx','TE','IR','te','ir','dSdt','dCdt','dXdt','dTEdt','dIRdt','Gx','rho','eta','eII','tII','dt','time','step','VolSrc','wx','wm');
    name = [opdir,'/',runID,'/',runID,'_cont'];
    save(name,'U','W','P','Pt','x','m','chi','mu','X','S','C','T','c','cm','cx','TE','IR','te','ir','dSdt','dCdt','dXdt','dTEdt','dIRdt','Gx','rho','eta','eII','tII','dt','time','step','VolSrc','wx','wm');
    name = [opdir,'/',runID,'/',runID,'_hist'];
    save(name,'hist');
    
    if step == 0
        logfile = [opdir,'/',runID,'/',runID,'.log'];
        if exist(logfile,'file'); delete(logfile); end
        diary(logfile)
    end
end