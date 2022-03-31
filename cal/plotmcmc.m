function [xMAP] = plotmcmc (x, P, x0, xbnds, count, BurnIn, VarNames)
%
% [xMAP] = PlotMCMCAnalytics (x, P, x0, xbnds, count, BurnIn, VarNames)
% plots the outputs of mcmc inversion
%
% INPUTS
% x         matrix of model parameters from mcmc [Niter x Nvar]
% P         vector of posterior probabilities [Niter x 1]
% x0        starting model [1 x Nvar, can be empty]
% xbnds     bnds of parameters [Nvar x 2]
% count     number of accepted models [scalar]
% BurnIn    burn in of markov chain [scalar]
% VarNames  name of variables
%
% OUTPUT
% xMAP      maximum a posteriori model [1 x Nvar]
%


[Niter, Nvar] = size(x);

Nrow = floor(sqrt(Nvar));
Ncol = ceil(Nvar/Nrow);

% find the maximum a posteriori model
[~, xMAPind] = max(P);
xMAP = x(xMAPind,:);

VarVary = diff(xbnds,[],2)>0;

if isempty(x0), x0 = x(1,:); end

figure(20);

% plot distributions of the model parameters
for mi = 1:Nvar
    subplot(Nrow,Ncol,mi);
    
    % plot the distribution
    histogram(x(BurnIn:end,mi), 40, 'EdgeColor', 'none');   hold on;
    
    % bounds
    plot(xbnds(mi,1)*ones(1,2), ylim, 'r-');
    plot(xbnds(mi,2)*ones(1,2), ylim, 'r-');
    
    % starting model
    plot(x0(mi)*ones(1,2), ylim, 'k:');
    
    % maximum a posteriori model
    plot(xMAP(mi)*ones(1,2), ylim, 'r:');
    hold off;
    
    if VarVary(mi); xlim(xbnds(mi,:)); end
    title(VarNames{mi});
end
sgtitle('Posterior PDFs (black: starting model; red: best-fit)');


% plot the chain to see if well mixed
figure(21); 
semilogy(-P);
ylimits = ylim;
hold on; plot(BurnIn*ones(1,2), ylimits, 'Color', 0.8*ones(1,3)); hold off;
text(BurnIn, ylimits(2), 'burn in', 'FontSize', 16, 'Rotation', 90, 'VerticalAlignment','bottom');
set(gca,'ydir','reverse');
yticklabs = get(gca,'YTickLabel');
set(gca,'YTickLabel',strcat('-',yticklabs));
xlabel('iteration number'); ylabel('log likelihood');
title(['Acceptance ratio = ' num2str(100*count/Niter,4)]);

end