function [xMAP] = plotpdfs (x, P, VarNames, xbnds, x0, BurnIn)
%
% [xMAP] = plotpdfs (x, P, xbnds, x0, VarNames, BurnIn)
% plots the posterior distributions
%
% INPUTS
% x         matrix of model parameters [Niter x Nvar]
% P         vector of posterior probabilities [Niter x 1]
% VarNames  name of variables [cell, Nvar x 1]
% xbnds     bnds of parameters [Nvar x 2]
% x0        starting model [1 x Nvar, can be empty]
% BurnIn    burn in of markov chain, for mcmc distribs [scalar]
%
% OUTPUT
% xMAP      maximum a posteriori model [1 x Nvar]
%

[~, Nvar] = size(x);

% check arguments
if nargin<3, VarNames = strcat({'v'},num2str((1:Nvar)','%02d'));    end
if nargin<4,    xbnds = [-inf,inf]*ones(Nvar,1);                    end
if nargin<5,       x0 = [];                                         end
if nargin<6,   BurnIn = 1;                                          end

Nrow = floor(sqrt(Nvar));
Ncol = ceil(Nvar/Nrow);

% find the maximum a posteriori model
[~, xMAPind] = max(P(BurnIn:end));
xMAP = x(BurnIn + xMAPind,:);

VarVary = diff(xbnds,[],2)>0;

if isempty(x0), x0 = nan(1,Nvar); end

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

end