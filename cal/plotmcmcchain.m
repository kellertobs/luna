function [] = plotmcmcchain (P, BurnIn, count)

Niter = length(P);

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