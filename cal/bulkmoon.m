function [c] = bulkmoon (model)
% 
% bulk composition of magma ocean for the moon
% taken from Table 1 of Schmidt & Kraettli 2022 JGR Planets
% 
% order:
% {'SiO$_2$','Al$_2$O$_3$','FeO','MgO','CaO','Na$_2$O'};
% 
% output is oxide composition of bulk moon, in wt%
% default uses Taylor whole moon composition

if nargin==0, model = "TWM"; end

switch model
    case "TWM"   % Taylor Whole Moon composition
        c = [44.5, 6.14, 10.9, 32.7, 4.60, 0.09];
    
    case "LPUM"  % lunar primitive upper mantle
        c = [46.1, 3.93, 7.62, 38.3, 3.18, 0.05];
end

c = c./sum(c);

end