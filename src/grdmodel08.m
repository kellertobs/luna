% MATLAB script to compute silicate melt viscosity.
%
%  Citation: Giordano D, Russell JK, & Dingwell DB (2008)
%  Viscosity of Magmatic Liquids: A Model. Earth & Planetary Science
%  Letters, v. 271, 123-134.
%
% ________________________________________________________________
% INPUT: Chemical compositions of silicate melts (wt% oxides) as:
% SiO2 TiO2 Al2O3 FeO(T) MnO MgO CaO Na2O K2O P2O5 H2O F2O-1
% One line for each melt composition
% _________________________________________________________________

% ________________________________________________________________
% OUTPUT: eta values at T(C) values (1 line per melt)
% ________________________________________________________________

% VFT Multicomponent-Model Coedfficients
% -----------------------------------------------------------------

% modified and optimised by Tobias Keller, 10. June, 2022

function   [eta] = grdmodel08(wtm,TC)

AT  =  -4.55;
bb  =  [159.56  -173.34 72.13 75.69 -38.98 -84.08 141.54 -2.43 -0.91 17.62];
cc  =  [2.75 15.72 8.32 10.2 -12.29 -99.54 0.3 ];

% Function molefrac_grd: converts wt % oxide basis to mole % oxide basis
[~,xmf_t]  =  molefrac(wtm);

% Load composition-basis matrix for multiplication against model-coefficients
% Result is two matrices bcf[nx by 10] and ccf[nx by 7]
siti  =  xmf_t(:,1) + xmf_t(:,2);
tial  =  xmf_t(:,2)+xmf_t(:,3);
fmm   =  xmf_t(:,4) + xmf_t(:,5) + xmf_t(:,6);
nak   =  xmf_t(:,8) + xmf_t(:,9);
b1    =  siti;
b2    =  xmf_t(:,3);
b3    =  xmf_t(:,4) + xmf_t(:,5) + xmf_t(:,10);
b4    =  xmf_t(:,6);
b5    =  xmf_t(:,7);
b6    =  xmf_t(:,8) + xmf_t(:,11) + xmf_t(:,12);
b7    =  xmf_t(:,11) + xmf_t(:,12) + log(1+xmf_t(:,11));
b12   =  siti.*fmm;
b13   =  (siti + xmf_t(:,3) + xmf_t(:,10)).*( nak + xmf_t(:,11) );
b14   =  xmf_t(:,3).*nak;

c1    =  xmf_t(:,1);
c2    =  tial;
c3    =  fmm;
c4    =  xmf_t(:,7);
c5    =  nak;
c6    =  log(1+xmf_t(:,11) + xmf_t(:,12));
c11   =  xmf_t(:,3) + fmm + xmf_t(:,7) - xmf_t(:,10);
c11   =  c11.*(nak + xmf_t(:,11) + xmf_t(:,12));
bcf   =  [b1 b2 b3 b4 b5 b6 b7 b12 b13 b14];
ccf   =  [c1 c2 c3 c4 c5 c6 c11];

BT    =  sum(bb.*bcf,2);
CT    =  sum(cc.*ccf,2);

TK    =  TC + 273.15;
eta   =  exp(AT + BT./(TK(:)-CT));

