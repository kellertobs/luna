clear all; close all;

addpath('../src');

% set run parameters
runID    =  '1D_luna_ref';       % run identifier
opdir    =  '../out';            % output directory
restart  =  0;                   % restart from file (0: new run; <1: restart from last; >1: restart from specified frame)
nop      =  500;                 % output frame plotted/saved every 'nop' time steps
plot_op  =  1;                   % switch on to live plot of results
save_op  =  1;                   % switch on to save output to file
plot_cv  =  0;                   % switch on to live plot iterative convergence

% set model domain parameters
D        =  1000e3;              % chamber depth [m]
L        =  1000e3/500;          % chamber width [m]
N        =  500 + 2;             % number of grid points in z-direction (incl. 2 ghosts)
h        =  D/(N-2);             % grid spacing (equal in both dimensions, do not set) [m]

% set model timing parameters
Nt       =  1e6;                 % number of time steps to take
hr       =  3600;                % conversion seconds to hours
yr       =  24*365.25*hr;        % conversion seconds to years
tend     =  1e3*yr;              % end time for simulation [s]
dt       =  1e4;                 % initial time step [s]
dtmax    =  1e6;                 % maximum time step [s]

% set initial thermo-chemical state
seed     =  15;                  % random perturbation seed
smth     =  (N/30)^2;            % regularisation of initial random perturbation
zlay     =  0.5;                 % layer thickness (relative to domain depth D)
wlay_T   =  4*h/D;               % thickness of smooth layer boundary (relative to domain depth D)
wlay_c   =  2*h/D;               % thickness of smooth layer boundary (relative to domain depth D)
T0       =  1800;                % temperature top layer [deg C]
T1       =  1800;                % temperature base layer [deg C]
dT       =  0;                   % amplitude of random noise [deg C]
c0       =  [0.44,0.10,0.20,0.13,0.12,0.01]; % major component top layer
c1       =  [0.44,0.10,0.20,0.13,0.12,0.01]; % major component base layer
dc       =  [1,-1,1,-1,0,0].*0e-4; % amplitude of random noise [wt SiO2]

% set model trace and isotope geochemistry parameters
te0      =  [1,1,1,1];           % trace elements top layer [wt ppm]
te1      =  [1,1,1,1];           % trace elements base layer [wt ppm]
dte      =  0e-3.*[-1,-1,1,1];   % trace elements random noise [wt ppm]
Kte      =  [0.01,0.1,3,10];     % trace elements partition coefficients
ir0      =  [0,-1];              % isotope ratios top layer [delta]
ir1      =  [0, 1];              % isotope ratios base layer [delta]
dir      =  [0, 0];              % isotope ratios random noise [delta]

% set thermo-chemical boundary parameters
Ptop     =  1e4;                 % top pressure [Pa]
bndmode  =  3;                   % boundary assimilation mode (0 = none; 1 = top only; 2 = bot only; 3 = top/bot;)
bndinit  =  0;                   % switch on (1) to initialise with already established boundary layers
dw       =  h;                   % boundary layer thickness for assimilation [m]
tau_T    =  5*yr;                % wall cooling/assimilation time [s]
Ttop     =  0;                   % wall temperature [degC] (nan = insulating)
Tbot     =  1950;                % wall temperature [degC] (nan = insulating)

% set thermo-chemical material parameters
cP       =  1200;                % heat capacity [J/kg/K]
kT0      =  4;                   % thermal conductivity [W/m/K]
calID    = 'luna';               % calibration ID

% set model rheology parameters
etax0    =  1e16;                % solid reference viscosity [Pas]
AA       = [ 0.6907, 0.1853; 0.1311, 0.1644; ];  % permission slopes
BB       = [ 0.6904, 0.3096; 0.9988, 0.0012; ];  % permission step locations
CC       = [ 0.5145, 0.1831; 0.6808, 1.8541; ];  % permission step widths

% set model buoyancy parameters
rhox0    =  [3270,4390,3000,3250,2730,2620];  % crystal phase ref. density [kg/m3] (at T0,cphs0,Ptop)
rhom0    =  [2710,3580,2580,2850,2530,2310];  % melt phase ref. density [kg/m3] (at T0,cphs0,Ptop)
aT       =  5e-5;                % thermal expansivity [1/K]
bPx      =  1e-11;               % solid compressibility [1/Pa]
bPm      =  3e-11;               % melt compressibility [1/Pa]
d0       =  1e-3;                % crystal size [m]
g0       =  1.62;                % gravity [m/s2]

% set numerical model parameters
TINT     =  'bd3i';              % time integration scheme ('bwei','cnsi','bd3i','bd3s')
ADVN     =  'weno5';             % advection scheme ('centr','upw1','quick','fromm','weno3','weno5','tvdim')
CFL      =  1.00;                % (physical) time stepping courant number (multiplies stable step) [0,1]
rtol     =  1e-4;                % outer its relative tolerance
atol     =  1e-7;                % outer its absolute tolerance
maxit    =  50;                  % maximum outer its
lambda   =  0.50;                % iterative lag parameter equilibration
etareg   =  1e11;                % regularisation factor for viscosity resisting convection
sgrreg   =  1e0;                 % regularisation factor for viscosity resisting segregation


%*****  RUN NAKHLA MODEL  *************************************************
run('../src/main')
%**************************************************************************

