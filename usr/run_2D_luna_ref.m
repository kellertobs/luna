clear all; close all;

addpath('../src');

% set run parameters
runID    =  '2D_luna_ref_eta12';       % run identifier
opdir    =  '../out';            % output directory
restart  =  0;                   % restart from file (0: new run; <1: restart from last; >1: restart from specified frame)
nop      =  200;                 % output frame plotted/saved every 'nop' time steps
plot_op  =  1;                   % switch on to live plot of results
save_op  =  1;                   % switch on to save output to file
plot_cv  =  0;                   % switch on to live plot iterative convergence
diseq    =  1;                   % switch on disequilibrium approach
entr_mth =  1;                   % switch on to use entropy equation for heat evolution, else temperature equation used
bnchm    =  0;                   % switch on to run manufactured solution benchmark on fluid mechanics solver

% set model domain parameters
D        =  1000e3;              % chamber depth [m]
L        =  1000e3;              % chamber width [m]
N        =  200 + 2;             % number of grid points in z-direction (incl. 2 ghosts)
h        =  D/(N-2);             % grid spacing (equal in both dimensions, do not set) [m]

% set model timing parameters
M        =  1e6;                 % number of time steps to take
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
T0       =  1700;                % temperature top layer [deg C]
T1       =  1700;                % temperature base layer [deg C]
dT       =  0;                   % amplitude of random noise [deg C]
c0       =  [0.40,0.10,0.21,0.13,0.12,0.04]; % major component top layer
cl       =  [0.40,0.10,0.21,0.13,0.12,0.04]; % major component base layer
dc       =  [1,-1,1,-1,1,-1].*1e-3; % amplitude of random noise [wt SiO2]

% set model trace and isotope geochemistry parameters
it0      =  1;                   % incompatible tracer top layer [wt ppm]
it1      =  1;                   % incompatible tracer base layer [wt ppm]
dit      =  0.00;                % incompatible tracer random noise [wt ppm]
KIT      =  1e-2;                % incompatible tracer partition coefficient
ct0      =  1;                   % compatible tracer top layer [wt ppm]
ct1      =  1;                   % compatible tracer base layer [wt ppm]
dct      = -0.00;                % compatible tracer random noise [wt ppm]
KCT      =  1e+2;                % compatible tracer partition coefficient
si0      =  0;                   % stable isotope ratio top layer [delta]
si1      =  0;                   % stable isotope ratio base layer [delta]
dsi      =  1.00;                % stable isotope ratio random noise [delta]
ri0      =  1;                   % radiogenic isotope top layer [wt ppm]
ri1      =  1;                   % radiogenic isotope base layer [wt ppm]
dri      = -0.00;                % radiogenic isotope random noise [wt ppm]
KRIP     =  10.;                 % radiogenic parent isotope partition coefficient
KRID     =  0.1;                 % radiogenic daughter isotope partition coefficient
HLRIP    =  1e3*yr;              % radiogenic parent isotope half-life [s]
HLRID    =  1e2*yr;              % radiogenic daughter isotope half-life [s]

% set thermo-chemical boundary parameters
Ptop     =  1e4;                 % top pressure [Pa]
bndmode  =  3;                   % boundary assimilation mode (0 = none; 1 = top only; 2 = bot only; 3 = top/bot;)
bndinit  =  0;                   % switch on (1) to initialise with already established boundary layers
dw       =  1*h;                 % boundary layer thickness for assimilation [m]
tau_T    =  5*yr;                % wall cooling/assimilation time [s]
Ttop     =  0;                   % wall temperature [degC] (nan = insulating)
Tbot     =  1900;                % wall temperature [degC] (nan = insulating)

% set thermo-chemical material parameters
cP       =  1300;                % heat capacity [J/kg/K]
kT0      =  4;                   % thermal conductivity [W/m/K]
calID    = 'luna6';              % calibration ID
Dsx      = -300;                 % entropy change of crystallisation [J/kg]

% set model rheology parameters
etam0    =  1e1;                 % melt  reference viscosity [Pas]
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
dx       =  1e-3;                % crystal size [m]
g0       =  1.62;                % gravity [m/s2]

% set numerical model parameters
CFL      =  0.75;                % (physical) time stepping courant number (multiplies stable step) [0,1]
ADVN     =  'FRM';               % advection scheme ('UPW2', 'UPW3', or 'FRM')
rtol     =  1e-3;                % outer its relative tolerance
atol     =  1e-6;                % outer its absolute tolerance
maxit    =  10;                  % maximum outer its
lambda   =  0.5;                 % iterative lag parameter equilibration
etareg   =  1e12;                % regularisation factor for viscosity resisting convection
sgrreg   =  1e0;                 % regularisation factor for viscosity resisting segregation
dffreg   =  etareg^0.5;          % regularisation factor for thermal, chemical, phase diffusion


%*****  RUN NAKHLA MODEL  *************************************************
run('../src/main')
%**************************************************************************

