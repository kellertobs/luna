clear all; close all;

addpath('../src');

% set run parameters
runID    =  '1D_luna_6cmp_S';    % run identifier
opdir    =  '../out/';           % output directory
restart  =  0;                   % restart from file (0: new run; <1: restart from last; >1: restart from specified frame)
nop      =  100;                 % output frame plotted/saved every 'nop' time steps
plot_op  =  1;                   % switch on to live plot of results
save_op  =  1;                   % switch on to save output to file
plot_cv  =  0;                   % switch on to live plot iterative convergence
react    =  1;                   % switch on reactive mode
diseq    =  1;                   % switch on disequilibrium approach
entr_mth =  1;                   % switch on to use entropy equation for heat evolution, else temperature equation used
bnchm    =  0;                   % switch on to run manufactured solution benchmark on fluid mechanics solver

% set model domain parameters
D        =  1000e3;              % chamber depth [m]
L        =  1000e3/200;          % chamber width [m]
N        =  200 + 2;             % number of grid points in z-direction (incl. 2 ghosts)
h        =  D/(N-2);             % grid spacing (equal in both dimensions, do not set) [m]

% set model timing parameters
M        =  1e5;                 % number of time steps to take
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
T0       =  1675;                % temperature top layer [deg C]
T1       =  1675;                % temperature base layer [deg C]
dT       =  0;                   % amplitude of random noise [deg C]
c0       =  [0.40,0.10,0.21,0.13,0.12,0.04]; % major component top layer [liquid fraction from catmip16 fig8]
cl       =  [0.40,0.10,0.21,0.13,0.12,0.04]; % major component base layer [liquid fraction from catmip16 fig8]
dc       =  [1,-1,1,-1,0,0].*0e-4; % amplitude of random noise [wt SiO2]

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
dsi      =  0.00;                % stable isotope ratio random noise [delta]
ri0      =  1;                   % radiogenic isotope top layer [wt ppm]
ri1      =  1;                   % radiogenic isotope base layer [wt ppm]
dri      = -0.00;                % radiogenic isotope random noise [wt ppm]
KRIP     =  10.;                 % radiogenic parent isotope partition coefficient
KRID     =  0.1;                 % radiogenic daughter isotope partition coefficient
HLRIP    =  1e3*yr;              % radiogenic parent isotope half-life [s]
HLRID    =  1e2*yr;              % radiogenic daughter isotope half-life [s]

% set thermo-chemical boundary parameters
Ptop     =  1e5;                 % top pressure [Pa]
bndmode  =  3;                   % boundary assimilation mode (0 = none; 1 = top only; 2 = bot only; 3 = top/bot;)
bndinit  =  0;                   % switch on (1) to initialise with already established boundary layers
dw       =  1*h;                 % boundary layer thickness for assimilation [m]
tau_T    =  3*yr;                % wall cooling/assimilation time [s]
Ttop     =  100;                 % wall temperature [degC] (nan = insulating)
Tbot     =  1800;                % wall temperature [degC] (nan = insulating)

% set thermo-chemical material parameters
cP       =  1300;                % heat capacity [J/kg/K]
kc0      =  1e-6;                % chemical diffusivity [kg/m/s]
kT0      =  1;                   % thermal conductivity [W/m/K]

% set phase diagram parameters
calID    = 'luna6';
tau_r    =  60;                  % reaction time [s]
Dsx      = -350;                 % entropy change of crystallisation [J/kg]

% set model rheology parameters
etam0    =  1e1;                 % melt  reference viscosity [Pas]
etax0    =  1e16;                % solid reference viscosity [Pas]
AA       = [ 0.6907, 0.1853; 0.1311, 0.1644; ];  % permission slopes
BB       = [ 0.6904, 0.3096; 0.9988, 0.0012; ];  % permission step locations
CC       = [ 0.5145, 0.1831; 0.6808, 1.8541; ];  % permission step widths

% set model buoyancy parameters
rhox0    =  [3270,4390,3500,3250,2730,2620];  % crystal phase ref. density [kg/m3] (at T0,cphs0,Ptop)
rhom0    =  rhox0 - 300;         % melt phase ref. density [kg/m3] (at T0,cphs0,Ptop)
aT       =  3e-5;                % thermal expansivity [1/K]
bPx      =  5e-12;               % solid compressibility [1/Pa]
bPm      =  1e-11;               % melt compressibility [1/Pa]
dx       =  1e-3;                % crystal size [m]
g0       =  1.62;                % gravity [m/s2]

% set numerical model parameters
CFL      =  0.1;                 % (physical) time stepping courant number (multiplies stable step) [0,1]
ADVN     =  'FRM';               % advection scheme ('UPW2', 'UPW3', or 'FRM')
theta    =  0.5;                 % time-stepping parameter (1 = 1st-order implicit; 1/2 = 2nd-order semi-implicit)
rtol     =  1e-4;                % outer its relative tolerance
atol     =  1e-7;                % outer its absolute tolerance
maxit    =  10;                  % maximum outer its
alpha    =  0.25;                % iterative lag parameter equilibration
delta    =  5;                   % smoothness of segregation speed
etareg   =  1e3;                 % bounds on viscosity resisting convection for regularisation/stabilisation
sgrreg   =  1e0;                 % bounds on viscosity resisting convection for regularisation/stabilisation
kcreg    =  1e2;                 % chemical diffusivity for regularisation [kg/m/s]
kTreg    =  1e2;                 % thermal conductivity for regularisation [W/m/K]

% create output directory
if ~isfolder([opdir,'/',runID])
    mkdir([opdir,'/',runID]);
end

% save input parameters and runtime options (unless restarting)
if restart == 0 
    parfile = [opdir,'/',runID,'/',runID,'_par'];
    save(parfile);
end

% run code
run('../src/main')
