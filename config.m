%% project: Zavorotny-Voronovich_model
%
%
% Configuration script - Simulation parameters for the Zavorotny-Voronovich
% model.

%% NECESSARY CONSTANTS
global C_light LAMBDA R0
C_light = 3e8;        % Speed of light
fL1 = 1575.42e6  ;    % L1 carrier frequency
LAMBDA = C_light/fL1; % Wavelength
R0 = 6378137;         % Mean Earth radius
kB = 1.38e-23;        % Boltzmann constant
T0 = 290;             % Reference temperature
KzdB = -123.6;
Kz = 10^(KzdB*.1);    % Minimum power constant: PtGt/(20200e3)^2

%% OCEAN MODEL
U10 = 3;               % Wind speed 10 meters above ocean surface
AGE = 1  ;             % Wave age: fully developed (1), mature (2), young (3)
kco = 2*pi/(3*LAMBDA); % Cut-off wavenumber
phi0 = 0;              % Wind direction relative to XZ plane

% Slopes statistical distribution: 'Gaussian' or 'Delta' for specular
% reflection.
dist = 'Gaussian';     
% dist = 'Delta';

%% WOODWARD AMBIGUITY FUNCTION

dt = 1/8; % Delay resolution
taumin = -3; % Minimum delay value
taumax = 3; % Maximum delay value
Ti = 1e-3; % Coherent integration time
DF = 1/Ti;
df = DF/4; % Doppler resolution
fmin = -5e3; % Minimum Doppler value
fmax = 5e3; % Maximum Doppler value

% Coordinate of the maximum (always at 0,0 to generate DDM)
ret = 0;
fd = 0;

%% GEOMETRY
e_deg = 70;   % Elevation angle [deg]
phit_deg = 0; % Angle between moving direction of the GPS satellite and the XZ plane
phir_deg = 0; % Angle between moving direction of the LEO satellite and the XZ plane
hr = 630e3;   % LEO altitude

%% SPATIAL MAP (ocean surface grid)

Xmax = 100e3;                 % Maximum X axis value (symmetrical) [m]
Ymax = 100e3;                 % Maximum Y axis value (symmetrical) [m]
Nx = 1e2;                    % Number of points in X axis
Ny = 1e2;                    % Number of points in Y axis
x = linspace(-Xmax,Xmax,Nx); % X axis
y = linspace(-Ymax,Ymax,Ny); % Y axis
dx = x(2)-x(1);              % X axis resolution [m]
dy = y(2)-y(1);              % Y axis resolution [m]
dA = dx*dy;                  % surface differential

%% RECEIVER
FndB = 3;          % Rx total noise figure [dB].
Fn = 10^(FndB*.1);
Ta = 200;          % Antenna temperature
T = (Fn-1)*T0+Ta;
GrdB = 7;          % Receiver antenna gain [dB]
Gr = 10^(GrdB*.1);

%% DDM
% Selection: DDM or WF
DDMWFflag = 'DDM';
% DDMWFflag = 'WF';

switch DDMWFflag
    case 'DDM'
        Fmax = 3e3;          % DDM Doppler axis maimum value (symmetrical)
        dff = DF/2;          % DDM Doppler axis resolution
        ff = -Fmax:dff:Fmax; % DDM Doppler axis
    case 'WF'
        ff = 0;              % Only 1 bin, WF
end

dtt = 1/4;         % DDM delay axis resolution
ttaum = -2:dtt:20; % DDM delay axis
