%% project: Zavorotny-Voronovich_model
%
%
% VelSat - calculates the speed of a satellite in orbit around the Earth
% with an altitude h

function [ vs,vr ] = VelSat( h )

global R0 % Mean Earth radius

G = 6.67384e-11; % Universal gravitational constant
M = 5.972e24;    % Earth mass

vs = sqrt(G*M./(h+R0)); % Absolute satellite speed
vT = 2*pi*R0/24/3600;   % Absolute speed of the Earth surface
vr = vs-vT;             % Satellite speed relative to the ocean surface
end