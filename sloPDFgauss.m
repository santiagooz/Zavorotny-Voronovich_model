%% project: Zavorotny-Voronovich_model
%
%
% sloPDFgauss - function the generates the slopes probability density
% function for the ocean surface affected by wind (gaussian)

function [P] = sloPDFgauss(s,C,phi0)

Rot = [cos(phi0) -sin(phi0);sin(phi0) cos(phi0)];
M = Rot*C*Rot';
P = 1/(2*pi*sqrt(det(M)))*exp(-1/2*s'/M*s);

end
