%% project: Zavorotny-Voronovich_model
%
%
% sloPDFdelta - function the generates the slopes probability density
% function for the specular reflection case (surface ideally flat).

function [P] = sloPDFdelta(s,dx,dy,Smin)

P = 1/dx/dy*(s(1)==Smin(1))*(s(2)==Smin(2));

end
