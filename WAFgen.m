%% project: Zavorotny-Voronovich_model
%
%
% WAFgen - function the generates the Wodward ambiguity function: WAF

function [ WAF,Rcm, S, taum, f ] = WAFgen( ret, fd, dt, taumin, taumax, Ti, df, fmin, fmax)

% Generation of the autocorrelation function of GPS Gold codes: Rc.
c1 = cacode(1);
c2 = cacode(1);
c2 = circshift(c2,[0 -ret]);
C1 = kron(c1,ones(1,1/dt));
C2 = kron(c2,ones(1,1/dt));
[Rc,tau] = xcorr(C1,C2);
Rc = Rc/length(C1);
tau = tau*dt;

% We keep the values for the input delay interval
at = find(tau==taumin);
bt = find(tau==taumax);
Rcm = Rc(at:bt);
taum = tau(at:bt);

% WAF in the Doppler axis
f = fmin:df:fmax;
S = sinc(f*Ti).*exp(-1i*pi*f*Ti);
ind0 = find(f==0);
indfd = find(f==fd);
deltafd = indfd-ind0;
S = circshift(S,[0 deltafd]);

WAF= S.'*Rcm; % Woodward ambiguity function

end

