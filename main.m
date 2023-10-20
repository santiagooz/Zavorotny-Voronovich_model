%% project: Zavorotny-Voronovich_model
% Santiago Ozafrain, Aug 2015. Updated: Oct 2023.
% UIDET-SENyT, Facultad de Ingenieria, UNLP.
%
% Main script - Generation of delay-Doppler maps (DDM) or correlation
% waveforms (WF) for ocean reflected GPS C/A signals according to the
% Zavorotny and Voronovich model from V. U. Zavorotny and A. G. Voronovich,
% "Scattering of GPS signals from the ocean with wind remote sensing
% application," in IEEE Transactions on Geoscience and Remote Sensing, vol.
% 38, no. 2, pp. 951-964, March 2000, doi: 10.1109/36.841977.

clear;clc;close all;

config % Loads simulation parameters

% Woodward Ambiguity Function generation
[WAF,Rcm,S,taum,f] = WAFgen(ret,fd,dt,taumin,taumax,Ti,df,fmin,fmax);

% Glistening Zone generation
[ Dtdif,fddif,sigma0,Xsp,Ysp,Rt,Rr ] = GZgen(x,y,dx,dy,e_deg,phit_deg,hr,phir_deg,U10,AGE,kco,phi0,dist);

%% DDM/WF generation
% Computation of the integrand in equation (27) in ZV_2000 paper.

INT=zeros(length(ff),length(ttaum));

wb = waitbar(0, '   Generating DDM/WF: 0% completed    ');
tic
for kk = 1:length(ff)
    for ll = 1:length(ttaum)
        for ii = 1:length(x)
            for jj = 1:length(y)
                r = [x(ii);y(jj);0];            % r point in space
                dT = ttaum(ll)-Dtdif(jj,ii);    % Difference between local replica delay and the corresponding to the point r
                dF = ff(kk)-fddif(jj,ii);       % Difference between local replica Doppler and the corresponding to the point r
                [~,IndT] = min(abs(dT-taum));
                [~,IndF] = min(abs(dF-f));
                Rts = (r-Rt);                   % Distance between GPS satellite and the r point
                Rsr = (Rr-r);                   % Distance between LEO satellite and the r point
                INT(kk,ll) = sigma0(jj,ii)*abs(WAF(IndF,IndT))^2*dA/norm(Rsr)^2+INT(kk,ll);
            end
        end

        % progression bar
        time_past = toc;
        time_per_loop = time_past/((kk-1)*length(ttaum)+ll);
        time_left = (length(ff)*length(ttaum)-(kk-1)*length(ttaum)-ll)*time_per_loop;
        time_left_h = floor(time_left/3600);
        time_left_m = floor((time_left - time_left_h*3600)/60);
        time_left_s = floor(time_left - time_left_h*3600 - time_left_m*60);
        msg = sprintf('   Generating DDM/WF: %i%% completed   \n%i:%i:%i remaining', floor(((kk-1)*length(ttaum)+ll)/(length(ff)*length(ttaum))*100), time_left_h, time_left_m, time_left_s);
        waitbar(((kk-1)*length(ttaum)+ll)/(length(ff)*length(ttaum)), wb, msg)
    end
end
close(wb)

% DDM scaled as SNR
SNR = Kz/(4*pi)*(LAMBDA/(4*pi))^2*Ti/(kB*T)*INT;
SNRdB = 10*log10(SNR);


%% Plot the generated DDM/WF

FontSize=12;
switch DDMWFflag
    case 'DDM'
        figure;surf(ttaum,ff*1e-3,SNR/max(SNR(:))),title('normalized DDM - ocean reflected signal','FontSize',FontSize),xlabel('\tau [chips]','FontSize',FontSize),ylabel('f [kHz]','FontSize',FontSize),% shading interp
        view([45 -45 45])
        figure;plot(ttaum,SNR/max(SNR(:))),title('normalized WF - ocean reflected signal','FontSize',FontSize),xlabel('\tau [chips]','FontSize',FontSize);axis([ttaum(1) 5 0 1.1]);grid on% shading interp

    case 'WF'
        figure;plot(ttaum,SNR/max(SNR)),title('normalized WF - ocean reflected signal'),xlabel('\tau [chips]')
end
set(gca,'FontSize', FontSize)

%% Plot iso-delay and iso-Doppler curves, and the BRCS in the spatial map

LineWidth=2;
FontSize=18;
figure;contour(x*1e-3,y*1e-3,Dtdif,min(Dtdif(:)):max(Dtdif(:)),'LineWidth',LineWidth);
h = colorbar;
ylabel(h, '\Delta\tau [chips]','FontSize',FontSize)
xlabel('x [km]','FontSize',FontSize)
ylabel('y [km]','FontSize',FontSize)
title('Iso-delay lines','FontSize',FontSize)
% set(gca,'FontSize',FontSize)

figure;contour(x*1e-3,y*1e-3,fddif*1e-3,min(fddif(:))*1e-3:max(fddif(:))*1e-3,'LineWidth',LineWidth);
h = colorbar;
ylabel(h, '\Delta f_D [kHz]','FontSize',FontSize)
xlabel('x [km]','FontSize',FontSize)
ylabel('y [km]','FontSize',FontSize)
title('Iso-Doppler lines','FontSize',FontSize)
% set(gca,'FontSize',FontSize)

sigma0=sigma0/max(sigma0(:));
figure;contour(x*1e-3,y*1e-3,sigma0,min(sigma0(:)):(1-min(sigma0(:)))/10:1,'LineWidth',LineWidth);
h = colorbar;
ylabel(h, '\sigma_0','FontSize',FontSize)
xlabel('x [km]','FontSize',FontSize)
ylabel('y [km]','FontSize',FontSize)
title('Normalized Bistatic Radar Cross Section','FontSize',FontSize)
% set(gca,'FontSize',FontSize)