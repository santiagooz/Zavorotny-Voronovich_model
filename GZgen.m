%% project: Zavorotny-Voronovich_model
%
%
% GZgen - function the generates the distribution of delay, Doppler and
% bistatic radar cross-section values over the glistening zone.

function [ Dtdif, fddif, sigma0, Xsp, Ysp, Rt, Rr ] = GZgen(x, y, dx, dy, e_deg, phit_deg, hr, phir_deg, U10, AGE, kco, phi0, distribution)

global LAMBDA C_light R0

e_water = 73;    % Constante dielectrica del agua marina.
K = 2*pi/LAMBDA; % Wavenumber.

% GPS sat
ht = 20200e3;     % GPS satellite altitude
Vt0 = VelSat(ht); % GPS satellite speed
phit = phit_deg*pi/180;

% LEO sat
Vr0 = VelSat(hr); % LEO satellite speed
phir=phir_deg*pi/180;

e = e_deg*pi/180;


CENTER = [0;0;-R0]; % Position of the center of the spheric Earth
Rt = zeros(3, 1); % GPS satellite position
Rr = zeros(3, 1); % LEO satellite position

while abs(norm(Rt-CENTER)-ht-R0)>1e2
    Rt = Rt + 100*[-cos(e);0;sin(e)];
end
while abs(norm(Rr-CENTER)-hr-R0)>1e2
    Rr = Rr + 100*[cos(e);0;sin(e)];
end

Rrx = Rr(1);
Rrz = Rr(3);
Vrx = Vr0/sqrt(1+(Rrx/Rrz)^2);
Vrz = -Rrx/Rrz*Vrx;
Vr = [Vrx 0 Vrz]'*(phir==0)+[0 Vr0 0]'*(phir~=0); % LEO satellite velocity

Rtx = Rt(1);
Rtz = Rt(3);
Vtx=-Vt0/sqrt(1+(Rtx/Rtz)^2);
Vtz=-Rtx/Rtz*Vtx;
Vt=[Vtx 0 Vtz]'*(phit==0)+[0 Vt0 0]'*(phit~=0); % GPS satellite velocity

Dt = zeros(length(y), length(x));     % Delay in spatial map
fd = zeros(length(y), length(x));     % Doppler in spatial map
sigma0 = zeros(length(y), length(x)); % Bistatic radar cross-section in spatial map

switch distribution
    case 'Delta'
        % Direction vector for the slope in the spatial map
        sx=zeros(length(y), length(x));
        sy=zeros(length(y), length(x));
        % Fresnel reflection coefficiente in the spatial map (Right2Left)
        rRL=zeros(length(y), length(x));
        for i=1:length(x);
            for j=1:length(y);
                r=[x(i);y(j);0]; % r point in space
                Rts=(r-Rt); % Distance between GPS satellite and r point
                Rsr=(Rr-r); % Distance between LEO satellite and r point
                m=Rts/norm(Rts); % Incidence direction unit vector
                n=Rsr/norm(Rsr); % Scattering direction unit vector
                fd(j, i)=1/LAMBDA*(m'*Vt-n'*Vr); % Signal Doppler relecting in point r
                q=K*(n-m); % Scattering vector
                sx(j, i)=-q(1)/q(3);
                sy(j, i)=-q(2)/q(3);
                Dt(j, i)=(norm(Rts)+norm(Rsr))/C_light; % Signal code delay relecting in point r
                theta_l=pi/2-acos((q'*n)/norm(q));
                rH=(e_water*sin(theta_l)-sqrt(e_water-cos(theta_l)^2))/(e_water*sin(theta_l)+sqrt(e_water-cos(theta_l)^2));
                rV=(sin(theta_l)-sqrt(e_water-cos(theta_l)^2))/(sin(theta_l)+sqrt(e_water-cos(theta_l)^2));
                rRL(j, i)=(rV-rH)/2;
            end
        end
        % Search for specular reflection point
        [a, b]=find(Dt==min(Dt(:))); 
        Ysp=a(1);
        Xsp=b(1);
        Ssp=[sx(Xsp, Ysp) sy(Xsp, Ysp)]; % Vector de dirección de pendiente en el punto especular.
        
        for i=1:length(x);
            for j=1:length(y);
                % BCRS generation
                P=sloPDFdelta([sx(j, i);sy(j, i)], dx, dy, Ssp);
                sigma0(j, i)=pi*abs(rRL(j, i))^2*(sx(j, i)^2+sy(j, i)^2+1)^2*P;
            end
        end
        
    case 'Gaussian'
        [~, MSSu, MSSc, ~]=WindElf(U10, AGE, kco); % Valores cuadrático medios de las pendientes en upwind y crosswind,  usando el modelo de Elfouhaily.
        MSSx=MSSu;
        MSSy=MSSc;
        C=[MSSx    0;
            0    MSSy]; % Covariance matrix for slope distribution
        for i=1:length(x)
            for j=1:length(y)
                r=[x(i);y(j);0];  % r point in space
                Rts=(r-Rt); % Distance between GPS satellite and r point
                Rsr=(Rr-r); % Distance between LEO satellite and r point
                Dt(j, i)=(norm(Rts)+norm(Rsr))/C_light; % Signal code delay relecting in point r
                
                m=Rts/norm(Rts); % Incidence direction unit vector
                n=Rsr/norm(Rsr); % Scattering direction unit vector
                fd(j, i)=1/LAMBDA*(m'*Vt-n'*Vr); % Signal Doppler relecting in point r
                
                q=K*(n-m); % Scattering vector
                s=-q(1:2)/q(3); % Direction vector of the slope in point r
                
                
                theta_l=pi/2-acos((q'*n)/norm(q));
                rH=(e_water*sin(theta_l)-sqrt(e_water-cos(theta_l)^2))/(e_water*sin(theta_l)+sqrt(e_water-cos(theta_l)^2));
                rV=(sin(theta_l)-sqrt(e_water-cos(theta_l)^2))/(sin(theta_l)+sqrt(e_water-cos(theta_l)^2));
                rRL=(rV-rH)/2; % Coeficiente de reflexion de Fresnel en el punto r (Right2Left).
                
                % BCRS generation
                P=sloPDFgauss(s, C, phi0);
                sigma0(j, i)=pi*abs(rRL)^2*norm(q)^4/q(3)^4*P;
            end
        end
        % Search for specular reflection point
        [a, b]=find(Dt==min(Dt(:)));
        Ysp=a(1);
        Xsp=b(1);
end


tc=1e-3/1023; % Chip time
Dtdif=(Dt-Dt(Ysp, Xsp))/tc; % Delay map, relative to the SP delay [chips]
fddif=(fd-fd(Ysp, Xsp)); % Doppler map, relative to the SP Doppler [Hz]

end

