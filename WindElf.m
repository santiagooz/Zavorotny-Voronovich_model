%% project: Zavorotny-Voronovich_model
%
%
% WindElf - generates the mean square slope values of the ocean surface
% according to the Elfouhaly wave model

function [ mss,mssc,mssu,mssT ] = WindElf( U10, AGE,kco )

k=linspace(1e-3,1e4,1000000);

alpha=0.0081;
beta=1.25;
g=9.82;

switch AGE
    case 1 % fully developed sea
        OMEGAc=0.84;
    case 2 % mature
        OMEGAc=1;
    case 3 % young, de 2 a 5 
        OMEGAc=2;
end

Cd=0.00144;
up=sqrt(Cd)*U10;
a0=0.1733;
ap=4;
km=370;
cm=0.23;
am=0.13*up/cm;
gamma=1.7*(OMEGAc<=1)+(1.7+6*log10(OMEGAc))*(OMEGAc>1);
sigma=0.08*(1+4*OMEGAc^(-3));
alphap=0.006*OMEGAc^0.55;
alpham=0.01*(1+log(up/cm))*(up<=cm)+0.01*(1+3*log(up/cm))*(up>cm);
k0=g/U10^2;
kp=k0*OMEGAc^2;
cp=sqrt(g/kp);
c=sqrt((g./k).*(1+(k/km).^2));

Lpm=exp(-1.25*(kp./k).^2);
GAMMA=exp(-1/(2*sigma^2)*(sqrt(k/kp)-1).^2);
Jp=gamma.^GAMMA;
Fp=Lpm.*Jp.*exp(-OMEGAc/sqrt(10)*((sqrt(k/kp)-1)));
Fm=Lpm.*Jp.*exp(-0.25*(k/km-1).^2);
Bl=0.5*alphap.*(cp./c).*Fp;
Bh=0.5*alpham.*(cm./c).*Fm;
S=(Bl+Bh)./k.^3;

dk=k(2)-k(1);
mssT=sum(k.^2.*S)*dk;

a= find(abs(k-kco)<1e-2);

mss=sum(k(1:a(1)).^2.*S(1:a(1)))*dk;
mssc=mss/1.7;
mssu=mss-mssc;
end