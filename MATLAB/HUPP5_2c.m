
close all
clear all
clc

N=256;
sidlaengd_Plan1=30e-6;
omega_in=4e-6;
lambda_noll=1550e-9;
k_noll=2*pi/lambda_noll;
a=sidlaengd_Plan1/N;

xvekt=-N/2*a:a:(N/2-1)*a;
yvekt=xvekt;
[xmat,ymat]=meshgrid(xvekt,yvekt);
rmat=sqrt(xmat.^2+ymat.^2);

% brytningsindexvariation i (x,y)-led
n_core=1.51; % i kärnan
n_clad=1.50; % i höljet
D_core=3e-6;
nmat=(rmat<=D_core/2)*n_core+(rmat>D_core/2)*n_clad;
figure(1)
imagesc(xvekt*1e6,yvekt*1e6,nmat)
xlabel('x [µm]')
ylabel('y [µm]')
colormap(jet)
colorbar
title('Brytningsindexvariation')
drawnow


% dämpmatris (behöver normalt inte ändras)
r_daemp_start=0.8*N/2*a; % ut till denna radie sker ingen dämpning
kantvaerde=0.8; % värdet på daempmat vid kanten av beräkningsfönstret (längs xy-axlarna)
daempmat=(rmat<=r_daemp_start)*1+(rmat>r_daemp_start).*(1-(1-kantvaerde)/(N/2*a-r_daemp_start)^2.*(rmat-r_daemp_start).^2);
daempmat(daempmat<0)=0;
figure(2)
mesh(xvekt*1e6,yvekt*1e6,daempmat)
xlabel('x [µm]')
ylabel('y [µm]')
title('Daempmat')
colormap(jet)
drawnow
disp('Tryck valfri tangent för att fortsätta!')
pause

% startfält "längst till vänster"
delta_y = 2.35e-6;
E_start=exp(-(xmat.^2+(ymat-delta_y).^2)/omega_in^2);

% total propagationssträcka och BPM-steg-storlek
L=1000e-6;
delta_z=2e-6;
Lvekt=delta_z:delta_z:L;

E1=E_start;
E_sida=zeros(N,length(Lvekt));
I_sida_norm=zeros(N,length(Lvekt));
steg_nummer=0;
for akt_L=Lvekt
    steg_nummer=steg_nummer+1;
    
    E2=BPM_steg(E1,delta_z,N,a,lambda_noll,nmat,daempmat); % Fältfördelning i Plan 2
    
    I2=abs(E2).^2; % Intensitetsfördelning i Plan 2
    
    % Vektorer/matriser för plottning sett "från sidan"
    E2_laengs_yaxeln=E2(:,N/2+1);
    E_sida(:,steg_nummer)=E2_laengs_yaxeln;
    I2_laengs_yaxeln=abs(E2_laengs_yaxeln).^2;
    I2_laengs_yaxeln_norm=I2_laengs_yaxeln/max(I2_laengs_yaxeln);
    I_sida_norm(:,steg_nummer)=I2_laengs_yaxeln_norm;
    
    if steg_nummer <50 | rem(steg_nummer,5)==0 % för att snabba på simuleringen plottas inte alla BPM-steg
        
        figure(10)
        imagesc(xvekt*1e6,yvekt*1e6,I2)
        hold on
        plot(D_core/2*cos(linspace(0,2*pi,50))*1e6,D_core/2*sin(linspace(0,2*pi,50))*1e6,'Color',[1 1 1]*0.6,'LineWidth',2)
        title(['Intensitet i tvärsnitt efter ' num2str(akt_L*1e3) ' mm propagation' ])
        axis('square')
        xlabel('x [µm]')
        ylabel('y [µm]')
        colormap(jet)
        hold off
        drawnow
        %pause
        
        
        figure(11)
        sc=imagesc(Lvekt*1e3,yvekt*1e6,I_sida_norm)
        title(['I_-sida: Normerad intensitet efter ' num2str(akt_L*1e3) ' mm propagation' ])
        xlabel('z [mm]')
        ylabel('y [µm]')
        colormap(jet)
        drawnow
        %pause
        
    end
    
    E1=E2;
end

dt = datatip(sc,0,0);

