set(0,'Defaultaxesfontsize',20)
set(0,'DefaultFigureWindowStyle','docked')
set(0,'DefaultLineLineWidth',2);
set(0,'Defaultaxeslinewidth',2)

set(0,'DefaultFigureWindowStyle','docked')

c_c = 299792358; % m/s speed of light
c_eps_0 = 8.8542149e-12; % F/m vacuum permittivity
c_eps_0_cm = c_eps_0/100; %F/cm
c_mu_0=1/c_eps_0/c_c^2;
c_q=1.60217653e-19;
c_hb=1.05457266913e-34; %Dirac constant
c_h=c_hb*2*pi;

InputParasL.E0 = 1e8;
InputParasL.we = 0;
InputParasL.t0 = 2e-12;
InputParasL.wg = 5e-13; %Guess weidth
InputParasL.phi = 0;
InputParasR = 0;

n_g = 3.5; %Time
vg = c_c/n_g*1e2;   %TWM cm/s group velocity
Lambda = 1550e-9;

plotN = 10; %define number of plot

L = 1000e-6*1e2; %cm
XL = [0,L]; %X axis value
YL = [-InputParasL.E0,InputParasL.E0]; %Y axis value

g_fwhm = 3.53e+012/10;
LGamma = g_fwhm*2*pi;
Lw0 = 0.0;
LGain = 0.1;

Nz = 500;
dz = L/(Nz-1);
dt = dz/vg;
fsync = dt*vg/dz;

Nt = floor(4*Nz); %Control total time escape
tmax = Nt*dt;
t_L=dt*Nz;  %time to travel length

z = linspace(0,L,Nz); %500 numbers from 0 to L L is 1000e-6*1e2
time = nan(1,Nt);
InputL = nan(1,Nt);
InputR = nan(1,Nt);
OutputL = nan(1,Nt);
OutputR = nan(1,Nt);

Ef = zeros(size(z)); %init Ef
Er = zeros(size(z)); %init Er

Pf = zeros(size(z));
Pr = zeros(size(z));

Efp = Ef;
Erp = Er;
Pfp = Pf;
Prp = Pr;

RL = 0; %add mirror Ef(0)= f(t)+RrEr(0)
RR=0; %add mirror Er(L) = Rl*Ef(L)

Ef1 = @SourceFct; %Calculate the Ef1 by function F
ErN = @SourceFct; %Calculate the ErN by function F

t=0;
time(1) = t;

InputL(1) = Ef1(t,InputParasL);
InputR(1) = ErN(t,InputParasR);

%OutputR(1) = Ef(Nz);
%OutputL(1) = Er(1);
OutputR(1)=Ef(Nz)*(1-RR); %Output of Of
OutputL(1)=Er(1)*(1-RL); %Output of Or

%Ef(1)= InputL(1);
%Er(Nz) = InputR(1);
Ef(1) = InputL(1)+RL*Er(1);
Er(Nz) = InputR(1)+RR*Ef(Nz);

beta_r = 0;
beta_i = 0;

beta = ones(size(z))*(beta_r+1i*beta_i);
exp_det = exp(-1i*dz*beta);

kappaStart = 1/3;
kappaStop =2/3;
kappa0=100;

figure('name','Fields')
subplot(3,1,1)
plot(z*10000,real(Ef),'r');
hold off
xlabel('z(\mum)')
ylabel('E_f')
subplot(3,1,2)
plot(z*10000,real(Er),'b');
xlabel('z(\mum)')
ylabel('E_r')
hold off
subplot(3,2,3)
plot(time*1e12,real(InputL),'r'); hold on
plot(time*1e12,real(OutputR),'r--'); hold on
plot(time*1e12,real(InputR),'b'); hold on
plot(time*1e12,real(OutputL),'b--');
xlabel('time(ps)')
ylabel('E')

hold off

for i = 2:Nt
    t = dt*(i-1);
    time(i) = t;

    InputL(i) = Ef1(t,InputParasL);
    InputR(i) = ErN(t,0);

    kappa = kappa0*ones(size(z));
    kappa(z<L*kappaStart) = 0;
    kappa(z>L*kappaStop) = 0;

    %Ef(1) = InputL(i);
    %Er(Nz) = InputR(i); Code without mirror
    Ef(1) = InputL(i)+RL*Er(1);
    Er(Nz) = InputR(i)+RR*Ef(Nz);
    Eftemp = Ef(1:Nz-1);
    % Ef(2:Nz) = fsync*Ef(1:Nz-1);
    % Er(1:Nz-1) = fsync*Er(2:Nz); Code with gain
    Ef(2:Nz) = fsync*exp_det(1:Nz-1).*Ef(1:Nz-1);%+1i*dz*kappa(2:Nz).*Er(2:Nz);
    Er(1:Nz-1) = fsync*exp_det(2:Nz).*Er(2:Nz);%+1i*dz*kappa(1:Nz-1).*Eftemp(1:Nz-1);

    %OutputR(i)=Ef(Nz);
    %OutputL(i)=Er(1); Code without mirror
    OutputR(i) = Ef(Nz)*(1-RR);
    OutputL(i) = Er(1)*(1-RL);

    Pf(1) = 0;
    Pf(Nz) = 0;
    Pr(1) = 0;
    Pr(Nz) = 0;
    Cw0 = -LGamma + 1i*Lw0;

    Tf = LGamma*Ef(1:Nz-2) + Cw0*Pfp(2:Nz-1) + LGamma*Efp(1:Nz-2);
    Pf(2:Nz-1) = (Pfp(2:Nz-1) + 0.5*dt*Tf)./(1-0.5*dt*Cw0);
    Tr = LGamma*Er(3:Nz) + Cw0*Prp(2:Nz-1) + LGamma*Erp(3:Nz);
    Pr(2:Nz-1) = (Prp(2:Nz-1)+0.5*dt*Tr)./(1-0.5*dt*Cw0);

    Ef(2:Nz-1) = Ef(2:Nz-1)-LGain*(Ef(2:Nz-1)-Pf(2:Nz-1));
    Er(2:Nz-1) = Er(2:Nz-1)-LGain*(Er(2:Nz-1)-Pr(2:Nz-1));

    if mod(i,plotN) == 0
        subplot(3,2,1)
        plot(z*10000,real(Ef),'r'); hold on
        plot(z*10000,imag(Ef),'r--'); hold off
        xlim(XL*1e4)
        ylim(YL)
        xlabel('z(\mum')
        ylabel('E_f')
        legend('\Re','\Im')
        hold off

        subplot(3,2,2)
        plot(z*10000,real(Er),'b');hold on
        plot(z*10000,imag(Er),'b--');hold off
        xlim(XL*1e4)
        ylim(YL)
        xlabel('z(\mum')
        ylabel('E_r')
        legend('\Re','\Im')
        hold off

        subplot(3,2,3);
        plot(time*1e12,real(InputL),'r');hold on
        plot(time*1e12,real(OutputR),'g');
        plot(time*1e12,real(InputR),'b');
        plot(time*1e12,real(OutputL),'m');
        xlim([0,Nt*dt*1e12])
        ylim(YL)
        xlabel('time(ps)')
        ylabel('0')
        legend('Left Input','Right Output','Right Input','Left Output','Location','east')
        hold off
        pause(0.01)

        subplot(3,2,4);
        plot(time*1e12,real(OutputR));hold on
        xlabel('time(ps)')
        ylabel('Right Output')
        hold off
           
    end
    Efp=Ef;
    Erp=Er;
    Pfp=Pf;
    Prp=Pr;
end

fftOutput = fftshift(fft(OutputR));
fftInput = fftshift(fft(InputL));
omega = fftshift(wspace(time));

subplot(3,2,5);
plot(omega,abs(fftOutput));hold on
plot(omega,abs(fftInput));hold on
xlim([-1,1]*1e13)
xlabel('THz')
ylabel('|E|')
legend('Output','Input')
hold off
    
subplot(3,2,6);
plot(omega,unwrap(angle(fftOutput)));hold on
plot(omega,unwrap(angle(fftInput)));hold on
xlabel('GHz')
ylabel('phase(E)')
legend('Output','Input')
hold off