clear all;close all;
% Caso Sintético - Cd manufacturado
% Datos
% Cd0
M = [0, 0.8, 0.85, 0.9, 0.95, 1.0, 1.05, 1.1, 1.2, 1.4, 1.6, 1.8, 2.0, 2.2, 2.5];
Cd0_exp = [0.14, 0.14, 0.142, 0.16, 0.24, 0.43, 0.449, 0.447, 0.434, 0.41, 0.385, 0.365, 0.35, 0.339, 0.32];

% Cdd2
Md2 = [0, 0.95, 1, 1.05, 1.1, 1.2, 1.4, 1.6, 1.8, 2.0, 2.2, 2.5];
Cdd2_exp = [2.9, 2.9, 3.0, 3.1, 3.6, 6.5, 7.6, 7.3, 6.8, 6.1, 5.4, 4.4];

% Cd = Cd0 + Cdd2 * delta^2

% Rango de Mach
mach = [2.5:-0.05:0]; % xq lo hace decreciente?
% Tiempo
t = [1/length(mach):1/length(mach):1];
% Ang de ataque max 10 grad
amax = 10*pi/180;% amax = deg2rad(10) convert degree2radiands
% Variación ang ataque vs tiempo +-10 grad
delta = amax*sin(4*t*pi); % delta rad
delta_deg = rad2deg(delta); % delta degree
sin_deg = rad2deg(delta/amax);
%figure
%plot(t,delta)

% Calculo Cd + una perturbación en los coef
Cd = (interp1(M, Cd0_exp, mach)+rand(1)*5e-4) + (interp1(Md2, Cdd2_exp, mach)+rand(1)*5e-4).*delta.^2;

% Interp del Cd0_exp(M) en func del vector mach, en realdiad tendria que ser Cd0 de otra corrida
Cd0_exp_interp =  interp1(M, Cd0_exp, mach);

x = Cd - Cd0_exp_interp;
Cddelta = zeros(1,length(x));

% Calculo de Cddelta solo si delta > 0.5 grad
for i = 1:length(x)
  if(x(i) > 0.5*pi/180)
    Cddelta(i) = (Cd(i) - Cd0_exp_interp(i)) / delta(i)^2;
  endif
endfor

% Figures
figure();subplot(2,1,1);plot(t,delta,'b',t,delta/amax,'r');title('Perturbacion radianes');legend('delta','sin');xlabel('t');subplot(2,1,2);plot(t,sin_deg,'b',t,delta_deg,'r');title('Perturbacion degree');legend('delta','sin');xlabel('t');print -dpng Figures/Perturbaciones_angulares.png
figure();plot(mach,Cd,'o-');title('Cd sintetico');print -dpng Figures/Cd_sintetico.png
figure();plot(mach,Cd0_exp_interp,'o-');title('Cd0 interpolado');print -dpng Figures/Cd0_interpolado.png
figure();plot(Md2,Cdd2_exp,mach,Cddelta,'+');title('Aprox Cddelta');legend('Real','Estimacion');xlabel('Mach');ylabel('Cddd2');print -dpng Figures/Aprox_Cdd2.png
