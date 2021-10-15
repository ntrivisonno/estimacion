clear all;close all;
% ex Caso Sintético - Estimacion Cd0 y Cdd2 adaptado a Resu_RBD
% todo armado para que guarde con nombres automaticamente
tStart = tic; %timer de multivariables
% Datos Cd McCoy
% Cd0
M = [0, 0.8, 0.85, 0.9, 0.95, 1.0, 1.05, 1.1, 1.2, 1.4, 1.6, 1.8, 2.0, 2.2, 2.5]';
Cd0_exp = [0.14, 0.14, 0.142, 0.16, 0.24, 0.43, 0.449, 0.447, 0.434, 0.41, 0.385, 0.365, 0.35, 0.339, 0.32]';

% Cdd2
Md2 = [0, 0.95, 1, 1.05, 1.1, 1.2, 1.4, 1.6, 1.8, 2.0, 2.2, 2.5]';
Cdd2_exp = [2.9, 2.9, 3.0, 3.1, 3.6, 6.5, 7.6, 7.3, 6.8, 6.1, 5.4, 4.4]';

% Cd = Cd0 + Cdd2 * delta^2

% Colocar nombre del caso a analziar el CD
caso = 'B06_alpha15';
% maquinaria para cargar datos en funcion del nombre del caso
folder_path=('/home/zeeburg/Documents/CIMEC/Tesis/estimacion/Resu_RBD/');
path_file = fullfile([folder_path],['Force_coef_proc_',caso,'.txt']);
fprintf('\nSe cargan los datos de CD del archivo:\n%s\n',path_file)

# Time,   Mach,     alfa,     beta,     delta2,     Cd,     CL_alfa,     Cn_p_alfa,     Cn_q_alfa 
data = load(path_file);
mach_cd = data(:,2);
alfa_cd = data(:,3);
delta2 = data(:,5);
Cd = data(:,6);

% Cd0 sale de corrida previa a delta 0, B05
data0 = load('/home/zeeburg/Documents/CIMEC/Tesis/estimacion/Resu_RBD/Force_coef_proc_B05.txt');
mach0 = data0(:,2);
Cd0 = data0(:,6);
% ver cual es el max de los min de los mach: mach o mach_cd
mach = [max(mach0):-0.05:min(mach0)];
%mach_cd = [max(mach_cd):-0.05:min(mach_cd)];

Cd0_interp =  interp1(mach0, Cd0, mach_cd);

x = Cd - Cd0_interp;
Cddelta = zeros(1,length(x));

% Calculo de Cddelta solo si delta > 0.5 grad
for i = 1:length(x)
  if(delta2(i) > 0.05*pi/180)
    Cddelta(i) = (Cd(i) - Cd0_interp(i)) / delta2(i);
  end
end

% Figures
figure();plot(mach_cd,Cd,'o-');
title_Cd = sprintf('Cd %s',caso);
title(title_Cd);xlabel('Mach','FontSize',14);ylabel('Cd','FontSize',14)
figname_Cd_interp = sprintf('Cd_%s',caso);xlim([0 2.5])
print(fullfile('Figures',[figname_Cd_interp,'.png']),'-dpng')
% gnumeric Cd
gnumeric_Cd_interp = [mach_cd Cd]; % se crea vector para exportar como csv
cvs_Cdinterp = sprintf('./csv/cvsgnum_%s',figname_Cd_interp);
csvwrite (cvs_Cdinterp,gnumeric_Cd_interp);


figure();plot(M,Cd0_exp,mach_cd,Cd0_interp,'o');
title_Cd0 = sprintf('Cd0 %s',caso);
title(title_Cd0);xlabel('Mach','FontSize',14);ylabel('Cd0','FontSize',14);legend('Real','Estimacion');
figname_Cd0_interp = sprintf('Aprox_Cd0_%s',caso);xlim([0 2.5])
print(fullfile('Figures',[figname_Cd0_interp,'.png']),'-dpng')
% gnumeric Cdd2
gnumeric_Cd0_interp = [mach_cd Cd0_interp]; % se crea vector para exportar como csv
cvs_Cd0interp = sprintf('./csv/cvsgnum_%s',figname_Cd0_interp);
csvwrite (cvs_Cd0interp,gnumeric_Cd0_interp);


figure();plot(Md2,Cdd2_exp,mach_cd,Cddelta,'+');
title_Cdd2=sprintf('Cdd2 %s',caso);
title(title_Cdd2);legend('Real','Estimacion');xlabel('Mach','FontSize',14);ylabel('Cdd2','FontSize',14);
figname_Cdd2 = sprintf('Aprox_Cdd2_%s',caso);xlim([0 2.5])
print(fullfile('Figures',[figname_Cdd2,'.png']),'-dpng')
% gnumeric Cdd2
gnumeric_Cdd2_interp = [mach_cd Cddelta']; % se crea vector para exportar como csv
cvs_Cdd2interp = sprintf('./csv/cvsgnum_%s',figname_Cdd2);
csvwrite (cvs_Cdd2interp,gnumeric_Cdd2_interp);


% gnumeric Cd0 McCoy
gnumeric_Cd0_McCoy = [M Cd0_exp]; % se crea vector para exportar como csv
%cvs_Cd0_McCoy = sprintf('cvsgnum_%s',figname_Cd0_interp);
csvwrite ('./csv/csvgnum_Cd0_McCoy',gnumeric_Cd0_McCoy);
% gnumeric Cdd2 McCoy
gnumeric_Cdd2_McCoy = [Md2 Cdd2_exp]; % se crea vector para exportar como csv
%cvs_Cdd2_McCoy = sprintf('cvsgnum_%s',figname_Cdd2);
csvwrite ('./csv/csvgnum_Cdd2_McCoy',gnumeric_Cdd2_McCoy);

time = toc(tStart);
fprintf('*-----------------------------------------------*\n')
fprintf('\n\nFIN! - OK - time = %d[s].\n',time)
