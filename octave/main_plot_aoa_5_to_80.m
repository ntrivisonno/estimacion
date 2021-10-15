% simil test_sint_adaptado_CFD.m plotea Cdd2 from aoa = 5 to 80 degrees
% guarda a todas las variables cargadas para poder plotearlas juntas en una misma grafica
clear all;close all;clc
tStart = tic; % timer de multivariables
% Datos
% Cd0
M = [0, 0.8, 0.85, 0.9, 0.95, 1.0, 1.05, 1.1, 1.2, 1.4, 1.6, 1.8, 2.0, 2.2, 2.5];
Cd0_exp = [0.14, 0.14, 0.142, 0.16, 0.24, 0.43, 0.449, 0.447, 0.434, 0.41, 0.385, 0.365, 0.35, 0.339, 0.32];
% Cdd2
Md2 = [0, 0.95, 1, 1.05, 1.1, 1.2, 1.4, 1.6, 1.8, 2.0, 2.2, 2.5];
Cdd2_exp = [2.9, 2.9, 3.0, 3.1, 3.6, 6.5, 7.6, 7.3, 6.8, 6.1, 5.4, 4.4];

data = load('/home/zeeburg/Documents/CIMEC/Tesis/estimacion/Resu_RBD/Force_coef_proc_B06_alpha5.txt');
[mach_cd_alpha5,Cddelta_alpha5] = calc_cdd2(data);
gnumeric_Cdd2_alpha5 = [mach_cd_alpha5,Cddelta_alpha5']; 
csvwrite ('./cvs_5_to_80/cvs_Cdd2interp_alpha5',gnumeric_Cdd2_alpha5);

data = load('/home/zeeburg/Documents/CIMEC/Tesis/estimacion/Resu_RBD/Force_coef_proc_B06_alpha15.txt');
[mach_cd_alpha15,Cddelta_alpha15] = calc_cdd2(data);
gnumeric_Cdd2_alpha15 = [mach_cd_alpha15,Cddelta_alpha15']; 
csvwrite ('./cvs_5_to_80/cvs_Cdd2interp_alpha15',gnumeric_Cdd2_alpha15);

data = load('/home/zeeburg/Documents/CIMEC/Tesis/estimacion/Resu_RBD/Force_coef_proc_B06_alpha30.txt');
[mach_cd_alpha30,Cddelta_alpha30] = calc_cdd2(data);
gnumeric_Cdd2_alpha30 = [mach_cd_alpha30,Cddelta_alpha30']; 
csvwrite ('./cvs_5_to_80/cvs_Cdd2interp_alpha30',gnumeric_Cdd2_alpha30);

data = load('/home/zeeburg/Documents/CIMEC/Tesis/estimacion/Resu_RBD/Force_coef_proc_B06_alpha40.txt');
[mach_cd_alpha40,Cddelta_alpha40] = calc_cdd2(data);
gnumeric_Cdd2_alpha40 = [mach_cd_alpha40,Cddelta_alpha40']; 
csvwrite ('./cvs_5_to_80/cvs_Cdd2interp_alpha40',gnumeric_Cdd2_alpha40);

data = load('/home/zeeburg/Documents/CIMEC/Tesis/estimacion/Resu_RBD/Force_coef_proc_B06_alpha50.txt');
[mach_cd_alpha50,Cddelta_alpha50] = calc_cdd2(data);
gnumeric_Cdd2_alpha50 = [mach_cd_alpha50,Cddelta_alpha50']; 
csvwrite ('./cvs_5_to_80/cvs_Cdd2interp_alpha50',gnumeric_Cdd2_alpha50);

data = load('/home/zeeburg/Documents/CIMEC/Tesis/estimacion/Resu_RBD/Force_coef_proc_B06_alpha60.txt');
[mach_cd_alpha60,Cddelta_alpha60] = calc_cdd2(data);
gnumeric_Cdd2_alpha60 = [mach_cd_alpha60,Cddelta_alpha60']; 
csvwrite ('./cvs_5_to_80/cvs_Cdd2interp_alpha60',gnumeric_Cdd2_alpha60);

data = load('/home/zeeburg/Documents/CIMEC/Tesis/estimacion/Resu_RBD/Force_coef_proc_B06_alpha70.txt');
[mach_cd_alpha70,Cddelta_alpha70] = calc_cdd2(data);
gnumeric_Cdd2_alpha70 = [mach_cd_alpha70,Cddelta_alpha70']; 
csvwrite ('./cvs_5_to_80/cvs_Cdd2interp_alpha70',gnumeric_Cdd2_alpha70);

data = load('/home/zeeburg/Documents/CIMEC/Tesis/estimacion/Resu_RBD/Force_coef_proc_B06_alpha80.txt');
[mach_cd_alpha80,Cddelta_alpha80] = calc_cdd2(data);
gnumeric_Cdd2_alpha80 = [mach_cd_alpha80,Cddelta_alpha80']; 
csvwrite ('./cvs_5_to_80/cvs_Cdd2interp_alpha80',gnumeric_Cdd2_alpha80);


%figure();plot(Md2,Cdd2_exp,mach_cd_alpha80,Cddelta_alpha80,'+');title('Aprox Cddelta alfa=80');legend('Real','Estimacion');xlabel('Mach');ylabel('Cddd2');%print -dpng Figures/Aprox_Cdd2_alfa80.png
figure();plot(Md2,Cdd2_exp,'b');hold on;
plot(mach_cd_alpha5,Cddelta_alpha5,'or');
plot(mach_cd_alpha15,Cddelta_alpha15,'+');
plot(mach_cd_alpha30,Cddelta_alpha30,'*');
plot(mach_cd_alpha40,Cddelta_alpha40,'x');
plot(mach_cd_alpha50,Cddelta_alpha50,'s');
plot(mach_cd_alpha60,Cddelta_alpha60,'d');
plot(mach_cd_alpha70,Cddelta_alpha70,'^');
plot(mach_cd_alpha80,Cddelta_alpha80,'d')
title('Cdd2 - Cuadratic drag');xlabel('Mach');ylabel('Cdd2')
legend('Real','alpha=5','alpha=15','alpha=30','alpha=40','alpha=50','alpha=60','alpha=70','alpha=80','location','SouthEast');
print -dpng Aprox_Cdd2_alpha5to80.png
%legend('Real','alpha=15','alpha=30','alpha=40','alpha=50','alpha=60','alpha=70','alpha=80','location','SouthEast')%Outside')
%xlim([0 1.3]);legend('location','NorthWest')

time = toc(tStart);
fprintf('*-----------------------------------------------*\n')
fprintf('\n\nFIN! - OK - time = %d[s].\n',time)
