function [mach_cd,Cddelta] = calc_cdd2(data)
  % funcion que calcula Cddelta segun Resu de RBD
  % Cd = Cd0 + Cddelta * delta^2
  % data el el archivo de la corrida a leer el Cd para calcular el Cdd2
  % ej:data = load('/home/zeeburg/Documents/CIMEC/Tesis/estimacion/Resu_RBD/Force_coef_proc_B06_alpha80.txt');
	 
  mach_cd = data(:,2);
  alfa = data(:,3);
  alfa_deg = rad2deg(alfa);
  delta2 = data(:,5);
  Cd = data(:,6);
  
    % Interp del Cd0_exp(M) en func del vec mach
  % Cd0_exp_interp =  interp1(M, Cd0_exp, mach);

  % Cd0_exp_interp sale de corrida previa a delta 0, B05
  data0 = load('/home/zeeburg/Documents/CIMEC/Tesis/estimacion/Resu_RBD/Force_coef_proc_B05.txt');
  mach0 = data0(:,2);
  Cd0 = data0(:,6);
  % ver cual es el max de los min de los mach: mach o mach_cd
  mach = [max(mach0):-0.05:min(mach0)];
  %mach_cd = [max(mach_cd):-0.05:min(mach_cd)];
  
  %Cd_interp =  interp1(mach_cd, Cd, mach);
  Cd0_interp =  interp1(mach0, Cd0, mach_cd);

  x = Cd - Cd0_interp;
  Cddelta = zeros(1,length(x));
  
  % Calculo de Cddelta solo si delta > 0.5 grad
  for i = 1:length(x)
    if(x(i) > 0.05*pi/180)
      Cddelta(i) = abs((Cd(i) - Cd0_interp(i))) / delta2(i);
    end
  end
  
end
