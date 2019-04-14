function Lab2()
clc
load('Data.mat')

%%% Initial Temps
Ti_FF = 34.85;
Ti_NI = 35.44;
Ti_NU = 36.34;

%%% Room temp
Tinf = 23.23;

%%% Parameters of air @300K
Pr_air = .707;
v_air = 15.89 *10^-6; % m^2/s
cp_air = 1007; %J/Kg*K
k_air = 26.23 *10^-3; %W/m*K

%%% Cup of water geometry
D = .074; %m
r = .037; %m
L = .085; %m
thickness = .0003; %m


%%% Theta values for the three cases
theta_FF = (forced_fanned-Tinf)./(Ti_FF-Tinf);
theta_NI = (nat_insulated-Tinf)./(Ti_NI-Tinf);
theta_NU = (nat_uninsulated-Tinf)./(Ti_NU-Tinf);
ln_theta_FF = log((forced_fanned-Tinf)./(Ti_FF-Tinf));
ln_theta_NI = log((nat_insulated-Tinf)./(Ti_NI-Tinf));
ln_theta_NU = log((nat_uninsulated-Tinf)./(Ti_NU-Tinf));


%%% Plot of raw temp data
figure('Name','Raw Data')
plot(time,forced_fanned,'.',time,nat_uninsulated,'.',time,nat_insulated,'.',time,room_temp,'.')
xlabel('Time (s)')
ylabel('Temp (\circC)')
title('Raw Cooling Curves')
legend('Forced Fanned','Natural Uninsulated','Natural Insulated','Room Temp')


%%% Forced fan Calcs
Print_Header('Forced Fan')
h_LCM_FF = LCM_Fit(Ti_FF,theta_FF,time);
fprintf('h_LCM = %.2f\n',h_LCM_FF)
figure('Name','Forced Fan LCM');
plot(time,theta_FF,'.',time,LCM_Temps(Ti_FF,time,h_LCM_FF));
xlabel('Time (s)')
ylabel('\theta')
title('\theta vs time for Forced Convection')


%%% Natural Insulated calcs
Print_Header('Natural Insulated')
h_LCM_NI = LCM_Fit(Ti_NI,theta_NI,time);
fprintf('h_LCM = %.2f\n',h_LCM_NI)
figure('Name','Natural Insulated LCM');
plot(time,theta_NI,'.',time,LCM_Temps(Ti_NI,time,h_LCM_NI));
xlabel('Time (s)')
ylabel('\theta')
title('\theta vs time for Insulated Natural Convection')

[beta_NI,Gr_NI,Ra_NI,Nu_condit_NI,Nu_NI,h_NI] = Natural_Convection_Values_Cylinder(Ti_NI,Tinf,Pr_air,v_air,k_air,D);
fprintf('Beta = %.3f\nGr = %.3e\nRa = %.3e\nD>Nu_condit = %d\n',beta_NI,Gr_NI,Ra_NI,D>=Nu_condit_NI)
fprintf('Nu = %.3f\nh_nat = %.3f\n',Nu_NI,h_NI)


%%% Natural Un-Insulated calcs
Print_Header('Natural Un-Insulated')
h_LCM_NU = LCM_Fit(Ti_NU,theta_NU,time);
fprintf('h_LCM = %.2f\n',h_LCM_NU)
figure('Name','Natural Un-Insulated LCM');
plot(time,theta_NU,'.',time,LCM_Temps(Ti_NU,time,h_LCM_NU));
xlabel('Time (s)')
ylabel('\theta')
title('\theta vs time for Uninsulated Natural Convection')

[beta_NU,Gr_NU,Ra_NU,Nu_condit_NU,Nu_NU,h_NU] = Natural_Convection_Values_Cylinder(Ti_NU,Tinf,Pr_air,v_air,k_air,D);
fprintf('Beta = %.3f\nGr = %.3e\nRa = %.3e\nD>Nu_condit = %d\n',beta_NU,Gr_NU,Ra_NU,D>=Nu_condit_NU)
fprintf('Nu = %.3f\nh_nat = %.3f\n',Nu_NU,h_NU)

[beta_NU_top,Gr_NU_top,Ra_NU_top,Nu_NU_top,h_NU_top] = ...
    Natural_Convection_Values_TopSurface(Ti_NU,Tinf,Pr_air,v_air,k_air,D);
fprintf('\nBeta_top = %.3f\nGr_top = %.3e\nRa_top = %.3e\n',beta_NU_top,Gr_NU_top,Ra_NU_top)
fprintf('Nu_top = %.3f\nh_nat_top = %.3f\n',Nu_NU_top,h_NU_top)

h_avg = mean([h_NU,h_NU_top]);
fprintf('\nAverage h = %.3f\n',h_avg)
end


function Print_Header(label)
fprintf('\n***********************************************************************************\n')
fprintf('%s\n',label)
fprintf('***********************************************************************************\n')
end

function theta = LCM_Temps(Ti,t,h)
Tinf = 23.23; % deg C
rho_water = 993.05; % kg/m^3
cp_water = 4178; % J/kg*K
V = 3.656 * 10^-4; % m^3
A = .0198; % m^2

theta = exp((-h.*A.*t)./(rho_water.*cp_water.*V));
end

function h_out = LCM_Fit(Ti,theta,time)
mse_T = [];
h = [];
for i = 1:.05:100
    T = LCM_Temps(Ti,time,i);
    mse = immse(theta,T);
    mse_T(end+1) = mse;
    h(end+1) = i;
end
h_out = h(find(mse_T == min(mse_T)));
end

function [beta,Gr,Ra,Nu_condition,Nu,h] = Natural_Convection_Values_Cylinder(Ti,Tinf,Pr,v,k_air,D)
g = 9.8; %m/s
L = .085; %m


Tfilm = (Ti+Tinf)/2;

beta = 1/Tfilm;

Gr = (g*beta*(Ti-Tinf)*L^3) / (v^2);

Ra = Gr * Pr;

Nu_condition = (35 * L) / (Gr^.25);

Nu = .59 * Ra^.25;

h = (k_air/D) * Nu;
end

function [beta,Gr,Ra,Nu,h] = Natural_Convection_Values_TopSurface(Ti,Tinf,Pr,v,k_air,D)
g = 9.8; %m/s
As = pi/4 * D^2;
P = pi * D;
L = As/P; %m

Tfilm = (Ti+Tinf)/2;

beta = 1/Tfilm;

Gr = (g*beta*(Ti-Tinf)*L^3) / (v^2);

Ra = Gr * Pr;

Nu = .54 * Ra^.25;

h = (k_air/D) * Nu;
end



