function Lab2()

load('Data.mat')

%%% Initial Temps
Ti_FF = 34.85;
Ti_NI = 35.44;
Ti_NU = 36.34;

%%% Room temp
Tinf = 23.23;

%%% Plot of raw temp data
figure('Name','Raw Data')
plot(time,forced_fanned,'.',time,nat_uninsulated,'.',time,nat_insulated,'.',time,room_temp,'.')
xlabel('Time (s)')
ylabel('Temp (\circC)')
title('Raw Cooling Curves')
legend('Forced Fanned','Natural Uninsulated','Natural Insulated','Room Temp')

%%% Theta values for the three cases
theta_FF = (forced_fanned-Tinf)./(Ti_FF-Tinf);
theta_NI = (nat_insulated-Tinf)./(Ti_NI-Tinf);
theta_NU = (nat_uninsulated-Tinf)./(Ti_NU-Tinf);


%%% Forced fan Calcs
h_T_FF = LCM_Fit(Ti_FF,theta_FF,time);
fprintf('h for forced fanned: %.2f\n',h_T_FF)
figure('Name','Forced Fan LCM'); plot(time,theta_FF,'.',time,LCM_Temps(Ti_FF,time,h_T_FF));

%%% Natural Insulated calcs
h_T_NI = LCM_Fit(Ti_NI,theta_NI,time);
fprintf('h for Natural Insulated: %.2f\n',h_T_NI)
figure('Name','Natural Insulated LCM'); plot(time,theta_NI,'.',time,LCM_Temps(Ti_NI,time,h_T_NI));

%%% Natural Un-Insulated calcs
h_T_NU = LCM_Fit(Ti_NU,theta_NU,time);
fprintf('h for Natural Un-Insulated: %.2f\n',h_T_NU)
figure('Name','Natural Un-Insulated LCM'); plot(time,theta_NU,'.',time,LCM_Temps(Ti_NU,time,h_T_NU));



end



function theta = LCM_Temps(Ti,t,h)
Tinf = 23.23; % deg C
rho = 993.05; % kg/m^3
cp = 4178; % J/kg*K
V = 3.656 * 10^-4; % m^3
A = .0198; % m^2

theta = exp((-h.*A.*t)./(rho.*cp.*V));
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


