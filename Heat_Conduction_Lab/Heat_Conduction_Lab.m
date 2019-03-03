function Heat_Conduction_Lab()
close all
clear all
clc

load('SphereData.mat')

%%%  --->  Data Analysis For all Four Trials  <---  %%%
% time_heat = [];
% theta_heat = [];
% lumped_cu_heat = [];
% T_heat = [];
% T_Center_Codel_heat = [];
% T_Surface_Model_heat = [];
% lumped_ss_heat = [];
% 
% time_cool = [];
% theta_cool = [];
% lumped_cu_cool = [];
% T_cool = [];
% T_Center_Codel_cool = [];
% T_Surface_Model_cool = [];
% lumped_ss_cool = [];


% for trial = 1:length(sphere_data)
%     
%     
% end

[t_heating,order_h] = sort([sphere_data(1).time;sphere_data(3).time]);
[t_cooling,order_c] = sort([sphere_data(2).time;sphere_data(4).time]);

T_heating = [sphere_data(1).ss_cent,sphere_data(1).ss_surf,sphere_data(1).cu_surf,sphere_data(1).cu_cent;
             sphere_data(3).ss_cent,sphere_data(3).ss_surf,sphere_data(3).cu_surf,sphere_data(3).cu_cent];
         
T_cooling = [sphere_data(2).ss_cent,sphere_data(2).ss_surf,sphere_data(2).cu_surf,sphere_data(2).cu_cent;
             sphere_data(4).ss_cent,sphere_data(4).ss_surf,sphere_data(4).cu_surf,sphere_data(4).cu_cent];
    
T_heating = T_heating(order_h,:);
T_cooling = T_cooling(order_c,:);

plot_data(t_heating,T_heating,'Heating')
plot_data(t_cooling,T_cooling,'Cooling')

end


function plot_data(t,T,trial)
% default plotting style
fs = 16;    % font size
ms = 10;    % marker size
lw = 2;     % line width

% load data from struct
% t = sphere_data(trial).time;    % Time (sec)
% T = [sphere_data(trial).ss_cent,sphere_data(trial).ss_surf,sphere_data(trial).cu_surf,sphere_data(trial).cu_cent];    % Sphere temperatures (C)


%%%  --->  Plot the temperature data  <---  %%%
figure
T_i = mean(T(1,:));%60;

fprintf("Ti = %.2f\n",T_i)

if trial == "Heating"
    T_inf = 60.0;
else
    T_inf = 29.2;
end

fprintf("T_inf = %.2f\n",T_inf)

theta = (T - T_inf)./(T_i - T_inf);
semilogy(t,theta,'linewidth',lw)
hold on

theta = log10(theta);

t_crop = t(40:80);
T1_crop = theta(40:80,1);
T2_crop = theta(40:80,2);
T3_crop = theta(40:80,3);
T4_crop = theta(40:80,4);

lin_T1 = fitlm(t_crop,T1_crop,'linear');
b_T1 = table2array(lin_T1.Coefficients(1,'Estimate'));
m_T1 = table2array(lin_T1.Coefficients(2,'Estimate'));

lin_T2 = fitlm(t_crop,T2_crop,'linear');
b_T2 = table2array(lin_T2.Coefficients(1,'Estimate'));
m_T2 = table2array(lin_T2.Coefficients(2,'Estimate'));

lin_T3 = fitlm(t_crop,T3_crop,'linear');
b_T3 = table2array(lin_T3.Coefficients(1,'Estimate'));
m_T3 = table2array(lin_T3.Coefficients(2,'Estimate'));

lin_T4 = fitlm(t_crop,T4_crop,'linear');
b_T4 = table2array(lin_T4.Coefficients(1,'Estimate'));
m_T4 = table2array(lin_T4.Coefficients(2,'Estimate'));


semilogy(t_crop,t_crop.*m_T1+b_T1,'k','LineWidth',2)%,t_crop,t_crop.*m_T2+b_T2,'k',t_crop,t_crop.*m_T3+b_T3,'k',t_crop,t_crop.*m_T4+b_T4,'k')


legend('Steel center','Steel surface','Copper center','Copper surface')
set(gca,'FontSize',fs-2);
xlabel('Time (sec)','Fontsize',fs)
ylabel('Temperature (\circC)','Fontsize',fs)
plot_title = sprintf('Sphere Temperatures : Trial %s',trial);
title(plot_title)


%%%  --->  Fit the copper sphere cooling curves  <---  %%%
figure
%     T_inf = 29.2;
%     T_i = 60.3;
r = 0.023875;   % radius (m)
V = (4/3)*pi*r^3;
As = 4*pi*r^2;
h = 2400; %initial value
rho = 8933; % kg/m^3
cp = 385;   % J/kgK
k = 401; % W/m*k

% Calculations
Bi = h*r/k;
fprintf('The Bi for copper, trial %s: %f\n',trial,Bi)
tau = rho*cp*V/(h*As);  % timne constant (sec)
T_copper = T_inf + (T_i-T_inf)*exp(-t/tau);

% Plot data
plot(t(1:10:end),T(1:10:end,3),'ro','linewidth',lw,'markersize',ms)
hold on
plot(t(1:10:end),T_copper(1:10:end),'linewidth',lw)
hold on

% Plot parameters
cu_title = sprintf('Copper sphere cooling : Trial %s',trial);
title(cu_title,'Fontsize',fs)
xlabel('Time (sec)','Fontsize',fs)
ylabel('Temperature (\circC)','Fontsize',fs)
legend('Center temperature','Lumped, h = 2400 W/m^2K')
set(gca,'FontSize',fs-2);


%h = determine_h(sphere_data(trial).cu_surf,T_i,T_inf,rho,V,As,cp,trial);

%log_dimensionless_temp(t,sphere_data(trial).cu_surf,sphere_data(trial).cu_cent,T_i,T_inf,trial,'Cu')

%lin_fit_logs(t,sphere_data(trial).cu_surf,sphere_data(trial).cu_cent,trial,'SS')

%determine_Bi(t,sphere_data(trial).cu_surf,T_i,T_inf,r,rho,cp,k,trial,'Cu')


%%%  --->  Steel sphere cooling curves  <---  %%%
figure
plot(t(1:10:end),T(1:10:end,1),'ro','linewidth',lw,'markersize',ms)    % plot the center data
hold on
plot(t(1:10:end),T(1:10:end,2),'bs','linewidth',lw,'markersize',ms)    % plot the edge data

h = 1550;   % adjustable guess for h --> 1250 initial value
rho = 8238; % kg/m^3
cp = 468;   % J/kgK
k = 13.4;   % W/mK
alpha = k./(rho*cp);
Bi = h*r/k;
r = 0.023875;   % radius (m)
V = (4/3)*pi*r^3;
As = 4*pi*r^2;
tau = rho*cp*V/(h*As);  % lumped time constant

% Bi number for steel
fprintf("Bi number for SS, trial %s: %f\n",trial,Bi)

% Find the first 15 series terms
clear lambda lambda_new
lambda = linspace(0,50,1000000);
% plot(lambda/pi,cot(lambda))
% hold on
% plot(lambda/pi,(1-Bi)./lambda)
% break

tol = 1e-4;
lambda = lambda(abs(cot(lambda) - (1-Bi)./lambda) < tol);  % find intersections
j = 1;  % index for temporary lambda vector
for i = 1:length(lambda) - 1    % remove duplicates
    if abs(lambda(i)-lambda(i+1)) > tol
        lambda_new(j) = lambda(i);
        j = j+1;
    end
end

lambda = lambda_new;
disp("Lambda:")
disp(lambda)
% length(lambda);      % check that there are 15 values
% plot(lambda,cot(lambda),'rs','linewidth',lw,'markersize',ms)
% ylim([-5 5])
% break


A = 2*(sin(lambda) - lambda.*cos(lambda))./(lambda - sin(lambda).*cos(lambda));
f = sin(0.9*lambda)./(0.9*lambda);    % position function f at the surface r = ro
Fo = alpha*t./r.^2;
T_center_model = ones(size(t)); % pre-allocate
T_surface_model = ones(size(t));

for i = 1:length(t)
    T_center_model(i) = T_inf+(T_i-T_inf)*sum(A.*exp(-lambda.^2*Fo(i)));
    T_surface_model(i) = T_inf+(T_i-T_inf)*sum(A.*exp(-lambda.^2*Fo(i)).*f);
end

% Plot series values
plot(t,T_center_model,'g','linewidth',lw)   % plot the model
hold on
plot(t,T_surface_model,'g','linewidth',lw)   % plot the model
hold on

% Plot lumped capacitance
T_steel = T_inf + (T_i-T_inf)*exp(-t/tau); % lumped
plot(t(1:10:end),T_steel(1:10:end),'b')

% figure properties
xlim([0 max(t)])
ss_title = sprintf('Steel sphere cooling : Trial %s',trial);
title(ss_title,'Fontsize',fs)
xlabel('Time (sec)','Fontsize',fs)
ylabel('Temperature (\circC)','Fontsize',fs)
legend('Center data','Surface data','Center Model','Surface Model','Lumped Capacitence')
set(gca,'FontSize',fs-2);


%log_dimensionless_temp(t,sphere_data(trial).ss_surf,sphere_data(trial).ss_cent,T_i,T_inf,trial,'SS')

%lin_fit_logs(t,sphere_data(trial).ss_surf,sphere_data(trial).ss_cent,trial,'SS')
end

function log_dimensionless_temp(t,T1,T2,T_i,T_inf,trial,material)
T1 = log((T1-T_inf) ./ (T_i-T_inf));
T2 = log((T2-T_inf) ./ (T_i-T_inf));

figure('Name','Dimensionless Temp')
plot(t(1:10:end),T1(1:10:end),'ro',t(1:10:end),T2(1:10:end),'bo')%,t,T_i,t,T_inf)
hold on

t_crop = t(40:200);
T1_crop = T1(40:200);
T2_crop = T2(40:200);

% linear fits
lin_br = fitlm(t_crop,T1_crop,'linear');
b_T1 = table2array(lin_br.Coefficients(1,'Estimate'));
m_T1 = table2array(lin_br.Coefficients(2,'Estimate'));

lin_br = fitlm(t_crop,T2_crop,'linear');
b_T2 = table2array(lin_br.Coefficients(1,'Estimate'));
m_T2 = table2array(lin_br.Coefficients(2,'Estimate'));

% plot(t_crop,T1_crop,'bo',t_crop,T2_crop,'go')
% hold on
if material == "SS"
    t_plot = [t_crop-20:t_crop+100];
else
    t_plot = [t_crop-20:t_crop+50];
end
plot(t_plot,t_plot*m_T1+b_T1,'k',t_plot,t_plot*m_T2+b_T2,'k')

plot_title = sprintf('Log of Dimensionless Temp: %s : Trial %d',material,trial);
title(plot_title)
xlabel('Time (s)','FontSize',16)
ylabel('Temperature (\circC)','FontSize',16)
legend('Surface','Center')
end


function lin_fit_logs(t,T1,T2,trial,material)
t = t(20:40);
T1 = T1(20:40);
T2 = T2(20:40);

% linear fit for brass
lin_br = fitlm(t,T1,'linear');
b_T1 = table2array(lin_br.Coefficients(1,'Estimate'));
m_T1 = table2array(lin_br.Coefficients(2,'Estimate'));

lin_br = fitlm(t,T2,'linear');
b_T2 = table2array(lin_br.Coefficients(1,'Estimate'));
m_T2 = table2array(lin_br.Coefficients(2,'Estimate'));

figure
plot_title = sprintf("Linear fir for %s, trial %d",material,trial);

plot(t,T1,'bo',t,T2,'go')
hold on
plot(t,t*m_T1+b_T1,'k',t,t*m_T2+b_T2,'k')

title(plot_title,'Fontsize',16)
xlabel('Time (s)','linewidth',2)
ylabel('ln(\theta)','linewidth',2)

end


function h = determine_h(T_sur,T_i,T_inf,rho,V,A_s,C,trial)

q = (rho*V*C.*(T_i-T_inf));

rhs = A_s.*(T_sur(1:20)-T_inf);

h = rhs\q;

fprintf('H for trial %d: %.2f\n',trial,h)

end


function determine_Bi(t,T_sur,T_i,T_inf,r,rho,C,k,trial,material)

theta = (T_sur-T_inf) ./ (T_i-T_inf);
Bi = log(theta) .* ((-(r^2)*rho*C) ./ (t*k));

Bi = mean(Bi);
fprintf('The Bi number for %s, trial %d: %.3f\n',material,trial,Bi)

end

















