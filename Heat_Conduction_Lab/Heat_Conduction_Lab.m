function Heat_Conduction_Lab()
close all
clear all
clc

% default plotting style
fs = 16;    % font size
ms = 10;    % marker size
lw = 2;     % line width

load('SphereData.mat')

%%%  --->  Data Analysis For all Four Trials  <---  %%%
for trial = 1:length(sphere_data)
    % load data from struct
    t = sphere_data(trial).time;    % Time (sec)
    T = [sphere_data(trial).ss_cent,sphere_data(trial).ss_surf,sphere_data(trial).cu_surf,sphere_data(trial).cu_cent];    % Sphere temperatures (C)
    
    
    %%%  --->  Plot the temperature data  <---  %%%
    figure
    T_i = mean(T(1,:));%60;
    
    fprintf("Ti = %.2f\n",T_i)
    
    if trial == 1 || trial ==3
        T_inf = 60.0;
    else
        T_inf = 29.2;
    end
    
    fprintf("T_inf = %.2f\n",T_inf)
    
    theta = (T - T_inf)./(T_i - T_inf);
    semilogy(t,theta,'linewidth',lw)
    legend('Steel center','Steel surface','Copper center','Copper surface')
    set(gca,'FontSize',fs-2);
    xlabel('Time (sec)','Fontsize',fs)
    ylabel('Temperature (\circC)','Fontsize',fs)
    plot_title = sprintf('Sphere Temperatures : Trial %d',trial);
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
    fprintf('The Bi for copper, trial %d: %f\n',trial,Bi)
    tau = rho*cp*V/(h*As);  % timne constant (sec)
    T_copper = T_inf + (T_i-T_inf)*exp(-t/tau);
    
    % Plot data
    plot(t(1:10:end),T(1:10:end,3),'ro','linewidth',lw,'markersize',ms)
    hold on
    plot(t(1:10:end),T_copper(1:10:end),'linewidth',lw)
    
    % Plot parameters
    cu_title = sprintf('Copper sphere cooling : Trial %d',trial);
    title(cu_title,'Fontsize',fs)
    xlabel('Time (sec)','Fontsize',fs)
    ylabel('Temperature (\circC)','Fontsize',fs)
    legend('Center temperature','Lumped, h = 2400 W/m^2K')
    set(gca,'FontSize',fs-2);
    
    
    %h = determine_h(sphere_data(trial).cu_surf,T_i,T_inf,rho,V,As,cp,trial);
    
    log_dimensionless_temp(t,sphere_data(trial).cu_surf,T_i,T_inf,trial,'Cu_Surf')
    log_dimensionless_temp(t,sphere_data(trial).cu_cent,T_i,T_inf,trial,'Cu_Surf')
    
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
    fprintf("Bi number for SS, trial %d: %f\n",trial,Bi)
    
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
    plot(t,T_surface_model,'g','linewidth',lw)   % plot the model
    
    % Plot lumped capacitance
    T_steel = T_inf + (T_i-T_inf)*exp(-t/tau); % lumped
    plot(t(1:10:end),T_steel(1:10:end))
    
    % figure properties
    xlim([0 max(t)])
    ss_title = sprintf('Steel sphere cooling : Trial %d',trial);
    title(ss_title,'Fontsize',fs)
    xlabel('Time (sec)','Fontsize',fs)
    ylabel('Temperature (\circC)','Fontsize',fs)
    legend('Center data','Surface data','Model Fits','Lumped')
    set(gca,'FontSize',fs-2);
    
    
    log_dimensionless_temp(t,sphere_data(trial).ss_surf,T_i,T_inf,trial,'SS Surf')
    log_dimensionless_temp(t,sphere_data(trial).ss_cent,T_i,T_inf,trial,'SS Cent')
    
end

end

function log_dimensionless_temp(t,T,T_i,T_inf,trial,material)
theta = (T-T_inf) ./ (T_i-T_inf);

figure('Name','Dimensionless Temp')
plot(t,log(theta),'o')%,t,T_i,t,T_inf)

plot_title = sprintf('Dimensionless Temp Plot : %s : Trial %d',material,trial);
title(plot_title)
xlabel('Time (s)')
ylabel('Temperature (\circC)')
%legend('Dimensionless Temp','Inner Temp','Surface Temp')
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

















