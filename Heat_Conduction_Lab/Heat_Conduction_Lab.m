function Heat_Conduction_Lab()
clc
close all
clear all

% default plotting style
fs = 16;    % font size
ms = 10;    % marker size
lw = 2;     % line width

load('SphereData.mat')

for trial = 1:length(sphere_data)
    t = sphere_data(trial).time;    % Time (sec)
    %T = s.data(:,2:11);  % Disc temperatures (C)
    T = [sphere_data(trial).ss_cent,sphere_data(trial).ss_surf,sphere_data(trial).cu_surf,sphere_data(trial).cu_cent];    % Sphere temperatures (C)
    
    %%% Plot the the data
    figure
    T_i = 60;
    T_inf = 29.2;
    theta = (T - T_inf)./(T_i - T_inf);
    semilogy(t,theta,'linewidth',lw)
    legend('Steel center','Steel surface','Copper center','Copper surface')
    set(gca,'FontSize',fs-2);
    xlabel('Time (sec)','Fontsize',fs)
    ylabel('Temperature (\circC)','Fontsize',fs)
    plot_title = sprintf('Sphere Temperatures : Trial %d',trial);
    title(plot_title)
    
    
    %%% Fit the copper sphere cooling curves
    figure
    T_inf = 29.2;
    T_i = 60.3;
    r = 0.023875;   % radius (m)
    V = (4/3)*pi*r^3;
    As = 4*pi*r^2;
    h = 2400;
    rho = 8933; % kg/m^3
    cp = 385;   % J/kgK
    tau = rho*cp*V/(h*As);  % timne constant (sec)
    T_copper = T_inf + (T_i-T_inf)*exp(-t/tau);
    
    plot(t(1:10:400),T(1:10:400,3),'ro','linewidth',lw,'markersize',ms)
    hold on
    plot(t(1:10:400),T_copper(1:10:400),'linewidth',lw)
    cu_title = sprintf('Copper sphere cooling : Trial %d',trial);
    title(cu_title,'Fontsize',fs)
    xlabel('Time (sec)','Fontsize',fs)
    ylabel('Temperature (\circC)','Fontsize',fs)
    legend('Center temperature','Lumped, h = 2400 W/m^2K')
    set(gca,'FontSize',fs-2);
    
    
    %%%
    % Steel sphere cooling curves
    figure
    plot(t(1:10:end),T(1:10:end,1),'ro','linewidth',lw,'markersize',ms)    % plot the center data
    hold on
    plot(t(1:10:end),T(1:10:end,2),'bs','linewidth',lw,'markersize',ms)    % plot the edge data
    
    h = 1250;   % adjustable guess for h
    rho = 8238; % kg/m^3
    cp = 468;   % J/kgK
    k = 13.4;   % W/mK
    alpha = k./(rho*cp);
    Bi = h*r/k;
    r = 0.023875;   % radius (m)
    V = (4/3)*pi*r^3;
    As = 4*pi*r^2;
    tau = rho*cp*V/(h*As);  % lumped time constant
    
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
    
    plot(t,T_center_model,'g','linewidth',lw)   % plot the model
    plot(t,T_surface_model,'g','linewidth',lw)   % plot the model
    % T_steel = T_inf + (T_i-T_inf)*exp(-t/tau); % lumped
    % plot(t(1:10:600),T_steel(1:10:600))
    xlim([0 max(t)])
    ss_title = sprintf('Steel sphere cooling : Trial %d',trial);
    title(ss_title,'Fontsize',fs)
    xlabel('Time (sec)','Fontsize',fs)
    ylabel('Temperature (\circC)','Fontsize',fs)
    legend('Center data','Surface data','Model Fits')
    set(gca,'FontSize',fs-2);
    
end

end





