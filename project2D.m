% project 2D
clear
close all
clc
data = readtable("clean_fermi_data.csv");



% Parameter settings
    mux = 0.5;   % Growth rate of biomass
    mup = 0.15; % growth rate of penicillin
    kp=  0.0002; % monod saturation constant
    ks = 0.1; % substrate inhibition
    kx = 0.15; % contois saturation
    yxs = 0.45;     % Biomass yield coefficient
    K = 0.04;  % decay rate of penicilin
    m = 0.022; % consumption of biomass to maintain life without production of penicilin
   
    Cs = 600;  % sugar feed concentration
    Coil = 1000; % oil feed concentration


% Volume stady state assumption
v = 80000;

% Initial conditions
x0 = 0.6;  % Initial biomass
p0 = 0.0;  % Initial penicillin

initial_conditions=[x0 p0];
tspan=[0 1600];
Fs = 80; % sugar feed rate
Foil = 30; % oil feed rate


[t, sol] = ode15s(@(t, vars) ode_system( vars,Fs, Foil, Cs, Coil, v, mux, yxs, kx, ks, mup, kp, K,m), tspan, initial_conditions);


figure(1)
plot(t,sol(:,1),LineWidth=3)
hold on
plot(t,sol(:,2),LineWidth=3)
legend('Biomass concentration','Penicillin concentration',fontsize=14) 
title('2D system, volume and substrate steady state',FontSize=14)
xlabel('Time [h]','FontSize',14)
ylabel('Concentration, [g/L]',FontSize=14)
ylim([0,150])



function dydt = ode_system( vars,Fs, Foil, Cs, Coil, v, mux, yxs, kx, ks, mup, kp, K,m)
    % Extract variables
    x = vars(1);  % Biomass
    p = vars(2);  % Penicillin
 

    % substrate steady state
    s = -kx*x*yxs*(v*m*x - Coil*Foil - Cs*Fs)/((v*m*x - Coil*Foil - Cs*Fs)*yxs + v*mux*x);

   
    % Define the non-linear equations with volume steady state assumption
    dXdt = mux * (s / (kx*x + s) )* x ;                          % Biomass equation
    dPdt = mup * s/(kp+s*(1+s/ks)) * x  - K*p;                   % Penicillin equation
    
    
    % Return the derivatives as a column vector
    dydt = [dXdt; dPdt];  % Ensure dydt is a column vector
end