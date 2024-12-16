% project 4D
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


% Initial conditions
x0 = 0.5;  % Initial biomass
p0 = 0.0;  % Initial penicillin
s0 = 1;   % Initial substrate
v0 = 60000; % initial volume 

initial_conditions=[x0 p0 s0 v0];
tspan=[0 230];
Fs = 80; % sugar feed rate
Foil = 30; % oil feed rate


[t, sol] = ode15s(@(t, vars) ode_system( vars,Fs, Foil, Cs, Coil, mux, yxs, kx, ks, mup, kp, K,m), tspan, initial_conditions);



subplot(1,2,1)
plot(t,sol(:,1),LineWidth=3)
hold on
plot(t,sol(:,2),LineWidth=3)
hold on
plot(t,sol(:,3),LineWidth=3)
legend('Biomass concentration','Penicillin concentration','Substrate concentration','Volume',fontsize=14) 
title('System solved with numerical solver')
xlabel('Time [h]',FontSize=14)
ylabel('Concentration, [g/L]',FontSize=14)
ylim([0,55])
%%
subplot(1,2,2)
semilogy(t, sol(:,4),LineWidth=3)
xlabel('Time [h]',FontSize=14)
ylabel('Volume [L]',FontSize=14)
title('Volume of vessel',FontSize=14)



function dydt = ode_system( vars,Fs, Foil, Cs, Coil, mux, yxs, kx, ks, mup, kp, K,m)
    % Extract variables
    x = vars(1);  % Biomass
    p = vars(2);  % Penicillin
    s = vars(3);  % Substrate
    v = vars(4);
 

    % substrate steady state

   
    % Define the non-linear equations with volume steady state assumption
    dVdt = Fs + Foil +24;
    dSdt = (Cs*Fs+Coil*Foil)/v - x * (mux/yxs * (s/(kx*x+s))*x) - m*x-s/v*dVdt; % substrate
    dXdt = mux * (s / (kx*x + s) )* x -x/v *dVdt;                          % Biomass equation
    dPdt = mup * s/(kp+s*(1+s/ks)) * x  - K*p - p/v * dVdt;                   % Penicillin equation
    
    
    % Return the derivatives as a column vector
    dydt = [dXdt; dPdt; dSdt; dVdt];  % Ensure dydt is a column vector
end
