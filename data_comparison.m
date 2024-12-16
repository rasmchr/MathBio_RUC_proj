%Data comparison

clear
close all
clc
data = readtable("clean_fermi_data.csv");


%first row of each batch is found and put into an array batchStartIndex
batchStartIndex = zeros(99,1);
% index of batchStartIndex array
j=1;
% row of table with batch ID to array
BatchID = table2array(data(:,'x2_PATControl_PAT_ref_PATRef_'));


for i=1:(length(BatchID)-1)
    
    if BatchID(i) ~= BatchID(i+1)
        batchStartIndex(j) = i+1;
        j = j+1;
    end

end

% choose batch you want to plot here
Batch = 47;

% theese will be automated
startID = batchStartIndex(Batch);
if Batch==99
    endID=length(data);
else
    endID = batchStartIndex(Batch+1)-1;
end

t = table2array(data(startID:endID,"Time_h_"));
p = table2array(data(startID:endID,"PenicillinConcentration_P_g_L_"));
Sc=table2array(data(startID:endID,"SubstrateConcentration_S_g_L_"));
Bc=table2array(data(startID:endID,"OfflineBiomassConcentratio_X_offline_X_gL___1___"));

% interpolating biomass concentration to have same length
original_vector = Bc(~isnan(Bc));
desired_length=length(p);
new_indices1 = linspace(1, length(original_vector), desired_length);
BcI = interp1(1:length(original_vector), original_vector, new_indices1, 'linear');

figure(1)
%subplot(1,2,1)
plot(t,p,'LineWidth',3)
hold on
plot(t,Sc,'LineWidth',3)
hold on
plot(t,BcI,'LineWidth',3)
xlabel('Time [h]','FontSize',15)
ylabel('Concentration [g/L]','FontSize',15)
title('Concentrations of batch 47 from data')
legend('Penicillin concentration','Substrate concentration','Biomass concentration')
ylim([-1 50])
hold on

%% 4D


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


[t, sol] = ode15s(@(t, vars) ode_system4D( vars,Fs, Foil, Cs, Coil, mux, yxs, kx, ks, mup, kp, K,m), tspan, initial_conditions);

%%
plot(t,sol(:,1),'LineWidth',3)
hold on
plot(t,sol(:,2),'LineWidth',3)
hold on
plot(t,sol(:,3),'LineWidth',3)
legend('Penicillin B47','Substrate B47','Biomass B47','Biomass M4D','Penicillin M4D','Substrate M4D','FontSize',12) 
title('Comparing data to 4D model')


%% 3D


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
x0 = 0.5;  % Initial biomass
p0 = 0.0;  % Initial penicillin
s0 = 1;   % Initial substrate

initial_conditions=[x0 p0 s0];
tspan=[0 230];
Fs = 80; % sugar feed rate
Foil = 30; % oil feed rate


[t, sol] = ode15s(@(t, vars) ode_system3D( vars,Fs, Foil, Cs, Coil, v, mux, yxs, kx, ks, mup, kp, K,m), tspan, initial_conditions);

%%
plot(t,sol(:,1),LineWidth=3)
hold on
plot(t,sol(:,2),LineWidth=3)
hold on
plot(t,sol(:,3),LineWidth=3)
legend('Penicillin B47','Substrate B47','Biomass B47','Biomass, M3D','Penicillin, M3D','Substrate, M3D','Fontsize',12) 
title('Comparing data to 3D model')
xlabel('Time [h]','FontSize',12)
ylabel('Concentration, [g/L]')

%% 2D

% Parameter settings
    mux = 0.092;   % Growth rate of biomass
    mup = 0.005; % growth rate of penicillin
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
tspan=[0 230];
Fs = 80; % sugar feed rate
Foil = 30; % oil feed rate


[t, sol] = ode15s(@(t, vars) ode_system2D( vars,Fs, Foil, Cs, Coil, v, mux, yxs, kx, ks, mup, kp, K,m), tspan, initial_conditions);

%%
plot(t,sol(:,2), LineWidth=3)
hold on
plot(t,sol(:,1),LineWidth=3)
legend('Penicillin B47','Biomass B47','Penicillin M2D','Biomass M2D', fontsize=12) 
title('Comparing data to 2D model',FontSize=12)
xlabel('Time [h]',FontSize=12)
ylabel('Concentration, [g/L]',FontSize=12)






%% difference 

% interpolating solution to have same length
s1=sol(:,1);
s2=sol(:,2);
%s3=sol(:,3);

time = table2array(data(startID:endID,"Time_h_"));


pd = [];
xd=[];
%sd=[];

for i=1:length(sol(:,1))
    % Find the closest value
    [~, idx] = min(abs(time - t(i)));
    pd(i) = p(idx)-s2(i);
    xd(i) = BcI(idx)-s1(i);
    %sd(i) = Sc(idx)-s3(i);

end


subplot(1,2,1)
plot(t,pd,'LineWidth',3);
hold on
%plot(t,sd,'LineWidth',3);
%hold on
plot(t,xd,'LineWidth',3);
title('Difference between data and model')
xlabel('Time, [h]',FontSize=15)
ylabel('Difference in concentration [g/L]',FontSize=15)
ylim([-30,40])
legend('Penicillin','Biomass',fontsize=12)


function dydt = ode_system4D( vars,Fs, Foil, Cs, Coil, mux, yxs, kx, ks, mup, kp, K,m)
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
    dydt = [dXdt; dPdt; dSdt;dVdt];  % Ensure dydt is a column vector
end
function dydt = ode_system3D( vars,Fs, Foil, Cs, Coil, v, mux, yxs, kx, ks, mup, kp, K,m)
    % Extract variables
    x = vars(1);  % Biomass
    p = vars(2);  % Penicillin
    s = vars(3);  % Substrate
 

    % substrate steady state

   
    % Define the non-linear equations with volume steady state assumption
    dSdt = (Cs*Fs+Coil*Foil)/v - x * (mux/yxs * (s/(kx*x+s))*x) - m*x; % substrate
    dXdt = mux * (s / (kx*x + s) )* x ;                          % Biomass equation
    dPdt = mup * s/(kp+s*(1+s/ks)) * x  - K*p;                   % Penicillin equation
    
    
    % Return the derivatives as a column vector
    dydt = [dXdt; dPdt; dSdt];  % Ensure dydt is a column vector
end
function dydt = ode_system2D( vars,Fs, Foil, Cs, Coil, v, mux, yxs, kx, ks, mup, kp, K,m)
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
