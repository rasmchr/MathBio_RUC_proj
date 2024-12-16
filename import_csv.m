% import cvs
clear
close all
clc
data = readtable("clean_fermi_data.csv");

%%

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

%% choose batch you want to plot here
Batch = 47;

% theese will be automated
startID = batchStartIndex(Batch);
if Batch==99
    endID=length(data);
else
    endID = batchStartIndex(Batch+1)-1;
end

% volume
v = table2array(data(startID:endID,"VesselVolume_V_L_"));
t = table2array(data(startID:endID,"Time_h_"));

figure(1)
eighty = repmat(80000,height(t),1);
%plot(t,eighty)
%hold on


for i=1:99

    Batch = i;
startID = batchStartIndex(Batch);
if Batch==99
    endID=width(data);
else
    endID = batchStartIndex(Batch+1)-1;
end


v = table2array(data(startID:endID,"VesselVolume_V_L_"));
t = table2array(data(startID:endID,"Time_h_"));

plot(t,v)
title('Volume of the vessel over time for all batches')
xlabel('Time [h]')
ylabel('Volume [L]')
ylim([0,10*10^4])
hold on

end




%%
volume = zeros(1,1150);
volume=volume./99

figure(1)
plot(t,v)
title('Volume of the vessel over time for batch 3')


%% penc conc over time
t = table2array(data(startID:endID,"Time_h_"));
p = table2array(data(startID:endID,"PenicillinConcentration_P_g_L_"));

figure(1)
scatter(t,p)
title('penicilin concentration over time')
xlabel('time[h]')
ylabel('penicilin concentration [g/L]')
%%
%hold on

FR = table2array(data(startID:endID,"SugarFeedRate_Fs_L_h_"));

figure(2)
subplot(3,1,1)
scatter(t,FR)
title('substrate rate')
xlabel('time[h]')
ylabel('substrate rate [L/h]')


% finding the penicilin g/L/t
% length of vector ie. number of timepoints in batch
NOit=endID-startID;
PR=zeros(NOit+1,1);
PR(1)=p(1);
for i=2:NOit

PR(i)=(p(i)-p(i-1))/(t(i)-t(i-1));

end
PR(end)=PR(end-1);

subplot(3,1,2)
plot(t,PR)
title('Penicilin rate')
xlabel('time[h]')
ylabel('penicilin rate [g/L/t]')



% biomass concentration
Bc=table2array(data(startID:endID,"OfflineBiomassConcentratio_X_offline_X_gL___1___"));


subplot(3,1,3)
scatter(t,Bc)
title('Biomass concentration pr time')
xlabel('time[h]')
ylabel('biomass concentration [g/L]')


%% substrate concentration
Sc=table2array(data(startID:endID,"SubstrateConcentration_S_g_L_"));

figure(3)
plot(t,Sc)
title('substrate concentration pr time')
xlabel('time [h]')
ylabel('Substrate concentration [g/L]')


%% substrate feed rate concentration
SRC=FR;


%% Simple model implementation
% batch 24

% allocating space
x=zeros(1,length(t));
pm=zeros(1,length(t));
s=zeros(1,length(t));

%% biomass modelling

% making biomass concentration dimensionless
CC = max(Bc);
BCD = Bc./CC;

% mu found experimentally
mu = 0.06;


% approximating mu
% IV: intial value is 0.5g/L and the carrying capacity is 22.7g/L
IV = 0.5/CC;

% simple model of biomass concentration (dimensionless proportion of
% carrying capacity)
x = (IV.*exp(mu*t))./(1-IV+IV.*exp(mu*t));


figure(4)
scatter(t,BCD)
hold on
plot(t,x)
title('biomass prediction(orange line) and actual data (dots)')




figure(5)
plot(t,Feedrate)
xlabel('time [h]')
ylabel('Feedrate [L/h]')
title('Feedrate modelled and from a single batch data')
hold on
plot(t,FR)

figure(6)
plot(t,volume)
xlabel('time [h]')
ylabel('Volume [L]')
title('Volume in liters pr hour')

%%
%Feedrate = zeros(1,1150) +80;
%Feedrate = Feedrate*100;

% penicilling concentration modelling
% Initial conditions
x0 = 0.1;  % Initial biomass
p0 = 0.0;  % Initial penicillin
s0 = 1;   % Initial substrate

% Parameter settings
    mu = 0.092;   % Growth rate
    ypx =  0.005;  % Penicillin yield coefficient
    yxs = 0.445;  % Biomass yield coefficient
    ks = 0.545;   % Substrate saturation constant
    dc = 0.002; % degredation of penicilin

time_vector = t;
%clear t;
initial_conditions = [x0, p0, s0];

tspan = [0 1150];
toutput = linspace(0,1150, 1150);

% Call ode45 solver
[t, sol] = ode15s(@(t, vars) ode_system(t, vars, time_vector, Feedrate, volume,mu,ypx,yxs,ks), toutput, initial_conditions);

% Optionally plot the results
figure(7)
plot( sol)
%xlim([0 200])
%ylim([-5 50])
xlabel('Time')
ylabel('concentration')
legend('Biomass (x)', 'Penicillin (p)', 'Substrate (s)')

%
%% plotting 1 batch all concentrations

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

figure(1)
plot(t,p)
hold on
plot(t,Sc)
hold on
scatter(t,Bc)
title('Penicillin, substrate and biomass concentration of batch 47')
xlabel('Time [h]')
ylabel('Concentration [g/L]')
legend('Penicillin B47','Substrate B47','Biomass B47')
ylim([-1 45])
hold on

%%
% parameter estimation
pn=p./max(p);
Scn=Sc./max(Sc);
Bcn=(~isnan((Bc./max(Bc))));

% penicilling concentration modelling
% Initial conditions
x0 = 0.1;  % Initial biomass
p0 = 0.0;  % Initial penicillin
s0 = 1;   % Initial substrate
v0=50000;

vars=[x0 p0 s0 v0];

% Parameter settings
    mu = 0.092;   % Growth rate
    ypx =  0.005;  % Penicillin yield coefficient
    yxs = 0.445;  % Biomass yield coefficient
    ks = 0.545;   % Substrate saturation constant
    dc = 0.002; % degredation of penicilin

% Parameter settings
    mu_opt = 0.092;   % Growth rate 0.130
    ypx_opt = 0.0050;  % Penicillin yield coefficient 0.015
    yxs_opt = 0.445;  % Biomass yield coefficient
    ks_opt = 0.545;   % Substrate saturation constant
    dc_opt = 0.002; %degradation costant

%%
% mu = linspace(0.2,0.7,15);
ypx = linspace(0,0.05,25);
yxs = linspace(0.3,0.8,25);
ks = linspace(0,0.8,25);
%dc = linspace(0,0.01,5);
    t0=0;
    tf=230;
    N=1150;
 tspan = linspace(t0,tf,N); % Custom time vector

% Call the ODE solver
options = odeset('RelTol',1e-6,'AbsTol',1e-9);

error=4000;


%%
clc
dydt=[];
Feedrate=zeros(1150,1)+0.4;
Fs=zeros(1150,1)+0.2;

mus=0.1;
mux=0.1;
yxs=0.4;
yps=0.4;
ks=0.2;


tspan=[0,230];
% Initial conditions
x0 = 0.1;  % Initial biomass
p0 = 0.0;  % Initial penicillin
s0 = 1;   % Initial substrate
v0=100;

initial_conditions = [x0,p0,s0,v0];


[t, sol] = ode45(@(t, vars) alt_ode_system(t, vars, time_vector, Feedrate, Fs, ks, ypx, yxs, mus, mux), tspan, initial_conditions);
figure()
%plot(t,pn)
%hold on
plot(sol)


%%


% Define the ODE system
function dydt = ode_systemF(vars,Feedrate,mu,ypx,yxs,ks)
    % Extract variables
    x = vars(1);  % Biomass
    p = vars(2);  % Penicillin
    s = vars(3);  % Substrate
    v = vars(4);

   
    
    % Define the non-linear equations
    dXdt = mu * (s / (ks + s) )* x - (x * Feedrate) / v;  % Biomass equation
    dPdt = ypx * x  - (p * Feedrate) / v;                  % Penicillin equation
    dSdt = 0.02*Feedrate - yxs * x - ((s * Feedrate) / v);  % Substrate equation
    dVdt = Feedrate;
    
    % Return the derivatives as a column vector
    dydt = [dXdt; dPdt; dSdt; dVdt];  % Ensure dydt is a column vector
end

function dydt = ode_system(t, vars, time_vector, Feedrate, volume,mu,ypx,yxs,ks)
    % Extract variables
    x = vars(1);  % Biomass
    p = vars(2);  % Penicillin
    s = vars(3);  % Substrate

    % Find the index in the time_vector that corresponds to the current time 't'
    [~, idx] = min(abs(time_vector - t));  % Find the closest index for t

    % Directly use the Feedrate and volume at this index
    current_Feedrate = Feedrate(idx);  % Get Feedrate at the current time
    current_volume = volume(idx);      % Get volume at the current time

   
    
    % Define the non-linear equations
    dXdt = mu * (s / (ks + s) )* x - (x * current_Feedrate) / current_volume;  % Biomass equation
    dPdt = ypx * x  - (p * current_Feedrate) / current_volume;                  % Penicillin equation
    dSdt = 0.02*current_Feedrate - yxs * x - ((s * current_Feedrate) / current_volume);  % Substrate equation
   
    
    % Return the derivatives as a column vector
    dydt = [dXdt; dPdt; dSdt];  % Ensure dydt is a column vector
end

function dydt = alt_ode_system(t, vars, time_vector, Feedrate, Fs, ks, ypx, yxs, mus, mux)
    % Extract variables
    x = vars(1);  % Biomass
    p = vars(2);  % Penicillin
    s = vars(3);  % Substrate
    v = vars(4);  % Volume

    % Find the index in the time_vector that corresponds to the current time 't'
    [~, idx] = min(abs(time_vector - t));  % Find the closest index for t

    % Directly use the Feedrate and volume at this index
    F_tot = Feedrate(idx)+Fs(idx);  % Get Feedrate at the current time
    %current_volume = volume(idx);      % Get volume at the current time
    
    
   
   
    
    % Define the non-linear equations
    %dXdt = mux*s/(ks+s)*x-(x*current_Feedrate)/current_volume;  % Biomass equation
    %dPdt = mupp*x-K*p-(p*current_Feedrate)/current_volume;               % Penicillin equation
    dXdt = mux *(s/(s+ks) )*x - (x * F_tot) /v;  % Biomass equation
    dPdt = ypx * x  - (p * F_tot) /v;    
    dSdt = Fs - mus*yxs*x - (s*F_tot)/v; % Substrate equation
    dVdt = F_tot + Fs(idx);
    
    % Return the derivatives as a column vector
    dydt = [dXdt; dPdt; dSdt; dVdt];  % Ensure dydt is a column vector
end
