% parameter estimation

% import cvs
clear
close all
clc

% getting real data to compare
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

% Choosing batch to fit to
Batch=3;
startID = batchStartIndex(Batch);
endID = batchStartIndex(Batch+1)-1;

% real data arrays
p = table2array(data(startID:endID,"PenicillinConcentration_P_g_L_"));
Sc=table2array(data(startID:endID,"SubstrateConcentration_S_g_L_"));
Bc=table2array(data(startID:endID,"OfflineBiomassConcentratio_X_offline_X_gL___1___"));


% interpolating biomass concentration to have same length
original_vector = Bc(~isnan(Bc));
desired_length=length(p);
new_indices = linspace(1, length(original_vector), desired_length);
BcI = interp1(1:length(original_vector), original_vector, new_indices, 'linear');


% normalizing vectors
%pn=p./max(p);
%Scn=Sc./max(Sc);
%Bcn=BcI./max(BcI);

%% parameter estimation

dydt=[];
Feedrate=80;
Vtot = 70000;

Cs_opt = 600;
kc_opt=0.145;
yxs_opt=0.445;
ypx_opt=0.9;
mux_opt = 0.103;
mx_opt = 0.029;
K_opt = 0.04;

Cs = linspace(400,800,10);
kc = linspace(0.05,0.25,10);
yxs = linspace(0.3,0.7,10);
ypx = linspace(0.8,1,10);
mux = linspace(0.02,0.2,10);
mx = linspace(0.01,0.05,10);
K = linspace(0.01,0.1,10);


tspan=[0,230];
% Initial conditions
x0 = 0.1;  % Initial biomass
p0 = 0.0;  % Initial penicillin
s0 = 1;   % Initial substrate
%v0=100;

initial_conditions = [x0,p0,s0];
desired_length = 1150;


error = 100000;

for Cs_i=1:10
    for kc_i=1:10
        for yxs_i=1:10
            for ypx_i=1:10
                for mux_i=1:10
                    for mx_i=1:10
                        for K_i=1:10

                            [t,sol] = ode45(@(t, vars) alt_ode_system(vars, Feedrate, kc(kc_i), ypx(ypx_i), yxs(yxs_i), mux(mux_i), Vtot, ...
                                K(K_i), Cs(Cs_i), mx(mx_i)), tspan, initial_conditions);
                            
                            %interpolating for the right length
                            original_X = sol(:,1);
                            new_indices = linspace(1, length(original_X), desired_length);
                            interpolated_X = interp1(1:length(original_X), original_X, new_indices, 'linear');
            
                            original_P = sol(:,3);
                            new_indices = linspace(1, length(original_P), desired_length);
                            interpolated_P = interp1(1:length(original_P), original_P, new_indices, 'linear');
            
                            original_S = sol(:,3);
                            new_indices = linspace(1, length(original_S), desired_length);
                            interpolated_S = interp1(1:length(original_S), original_S, new_indices, 'linear');
            
                            % error meassure
                            MSE(1)= abs((sum(interpolated_X)-sum(BcI))^2/1150);
                            MSE(2)= abs((sum(interpolated_P)-sum(p))^2/1150);
                            MSE(3)= abs((sum(interpolated_S)-sum(Sc))^2/1150);
                            errorMean = sum(MSE)/3;
            
                            if errorMean < error
                            kc_opt = kc(kc_i);
                            ypx_opt = ypx(ypx_i);
                            yxs_opt = yxs(yxs_i);
                            mux_opt = mux(mux_i);
                            K_opt = K(K_i);
                            Cs_opt = Cs(Cs_i);
                            mx_opt = mx(mx_i);
                            error = errorMean;
                            end
                        end 
                    end
                end
            end
        end
    end
end
%%


[t, sol] = ode45(@(t, vars) alt_ode_system(vars, Feedrate, kc_opt, ypx_opt, yxs_opt, mux_opt, Vtot, K_opt, Cs_opt, mx_opt), tspan, initial_conditions);

figure()
plot(t,sol(:,1))
hold on
plot(t,sol(:,2))
hold on
plot(t,sol(:,3))
legend("biomass","penicilin", "substrate")


%%





function dydt = alt_ode_system(vars, Feedrate, kc, ypx, yxs, mux,Vtot,K,Cs,mx)
    % Extract variables
    x = vars(1);  % Biomass
    p = vars(2);  % Penicillin
    s = vars(3);  % Substrate

    % Find the index in the time_vector that corresponds to the current time 't'
    %[~, idx] = min(abs(time_vector - t));  % Find the closest index for t

    % Directly use the Feedrate and volume at this index
    %F_tot = Feedrate; %(idx)+Fs(idx);  % Get Feedrate at the current time
   
    % contois
    mu=(mux*s)/(kc*x+s);

    
    % Define the non-linear equations
    dXdt = mu *x ;  % Biomass equation
    dPdt = ypx * x  - K*p;    
    dSdt = Cs*Feedrate/Vtot - mu*x*(1/yxs) - mx*x; % Substrate equation
    
    
    % Return the derivatives as a column vector
    dydt = [dXdt; dPdt; dSdt];  % Ensure dydt is a column vector
end
