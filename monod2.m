% Clear workspace and command window
clear; clc;

% Parameters
mu_max = 0.4;    % Maximum specific growth rate (1/h)
K_s = 0.1;       % Half-saturation constant (g/L)
Y_xs = 0.5;      % Biomass yield coefficient (g biomass/g substrate)
Y_px = 0.1;      % Product yield coefficient (g product/g biomass)
F_in = 0.05;     % Substrate feed rate (g/L/h)
V = 10;          % Volume of the reactor (L)

% Initial conditions
X0 = 0.1;        % Initial biomass concentration (g/L)
S0 = 1.0;        % Initial substrate concentration (g/L)
P0 = 0.0;        % Initial product concentration (g/L)
initial_conditions = [X0; S0; P0];  % Vector of initial conditions

% Time span for the simulation
tspan = [0 10];  % Simulate from time 0 to 10 hours

% Call ode45 solver
[t, sol] = ode45(@(t, y) monod_ODE(t, y, mu_max, K_s, Y_xs, Y_px, F_in), tspan, initial_conditions);

% Plot results
figure;

subplot(3, 1, 1);
plot(t, sol(:, 1), 'b', 'LineWidth', 2); % Biomass concentration
xlabel('Time (hours)');
ylabel('Biomass Concentration (g/L)');
title('Biomass Growth Over Time');
grid on;

subplot(3, 1, 2);
plot(t, sol(:, 2), 'r', 'LineWidth', 2); % Substrate concentration
xlabel('Time (hours)');
ylabel('Substrate Concentration (g/L)');
title('Substrate Consumption Over Time');
grid on;

subplot(3, 1, 3);
plot(t, sol(:, 3), 'g', 'LineWidth', 2); % Product concentration
xlabel('Time (hours)');
ylabel('Product Concentration (g/L)');
title('Product Formation Over Time');
grid on;

% Define the Monod ODE system
function dydt = monod_ODE(t, y, mu_max, K_s, Y_xs, Y_px, F_in)
    % Extract variables
    X = y(1);  % Biomass concentration
    S = y(2);  % Substrate concentration
    P = y(3);  % Product concentration

    % Monod equation for specific growth rate
    mu = (mu_max * S) / (K_s + S);  % Specific growth rate

    % ODEs
    dXdt = mu * X;                   % Change in biomass
    dSdt = F_in - (1/Y_xs) * mu * X; % Change in substrate (with feedrate)
    dPdt = Y_px * mu * X;            % Change in product

    % Return derivatives as a coloumn vector
    dydt = [dXdt; dSdt; dPdt];
end
