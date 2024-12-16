

tspan=[0 200];
initial_conditions = [10 0];

[t,sol] = ode15s(@(t, vars) odesystem(initial_conditions), tspan, initial_conditions);

figure(1)
subplot(1,2,1)
plot(t,sol(:,2))
subplot(1,2,2)
plot(t,sol(:,1))


function sol = odesystem(initial_conditions)

 x = initial_conditions(1);  % Biomass
 p = initial_conditions(2);  % Penicillin

dxdt = -7.342500006*10^8*x*(x - 50.64935065)/((5.532065218*10^11*x - 2.694875777*10^12)*(0.15*x - 7.980978267*10^9*x*(x - 50.64935065)/(5.532065218*10^11*x - 2.694875777*10^12)));
dpdt = -3.990489134*10^7*x^2*(x - 50.64935065)/((5.532065218*10^11*x - 2.694875777*10^12)*(0.0002 - 7.980978267*10^9*x*(x - 50.64935065)*(1 - 7.980978267*10^9*x*(x - 50.64935065)/((5.532065218*10^11*x - 2.694875777*10^12)*Kl))/(5.532065218*10^11*x - 2.694875777*10^12))) - 0.04*p;

sol = [dxdt ; dpdt];
end