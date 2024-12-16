% volume evaportation

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
Batch=47;
startID = batchStartIndex(Batch);
endID = batchStartIndex(Batch+1)-1;

time = table2array(data(startID:endID,"Time_h_"));

Fs = table2array(data(startID:endID,"SugarFeedRate_Fs_L_h_"));
Fa = table2array(data(startID:endID,"AcidFlowRate_Fa_L_h_"));
Fb = table2array(data(startID:endID,"BaseFlowRate_Fb_L_h_"));
Fw = table2array(data(startID:endID,"WaterForInjection_dilution_Fw_L_h_"));
Fpaa = table2array(data(startID:endID,"PAAFlow_Fpaa_PAAFlow_L_h__"));
Foil = table2array(data(startID:endID,"OilFlow_Foil_L_hr_"));
Fdis = table2array(data(startID:endID,"DumpedBrothFlow_Fremoved_L_h_"));

figure(1)
plot(time,Fa+Fb+Fw+Fpaa,LineWidth=2)
hold on
plot(time,Fdis,LineWidth=2)
hold on
plot(time, Fs,LineWidth=2)
hold on 
plot(time, Foil,LineWidth=2)
legend('F_{other}','F_{dis}','F_{s}','F_{oil}')
title('All thing entering and leaving the vessel')
ylabel('Feed rates [L/h]')
xlabel('Time [h]')

%
volume = table2array(data(startID:endID,"VesselVolume_V_L_"));


%%


Fevp = zeros(1,length(dvdt));
Fevp = Fs/5 + Foil/5 + Fpaa/5 + Fa/5 + Fb/5 + Fw/5 + Fdis/5 -dvdt';

for i=1:length(Fevp)
    
    if Fevp(i)<-60
        Fevp(i) = mean(Fevp(1:150));
    end

end

plot(time,Fevp*5)
title('Evaporation of liquid in tank based on data')
xlabel('Time [h]')
ylabel('F_{evp} [L/h]')
hold on
%plot(volume)
hold on



dvdt = zeros(1,length(Fs));
dvdt(1)=0;

for i=2:length(Fs)
   
    dvdt(i) = volume(i)-volume(i-1);

end



%% approximatin (modleling) v
volumeest = zeros(1,1150);
volumeest(1:125) = 57800.*exp(0.0003652367166.*time(1:125));
volumeest(126:248) = 58292.8+17.7.*time(126:248);
volumeest(249:375)= 31088 + 127.*time(249:375);
volumeest(376:500)=59167.1+51.9.*time(376:500);
volumeest(501:510) = 60625+48.*time(501:510);
volumeest(511:650) = 55612+44.8.*time(511:650);
volumeest(651:660) = 564817-737.3.*time(651:660);
volumeest(661:750) = 51678.8+39.9.*time(661:750);
volumeest(751:761) = 593543-681.6.*time(751:761);
volumeest(762:800) = 59391.8+20.*time(762:800);
volumeest(801:850) = 1210.4+92.9.*time(801:850);
volumeest(851:860) = 693856-721.*time(851:860);
volumeest(861:950) = 36658+43.*time(861:950);
volumeest(951:960) = 775035-733.*time(951:960);
volumeest(961:1049) = 18190.3+54.*time(961:1049);
volumeest(1050:1060) = 744963-637.*time(1050:1060);
volumeest(1061:1150) = 2404+63.*time(1061:1150);



figure()


%volumeest(501:1150) = 80000;
%f = 0.005;
%volumeest(401:1150) =  + 30000*sin(2*pi*f*time(401:1150)) - time(401:1150).*0.2;

plot(volumeest)
hold on 
plot(volume)




%%
dvdt2 = Fs+Foil+Fpaa+Fa+Fb+Fw-60+Fdis;
plot(time,dvdt2,LineWidth=3)
hold on 
plot(time,dvdt*5,LineWidth=3)
legend('dV/dt with estimated evaporation of 60 L/h','dV/dt from data')
xlabel('Time [h]')
ylabel('Change in volume [L/h]')
title('dV/dt estimated from data and from data with approximated evaporation feed rate')
%%

plot(time,Fs+Foil)
title('Feed of substrate, F_{s}+F_{oil}')
ylabel('F_{s}+F_{oil} [L/h]')
xlabel('Time [h]')

