% statistics on data
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

% dividing up batches with their starting ID
for i=1:(length(BatchID)-1)
    
    if BatchID(i) ~= BatchID(i+1)
        batchStartIndex(j) = i+1;
        j = j+1;
    end

end

%start and end id of each batch
startID=batchStartIndex;
endID=zeros(99,1);
for Batch=1:99

     % theese will be automated
    if Batch==99
        endID(Batch)=height(data);
    else
        endID(Batch) = batchStartIndex(Batch+1)-1;
    end
    
end







%% all of breaking type (penicilin concentration)
pt1 = zeros(1,1150);
NOt1=1;

% time vector of the first 1150 timestamps
tarr=table2array(data(batchStartIndex(3):(batchStartIndex(4)-1),"Time_h_"));

for Batch = 1:99

    
    parr = table2array(data(startID(Batch):endID(Batch),"PenicillinConcentration_P_g_L_"));


    if length(parr)==1150 && parr(1150)<max(parr)-1
        % penc conc over time
        pt1(NOt1,:) = parr;
        
        % plotting each individual
            figure(1)
            subplot(1,2,1)
            plot(tarr,pt1(NOt1,:),'LineWidth',2)
            hold on
            title('penicilin concentration over time, type 1')
            xlabel('time[h]')
            ylabel('penicilin concentration [g/L]')

            NOt1=NOt1+1;
    
    end
end

%mean and standart deviation
meanarray = mean(pt1);
stdarr = std(pt1);

% plotting with std
figure(2)
errorbar(tarr,meanarray,stdarr,'b')
hold on
plot(tarr,meanarray,'r')
title('Mean of the pencilin concentration and standard deviation')
xlabel('time[h]')
ylabel('penicilin concentration [g/L]')
%% all logistic growth (penicilin concentration) evt MM
% counter
NOt2=1;


% time vector of the first 1150 timestamps
tarr=table2array(data(batchStartIndex(3):(batchStartIndex(4)-1),"Time_h_"));

% all logistic 61-89
for Batch = 1:99

     
    % penc conc over time
    parr = table2array(data(startID(Batch):endID(Batch),"PenicillinConcentration_P_g_L_"));
    
 
    if length(parr)==1150 &&  parr(end)>max(parr)-2
        % putting into matrix
            p(NOt2,:) = parr;
           
            figure(1)
            subplot(1,2,2)
            plot(tarr,p(NOt2,:),'LineWidth',2)
            hold on
            title('penicilin concentration over time, type 2')
            xlabel('time[h]')
            ylabel('penicilin concentration [g/L]')
            NOt2=NOt2+1;
    end
end
%mean and standart deviation
meanarray = mean(p);
stdarr = std(p);

figure(4)
errorbar(tarr,meanarray,stdarr,'b')
hold on
plot(tarr,meanarray,'r')
title('Mean of the pencilin concentration and standart deviation, type 2')
xlabel('time[h]')
ylabel('penicilin concentration [g/L]')



%% biomass concentration
NObc=1;
tbc=[];

% time vector
Bc=table2array(data(startID(1):endID(1),"OfflineBiomassConcentratio_X_offline_X_gL___1___"));
tbc=find(~isnan(Bc));

for Batch=1:99
    % biomass concentration
    Bc=table2array(data(startID(Batch):endID(Batch),"OfflineBiomassConcentratio_X_offline_X_gL___1___"));
        
    

          % allocating space for new vectors for line plot instead of scatter
          Bcn=(Bc(~isnan(Bc)));
        
         if length(Bcn)==21

            BCm(NObc,:)=Bcn;
            figure(5)
            plot(tbc./5,Bcn)
            hold on
            title('Biomass concentration pr time')
            xlabel('time[h]')
            ylabel('biomass concentration [g/L]')


            NObc= NObc+1;

        end

end

%mean and standart deviation
meanarray = mean(BCm);
stdarr = std(BCm);

figure(6)
errorbar(tbc./5,meanarray,stdarr,'b')
hold on
plot(tbc./5,meanarray,'r')
title('Mean of the biomass concentration and standart deviation')
xlabel('time[h]')
ylabel('biomass concentration [g/L]')


%% substrate concentration
NOs=1;
Sc=[];

% time vector of the first 1150 timestamps
tarr=table2array(data(batchStartIndex(3):(batchStartIndex(4)-1),"Time_h_"));


for Batch=1:99
    

    
     sca = table2array(data(startID(Batch):endID(Batch),"SubstrateConcentration_S_g_L_"));
            
        

        if length(sca)==1150

            Sc(NOs,:)=sca;

            figure(7)
            plot(tarr,Sc(NOs,:))
            hold on
            title('substrate concentration pr time')
            xlabel('time [h]')
            ylabel('Substrate concentration [g/L]')
            ylim([0 0.01])
            xlim([0 230])

              NOs=NOs+1;
        end


end
%mean and standart deviation
meanarray = mean(Sc);
stdarr = std(Sc);

figure(8)
errorbar(tarr,meanarray,stdarr,'b')
hold on
plot(tarr,meanarray,'r')
title('Mean of the substrate concentration and standard deviation')
xlabel('time [h]')
ylabel('Substrate concentration [g/L]')
