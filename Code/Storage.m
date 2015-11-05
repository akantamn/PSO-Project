%function [Costs_Struct]= ProjectFinal(RenewPercent)

clear all
clc

%%
% firsthour = 3811; % Day 158, hour 19
% lasthour = 4192; % Day 174, hour 16; 15 days, 21 hours
firsthour = 1;
lasthour = 8760;
coalunit = 5; % Bus 5
% Battery 2 MW, 2 MWhr
% batterypower = 2000; % MW/1000
batterycapacity = 2000; % MWhr/1000
batterycharge = batterycapacity;
batterygen = 0;
batterytuning = 0.6;
timestep = 1; % Number of hours in timestep (for battery)
%curtailment = 'off';
curtailment = 'on';

%%
load LoadProfileData

addpath(genpath(pwd))


%%
% Converting 365 X 24 matrices into hourly vectors
% In our simulation, we have 700 residences at each bus, 100 Commercial
% units and two Industries
Residence1Row = reshape((700*Residence1)',[],1);
Residence2Row = reshape((700*Residence2)',[],1);
Residence3Row = reshape((700*Residence3)',[],1);
CommercialRow = reshape((100*Commercial)',[],1);
IndustrialRow = reshape((2*Industrial)',[],1);

ConsolidatedRow = Residence1Row+Residence2Row+Residence3Row+CommercialRow+IndustrialRow;
temp = [Residence1Row Residence2Row Residence3Row CommercialRow IndustrialRow];
%Maximum solar and wind generation capacity = total energy of all
%Residenca1 and Recidence2 respectively
SolarRow = reshape(700*Solar24',[],1);
WindRow = reshape(700*Wind24',[],1);


% Set range by size of data, unless specified as smaller than dataset
if and(exist('firsthour'), exist('lasthour'))
    temp1 = firsthour; temp2 = lasthour;
end
firsthour = 1;
lasthour = size(IndustrialRow,1);
if and(exist('firsthour'), exist('lasthour'))
    firsthour = temp1; lasthour = temp2;
end


%% Begin Matpower Simulation
%
clear GenAlloc GenSum results Cost mpc

mpc=loadcase('GenPlusRenew2'); %Loads the cost functions, generator branch and load data from file
RenewPercent=1; % Percentage Residential Net Metering


%running only for first 10 days, modify the iteration count to run for
%whole year
for i=firsthour:lasthour %size(IndustrialRow,1)
    i
    %Modify load at each bus according to load profiles
    % Pd
    mpc.bus(:,3) = [Residence1Row(i,1); Residence2Row(i,1); Residence3Row(i,1); CommercialRow(i,1); IndustrialRow(i,1)];
    %Modify max gen capacity for each generator based on available Solar
    %and Wind Potential. Oil peaker and Coal base load are fixed
%     if sum(mpc.bus(:,3))-sum([RenewPercent*SolarRow(i);RenewPercent*WindRow(i);20000])<= 18000
%         mpc.gen(:,8)=[1;1;1;0];
%     end
%    mpc.gen(:,9)=[RenewPercent*SolarRow(i);RenewPercent*WindRow(i);50000;20000];
    
    % Battery method simply alternates charge, discharge
%     if batterycharge > 0
%         charging = 0; % Discharging (Generation)
%     elseif batterycharge <= 0
%         charging = 1; % Charging (Load)
%     elseif batterycharge == -1
%         charging = 2; % Neither discharge nor charge
%     end
    
    %Battery method charges when renewables ramp up, vice versa
    if i==1
        prevgen(1) = mpc.gen(1,9);
        prevgen(2) = mpc.gen(2,9);
    end
    
    if (RenewPercent*SolarRow(i)+RenewPercent*WindRow(i)) < (prevgen(1)+prevgen(2)) % Renewables ramping down
        charging = 0; % Discharging (Generation)
    elseif (RenewPercent*SolarRow(i)+RenewPercent*WindRow(i)) > (prevgen(1)+prevgen(2)) % Renewables ramping up
        charging = 1; % Charging (Load)
    else
        charging = 2; % Neither discharge nor charge
    end
    charging;
    
    if charging == 0
        batterygen = batterycharge*batterytuning;
    elseif charging == 1
        batterygen = 0;
        mpc.bus(3,3) = mpc.bus(3,3) + batterycapacity*batterytuning;
    elseif charging == 2
        batterygen = 0;
    end
    
    % Pmax
    mpc.gen(:,9)=[RenewPercent*SolarRow(i);RenewPercent*WindRow(i); batterygen; 50000; 20000];
    % Save the results of dcopf in a struct
        % Status
        mpc.gen(:,8)=[1; 1; 1; 1; 1];
    results(i,1)=dcopf(mpc);
    
    if strcmp(curtailment, 'on')
%        if results(i,1).gen(coalunit,2)<18000
        if results(i,1).gen(coalunit,2) < 18000
            clear results(i,1);
            mpc.gen(:,8)=[0; 0; 0; 1; 1];
            results(i,1)=dcopf(mpc);
        end
    end
    
    if charging == 0
        batterycharge = batterycharge - batterycharge*batterytuning;
    elseif charging == 1
        batterycharge = batterycharge + batterycapacity*batterytuning;
    elseif charging == 2
    end
    
    results(i,1).gen(:,2)
%     prevgen(1) = mpc.bus(1,3);
%     prevgen(2) = mpc.bus(2,3);
    prevgen(1) = results(i,1).gen(1,2);
    prevgen(2) = results(i,1).gen(2,2);
    
    %Results of DCOPF Gen Allocation
    GenAlloc(i,:)=results(i,1).gen(:,2);
    %Load at this instant
    GenSum(i)=sum(GenAlloc(i,:));
    %Cost(i,:)=GenCost(GenAlloc(i),[5],[40], [5.1]);
    %Total operation cost
    Cost(i,1)=results(i,1).f;
end

 
fprintf('\n Total Cost of Operation is \t $%0.2f', sum(Cost))
fprintf('\n Cost per kWh is \t $%f', sum(Cost)/(10*sum(GenSum)))
fprintf('\n Annual Peak Load is \t %0.2f MW',max(GenSum)/1000)
%fprintf('\n Peak Load is \t %0.2f MW Solar \t %0.2f MW Gen  ',)

fprintf('\n Energy Generated is \t %0.2f Solar \t %0.2f Wind %0.2f Storage %0.2f Gas %0.2f Coal',sum(GenAlloc)/1000 ) ...
    %sum(Cost),sum(Cost)/(1000*sum(GenSum)),max(GenAlloc(:,1)),max(GenAlloc(:,2)),mean(GenAlloc)*100./max(GenAlloc))
%clear GenAlloc Gensum results Cost

tempAlloc(:,:)=GenAlloc(firsthour:lasthour,1:5);
tempAlloc(:,5)=GenSum(firsthour:lasthour);

figure;
plot(GenSum,'LineWidth',1)
hold all
plot(GenAlloc,'LineWidth',1)
title(sprintf('The percentage of Renewables is %f',RenewPercent))
legend('Load','Solar','Wind','Storage','Gas','Coal')
xlabel('Hours')
ylabel('MW')
xlim([firsthour,lasthour])





% LoadProfileData.mat
% Residence1Row, Residence2Row, Residence3Row, CommercialRow, IndustrialRow
%   ConsolidatedRow is the sum
%   temp holds 5 above rows
% SolarRow, WindRow
%   Loads the cost functions, generator branch and load data from file
% mpc = loadcase('GenPlusRenew')
% Percentage of residents that install net-metered generation: RenewPercent=1

% 3811:4192
% Starts at day 158, hour
% Size 381

% Output
%   Total Cost of Operation = $8M
%   Cost/kWh = $0.09/kWh
%   Annual Peak Load = 32 MW
%   Energy Generated = 800 Solar, 950 Wind, 750 Gas, 7000 Coal
%    total 10046


