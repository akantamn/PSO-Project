%function [Costs_Struct]= ProjectFinal(RenewPercent)

clear all

load LoadProfileData

addpath(genpath(pwd))
%%
%Converting 365 X 24 matrices into hourly vectors
% In our simulation, we have 700 residences at each bus, 100 Commercial
% units and two Industries
Residence1Row=reshape((700*Residence1)',[],1);
Residence2Row=reshape((700*Residence2)',[],1);
Residence3Row=reshape((700*Residence3)',[],1);
CommercialRow=reshape((100*Commercial)',[],1);
IndustrialRow=reshape((2*Industrial)',[],1);

ConsolidatedRow=Residence1Row+Residence2Row+Residence3Row+CommercialRow+IndustrialRow;
temp=[Residence1Row Residence2Row Residence3Row CommercialRow IndustrialRow];
%Maximum solar and wind generation capacity = total energy of all
%Residenca1 and Recidence2 respectively
SolarRow=reshape(700*Solar24',[],1);
WindRow=reshape(700*Wind24',[],1);

%% Begin Matpower Simulation
%
 clear GenAlloc GenSum results Cost mpc
 
mpc=loadcase('GenPlusRenew'); %Loads the cost functions, generator branch and load data from file
RenewPercent=1; %What percentage of residents install net metered renewable?

%running only for first 10 days, modify the iteration count to run for
%whole year
for i=3811:4192 %size(IndustrialRow,1)
    %Modify load at each bus according to load profiles
    mpc.bus(:,3) = [Residence1Row(i,1);Residence2Row(i,1);Residence3Row(i,1);CommercialRow(i,1);IndustrialRow(i,1)];
    %Modify max gen capacity for each generator based on available Solar
    %and Wind Potential. Oil peaker and Coal base load are fixed
%     if sum(mpc.bus(:,3))-sum([RenewPercent*SolarRow(i);RenewPercent*WindRow(i);20000])<= 18000
%         mpc.gen(:,8)=[1;1;1;0];
%     end
    mpc.gen(:,9)=[RenewPercent*SolarRow(i);RenewPercent*WindRow(i);50000;20000];
    %save the results of dcopf in a struct
            mpc.gen(:,8)=[1;1;1;1];

    results(i,1)=dcopf(mpc);
    
    if results(i,1).gen(4,2)<18000
        clear results(i,1)
    mpc.bus(:,3) = [Residence1Row(i,1);Residence2Row(i,1);2*Residence3Row(i,1);CommercialRow(i,1);IndustrialRow(i,1)];
        results(i,1)=dcopf(mpc);
        chargepoint(i)=20000;
    end
        
    %Results of DCOPF Gen Allocation
    GenAlloc(i,:)=results(i,1).gen(:,2);
    %Load at this instant
    GenSum(i)=sum(GenAlloc(i,:));
    %Cost(i,:)=GenCost(GenAlloc(i),[5],[40], [5.1]);
    %Total operation cost
    Cost(i,1)=results(i,1).f;
    i
end

 
fprintf('\n Total Cost of Operation is \t $%0.2f', sum(Cost))
fprintf('\n Cost per kWh is \t $%f', sum(Cost)/(10*sum(GenSum)))
fprintf('\n Annual Peak Load is \t %0.2f MW',max(GenSum)/1000)
%fprintf('\n Peak Load is \t %0.2f MW Solar \t %0.2f MW Gen  ',)

fprintf('\n Energy Generated is \t %0.2f Solar \t %0.2f Wind %0.2f Gas %0.2f Coal',sum(GenAlloc)/1000 ) ...
    %sum(Cost),sum(Cost)/(1000*sum(GenSum)),max(GenAlloc(:,1)),max(GenAlloc(:,2)),mean(GenAlloc)*100./max(GenAlloc))
 %clear GenAlloc Gensum results Cost
 
 tempAlloc(:,:)=GenAlloc(3811:4192,1:4);
tempAlloc(:,5)=GenSum(3811:4192);
 
 figure;
 plot(GenSum,'LineWidth',5)
 hold all
 plot(GenAlloc,'LineWidth',5)
 title(sprintf('The percentage of Renewables is %f',RenewPercent))
 legend('Load','Solar','Wind','Gas','Coal')
 xlabel('Hours')
 ylabel('MW')
 xlim([3811,4192])
 plot(chargepoint,'s','MarkerSize',5,'MarkerFaceColor','b')
 