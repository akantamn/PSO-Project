load LoadProfileData

addpath(genpath(pwd))
%%
%Converting 365 X 24 matrices into hourly vectors
% In our simulation, we have 700 residences at each bus, 100 Commercial
% units and two Industries
Residence1Row=reshape((700*Residence1/1000)',[],1);
Residence2Row=reshape((700*Residence2/1000)',[],1);
Residence3Row=reshape((700*Residence3/1000)',[],1);
CommercialRow=reshape((100*Commercial/1000)',[],1);
IndustrialRow=reshape((2*Industrial/1000)',[],1);

%Maximum solar and wind generation capacity = total energy of all
%Residenca1 and Recidence2 respectively
SolarRow=reshape(700*Solar24'/1000,[],1);
WindRow=reshape(700*Wind24'/1000,[],1);

%% Begin Matpower Simulation
%
 clear GenAlloc GenSum results Cost
 
mpc=loadcase('GenPlusRenew'); %Loads the cost functions, generator branch and load data from file
RenewPercent=0.1; %What percentage of residents install net metered renewable; 1 = 100%

%running only for first 10 days, modify the iteration count to run for
%whole year
for i=1:size(IndustrialRow,1)
    %Modify load at each bus according to load profiles
    mpc.bus(:,3) = [Residence1Row(i,1);Residence2Row(i,1);Residence3Row(i,1);CommercialRow(i,1);IndustrialRow(i,1)];
    %Modify max gen capacity for each generator based on available Solar
    %and Wind Potential. Oil peaker and Coal base load are fixed
    if sum(mpc.bus(:,3))-sum([RenewPercent*SolarRow(i);RenewPercent*WindRow(i);20000])<= 18000
        mpc.gen(:,8)=[0;0;1;1];
    end
    mpc.gen(:,9)=[RenewPercent*SolarRow(i);RenewPercent*WindRow(i);25;20];
    %save the results of dcopf in a struct
    results(i,1)=dcopf(mpc);
   
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
fprintf('\n Annual Peak Load is \t %0.2f MW',max(GenSum))
%fprintf('\n Peak Load is \t %0.2f MW Solar \t %0.2f MW Gen  ',)

fprintf('\n Energy Generated is \t %0.2f Solar \t %0.2f Wind %0.2f Oil %0.2f Gas',sum(GenAlloc) ) ...
    %sum(Cost),sum(Cost)/(1000*sum(GenSum)),max(GenAlloc(:,1)),max(GenAlloc(:,2)),mean(GenAlloc)*100./max(GenAlloc))
 %clear GenAlloc Gensum results Cost
 
 figure;
 plot(GenSum)
 hold all
 plot(GenAlloc)
 legend('Load','Solar','Wind','Oil','Gas')
 xlabel('Hours')
 ylabel('MW')
 