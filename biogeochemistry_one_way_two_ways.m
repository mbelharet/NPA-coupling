
clearvars;

Version = {'one-way','two-ways'};

load('/media/belharet/HD_belharet/Optimization_admb/outputs/MALASPINA/global/olevel');
id = find(olevel<=1000);

vars = {'O2','GOC','POC','PHY','PHY2','ZOO','ZOO2','NO3','NH4','PO4','Si','NCHL','DCHL','TCHL','TPOC'};
units = {'mmol m^{-3}','mmol m^{-3}','mmol m^{-3}','mmol m^{-3}','mmol m^{-3}','mmol m^{-3}','mmol m^{-3}','mmol m^{-3}','mmol m^{-3}','mmol m^{-3}','mmol m^{-3}','mg m^{-3}','mg m^{-3}','mg m^{-3}','mmol m^{-3}'};

DEP = repmat(olevel(id),1,365);
X = repmat(1:365,length(id),1);
for i_var = 4:4%length(vars)

var = vars{i_var};

%two-ways
V = ncread('/home/belharet/Bureau/APECOSM_COUPLED/mokrane-coupling/two-ways/output_pisces/DYFAMED_1d_20000101_20291231_ptrc_T.nc',var);
V_ = squeeze(V(2,2,id,:));
%one-way
V_1 = ncread('/home/belharet/Bureau/APECOSM_COUPLED/mokrane-coupling/one-way/output_pisces/DYFAMED_1d_20000101_20291231_ptrc_T.nc',var);
V_1_ = squeeze(V_1(2,2,id,:));

%%
id_dep_200 = find(olevel<=100);% & olevel<=500);

figure;
plot(nanmean(V_(id_dep_200,end-365+1:end)),'linewidth',2);
hold on
plot(nanmean(V_1_(id_dep_200,end-365+1:end)),'linewidth',2);

legend('two\_ways','one\_way')
title(var)

grid minor
xlabel('Time (days)')
ylabel(units{i_var})
set(gca,'fontsize',8,'fontweight','bold')

%%
figure;
VAR = 100*(V_(:,end-365+1:end) - V_1_(:,end-365+1:end))./V_(:,end-365+1:end);
%VAR = V_(:,end-365+1:end) - V_1_(:,end-365+1:end);

pcolor(X,-DEP,VAR); shading flat; 
colorbar;
title([var ' : (2ways - 1way)'])
xlabel('Time (days)')
set(gca,'fontsize',8,'fontweight','bold')

%%
figure;

pcolor(X,-DEP,V_(:,end-365+1:end)); shading flat; 
colorbar;
title([var ' : (2ways)'])
xlabel('Time (days)')
set(gca,'fontsize',8,'fontweight','bold')

figure;

pcolor(X,-DEP,V_1_(:,end-365+1:end)); shading flat; 
colorbar;
title([var ' : (1way)'])
xlabel('Time (days)')
set(gca,'fontsize',8,'fontweight','bold')

end
%%
% figure;
% pcolor(V_1_(:,end-365+1:end)); shading flat;
% 
% figure;
% pcolor(V_(:,end-365+1:end)); shading flat;
    
    
