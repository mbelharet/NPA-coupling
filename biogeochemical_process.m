%clearvars;

conf = 'one-way';%'two-ways';%
suffix = '_mod';
pathway = ['/home/belharet/Bureau/APECOSM_COUPLED/mokrane-coupling/' conf '/output_pisces/'];
filaname = 'DYFAMED_1d_20000101_20291231_ptrc_T.nc';

variables = {'O2','GOC','POC','PHY2','PHY','ZOO2','ZOO','PAR','PAR_APE','PAR_APE_NDCY','NO3','NH4','PO4','Alkalini','Si','NCHL','DCHL','TCHL'};


for i = 1:length(variables)
   var = variables{i};
   
   eval([var ' = ncread([pathway filaname],var);']) 
   
   eval([var '= squeeze(' var '(2,2,:,:));'])
   
end

variables_ = {'tem','sal'};
filaname_ = 'DYFAMED_1d_20000101_20191231_grid_T.nc';

for i = 1:length(variables_)
   var = variables_{i};
   
  eval([var ' = ncread([''/home/belharet/Bureau/APECOSM_COUPLED/mokrane-coupling/' conf '/output_pisces' suffix '/'' filaname_],var);']) 
  
   %eval([var ' = ncread([pathway filaname_],var);']) 
   
   eval([var '= squeeze(' var '(2,2,:,:));'])
   
end

variables = [variables, variables_];
%%
load('/media/belharet/HD_belharet/Optimization_admb/outputs/MALASPINA/global/olevel');
id = find(olevel<=200);
depth = olevel(id);

dep = repmat(depth,1,size(O2,2));
x = repmat(1:size(O2,2),length(depth),1);

%% by season

for i = 1:length(variables)
   eval(['v ='  variables{i} ';']);
   
   eval([variables{i} '_season(:,1) = nanmean(v(:,end-365-30:end-365+60),2);'])
   eval([variables{i} '_season(:,2) = nanmean(v(:,end-365+61:end-365+150),2);'])
   eval([variables{i} '_season(:,3) = nanmean(v(:,end-365+151:end-365+240),2);'])
   eval([variables{i} '_season(:,4) = nanmean(v(:,end-365+241:end-365+334),2);'])
   
   eval([variables{i} '_season_std(:,1) = nanstd(v(:,end-365-30:end-365+60),1,2);'])
   eval([variables{i} '_season_std(:,2) = nanstd(v(:,end-365+61:end-365+150),1,2);'])
   eval([variables{i} '_season_std(:,3) = nanstd(v(:,end-365+151:end-365+240),1,2);'])
   eval([variables{i} '_season_std(:,4) = nanstd(v(:,end-365+241:end-365+334),1,2);'])
   
end

%%



% figure;
% subplot(211)
% pcolor(x,-dep,v(id,:)); shading flat;
% colorbar
 

% v = 'NO3';
% 
% id = find(olevel<=2500);
% depth = olevel(id);
% 
% figure;
% eval(['plot(' v '_season(id,:),-depth,''linewidth'',2)'])

