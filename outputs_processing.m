
clearvars;

Version = {'one-way','two-ways'};%'apecosm-forced'
variable = 'OOPE'; %'FORAGE';% 

% parameter to convert joules into kg
psi = 4; %j/mg

load('/media/belharet/HD_belharet/Optimization_admb/outputs/MALASPINA/global/olevel');
id = find(olevel<=1000);

% les classes de taille
length_ = ncread('/media/belharet/HD_belharet/Milestone/orca1_REA_REF_OOPE_Y1958D364.nc','length');
L = length_(:,1);


for i_ver=1:length(Version)

pathway = ['/home/belharet/Bureau/APECOSM_COUPLED/mokrane-coupling/' Version{i_ver} '/output_apecosm/'];

list_forage = dir([pathway 'APECOSM_DYFAMED_FORCED_FORAGE_Y*.nc']); 

if(isempty(list_forage))
    list_forage = dir([pathway 'APECOSM_DYFAMED_FORCED_FORAGE_Y*.nc.*']);
end

list_oope = dir([pathway 'APECOSM_DYFAMED_FORCED_OOPE_Y*.nc']); 

if(isempty(list_oope))
    list_oope = dir([pathway 'APECOSM_DYFAMED_FORCED_OOPE_Y*.nc.*']);
end


var = [];
var_oope = [];
for i =1:length(list_forage)
    filename = [pathway list_forage(i).name];
    v = ncread(filename,'FORAGE');
    v_ = squeeze(v(:,:,id,2,2,:,:)) / psi; % kg/m3
    var = cat(5,var,v_); %(weight,com,depth,dn,time)
    
    filename = [pathway list_oope(i).name];
    v = ncread(filename,'OOPE');
    v_ = squeeze(v(:,:,2,2,:)) / psi;%(weight,com,time) % kg/m2
    var_oope = cat(3,var_oope,v_);
    
    
end

%%
eval(['epi_' num2str(i_ver) '= squeeze(var(:,1,:,:,:));'])% (weight,depth,dn,time)
eval(['mig_' num2str(i_ver) '= squeeze(var(:,2,:,:,:));'])% (weight,depth,dn,time)
eval(['mes_' num2str(i_ver) '= squeeze(var(:,3,:,:,:));'])% (weight,depth,dn,time)

eval(['epi_' num2str(i_ver) '_oope= squeeze(var_oope(:,1,:));'])% (weight,time)
eval(['mig_' num2str(i_ver) '_oope= squeeze(var_oope(:,2,:));'])% (weight,time)
eval(['mes_' num2str(i_ver) '_oope= squeeze(var_oope(:,3,:));'])% (weight,time)

end

%%
Sizes = [5.e-3, 30.e-3, 250.e-3, 600.e-3, 2.]; % en m
dn = 1;
x_scale = 1:size(var_oope,3);
sz1= 4;%'3';
sz2= 5;
sz1_ = num2str(find(L<=0.1,1,'last'));
sz2_ = num2str(find(L<=0.2,1,'last'));

%%
figure;
subplot(321)
eval(['plot(squeeze(mig_1(1,:,dn,x_scale(end))),-olevel(id))'])
hold on
eval(['plot(squeeze(mig_2(1,:,dn,x_scale(end))),-olevel(id))'])
legend(Version{1},Version{2})
title({['Com : Mig / Size : '  num2str(Sizes(1)) ' m']},'fontsize',6,'fontweight','bold')
grid minor;
set(gca,'fontsize',8,'fontweight','bold')

subplot(323)
eval(['plot(squeeze(mig_1(' num2str(sz1) ',:,dn,x_scale(end))),-olevel(id))'])
hold on
eval(['plot(squeeze(mig_2(' num2str(sz1) ',:,dn,x_scale(end))),-olevel(id))'])
title({['Com : Mig / Size : '  num2str(Sizes(sz1)) ' m']},'fontsize',6,'fontweight','bold')
grid minor;
set(gca,'fontsize',8,'fontweight','bold')

subplot(322)
eval(['plot(squeeze(mes_1(1,:,dn,x_scale(end))),-olevel(id))'])
hold on
eval(['plot(squeeze(mes_2(1,:,dn,x_scale(end))),-olevel(id))'])
title({['Com : Res / Size : '  num2str(Sizes(1)) ' m']},'fontsize',6,'fontweight','bold')
grid minor;
set(gca,'fontsize',8,'fontweight','bold')


subplot(324)
eval(['plot(squeeze(mes_1(' num2str(sz1) ',:,dn,x_scale(end))),-olevel(id))'])
hold on
eval(['plot(squeeze(mes_2(' num2str(sz1) ',:,dn,x_scale(end))),-olevel(id))'])
title({['Com : Res / Size : '  num2str(Sizes(sz1)) ' m']},'fontsize',6,'fontweight','bold')
grid minor;
set(gca,'fontsize',8,'fontweight','bold')

subplot(325)
eval(['plot(squeeze(mig_1(' num2str(sz2) ',:,dn,x_scale(end))),-olevel(id))'])
hold on
eval(['plot(squeeze(mig_2(' num2str(sz2) ',:,dn,x_scale(end))),-olevel(id))'])
title({['Com : Mig / Size : '  num2str(Sizes(sz2)) ' m']},'fontsize',8,'fontweight','bold')
grid minor;
set(gca,'fontsize',8,'fontweight','bold')


subplot(326)
eval(['plot(squeeze(mes_1(' num2str(sz2) ',:,dn,x_scale(end))),-olevel(id))'])
hold on
eval(['plot(squeeze(mes_2(' num2str(sz2) ',:,dn,x_scale(end))),-olevel(id))'])
title({['Com : Res / Size : '  num2str(Sizes(sz2)) ' m']},'fontsize',8,'fontweight','bold')
grid minor;
set(gca,'fontsize',8,'fontweight','bold')

%print('Figures/compare_forced_coupled_vertical_profiles','-dpng')

figure;
subplot(321)
eval(['plot(squeeze(epi_1(1,:,dn,x_scale(end))),-olevel(id))'])
hold on
eval(['plot(squeeze(epi_2(1,:,dn,x_scale(end))),-olevel(id))'])
legend(Version{1},Version{2})
title({['Com : Epi / Size : '  num2str(Sizes(1)) ' m']},'fontsize',6,'fontweight','bold')
grid minor;
set(gca,'fontsize',8,'fontweight','bold')

subplot(323)
eval(['plot(squeeze(epi_1(' num2str(sz1) ',:,dn,x_scale(end))),-olevel(id))'])
hold on
eval(['plot(squeeze(epi_2(' num2str(sz1) ',:,dn,x_scale(end))),-olevel(id))'])
title({['Com : Epi / Size : '  num2str(Sizes(sz1)) ' m']},'fontsize',6,'fontweight','bold')
grid minor;
set(gca,'fontsize',8,'fontweight','bold')

subplot(325)
eval(['plot(squeeze(epi_1(' num2str(sz2) ',:,dn,x_scale(end))),-olevel(id))'])
hold on
eval(['plot(squeeze(epi_2(' num2str(sz2) ',:,dn,x_scale(end))),-olevel(id))'])
title({['Com : Epi / Size : '  num2str(Sizes(sz2)) ' m']},'fontsize',8,'fontweight','bold')
grid minor;
set(gca,'fontsize',8,'fontweight','bold')
%%
%%
com = 'epi';%'mes';%'mig';%
com_lab = 'epi';%'Res'; %'Mig' ;%
x = repmat(x_scale*5/365,length(id),1);
Dep = repmat(-olevel(id),1,length(x_scale));

%%

% surface
figure;
subplot(311)
eval(['plot(x(1,:),squeeze(' com '_1(1,1,dn,x_scale)))'])
hold on
eval(['plot(x(1,:),squeeze(' com '_2(1,1,dn,x_scale)))'])

title({['Com : ' com_lab ' / Size : '  num2str(Sizes(1)) ' m'], 'Surface'},'fontsize',8,'fontweight','bold')
grid minor;
set(gca,'fontsize',8,'fontweight','bold')

subplot(312)
eval(['plot(x(1,:),squeeze(' com '_1(' num2str(sz1) ',1,dn,x_scale)))'])
hold on
eval(['plot(x(1,:),squeeze(' com '_2(' num2str(sz1) ',1,dn,x_scale)))'])

title({['Com : ' com_lab ' / Size : '  num2str(Sizes(sz1)) ' m'], 'Surface'},'fontsize',8,'fontweight','bold')
%print(['Figures/time_series_surface_' com],'-dpng')
grid minor;
set(gca,'fontsize',8,'fontweight','bold')

subplot(313)
eval(['plot(x(1,:),squeeze(' com '_1(' num2str(sz2) ',1,dn,x_scale)))'])
hold on
eval(['plot(x(1,:),squeeze(' com '_2(' num2str(sz2) ',1,dn,x_scale)))'])

title({['Com : ' com_lab ' / Size : '  num2str(Sizes(sz2)) ' m'], 'Surface'},'fontsize',8,'fontweight','bold')

grid minor;
set(gca,'fontsize',8,'fontweight','bold')

%print(['Figures/compare_forced_coupled_time_series_surface_' com],'-dpng')
%%

%profondeur
dep_limit = 1000;
id_ = find(olevel<=dep_limit,1,'last');

figure;
subplot(311)
eval(['plot(x(1,:),squeeze(' com '_1(1,id_,dn,x_scale)))'])
hold on
eval(['plot(x(1,:),squeeze(' com '_2(1,id_,dn,x_scale)))'])
title({['Com : ' com_lab ' / Size : '  num2str(Sizes(1)) ' m'], [num2str(dep_limit)  ' m']},'fontsize',8,'fontweight','bold')
grid minor;
set(gca,'fontsize',8,'fontweight','bold')

subplot(312)
eval(['plot(x(1,:),squeeze(' com '_1(' num2str(sz1) ',id_,dn,x_scale)))'])
hold on
eval(['plot(x(1,:),squeeze(' com '_2(' num2str(sz1) ',id_,dn,x_scale)))'])

title({['Com : ' com_lab  ' / Size : '  num2str(Sizes(sz1)) ' m'], [num2str(dep_limit)  ' m']},'fontsize',8,'fontweight','bold')
grid minor;
set(gca,'fontsize',8,'fontweight','bold')

subplot(313)
eval(['plot(x(1,:),squeeze(' com '_1(' num2str(sz2) ',id_,dn,x_scale)))'])
hold on
eval(['plot(x(1,:),squeeze(' com '_2(' num2str(sz2) ',id_,dn,x_scale)))'])

title({['Com : ' com_lab ' / Size : '  num2str(Sizes(sz2)) ' m'], [num2str(dep_limit)  ' m']},'fontsize',8,'fontweight','bold')

grid minor;
set(gca,'fontsize',8,'fontweight','bold')

%print(['Figures/compare_forced_coupled_time_series_' num2str(dep_limit)  '_m_' com],'-dpng')


%%
%  figure;
%  subplot(211)
%  eval(['pcolor(x,Dep,squeeze(' com '_1(1,:,1,x_scale)));']) 
%  shading flat
%  colorbar
%  title([com ' : w1  Forced'])
%  
%  subplot(212)
%  eval(['pcolor(x,Dep,squeeze(' com '_2(1,:,1,x_scale)));']) 
%  shading flat
%  colorbar
%  title([com ' : w1  one-way'])
 %%
 %x_scale = 1:size(x,2)/2;
 figure;
subplot(311)
eval(['plot(x(1,:), squeeze( ' com '_1_oope(1,x_scale)))'])
hold on
eval(['plot(x(1,:), squeeze(' com '_2_oope(1,x_scale)))'])
title({['Com : ' com_lab ' / Size-class : 1' ], 'OOPE'},'fontsize',6,'fontweight','bold')
grid minor;
set(gca,'fontsize',8,'fontweight','bold')

subplot(312)
eval(['plot(x(1,:), squeeze(' com '_1_oope(' sz1_ ',x_scale)))'])
hold on
eval(['plot(x(1,:), squeeze(' com '_2_oope(' sz1_ ',x_scale)))'])
title({['Com : ' com_lab ' / Size-class : ' sz1_  ], 'OOPE'},'fontsize',6,'fontweight','bold')
grid minor;
set(gca,'fontsize',8,'fontweight','bold')

subplot(313)
eval(['plot(x(1,:), squeeze(' com '_1_oope(' sz2_ ',x_scale)))'])
hold on
eval(['plot(x(1,:), squeeze(' com '_2_oope(' sz2_ ',x_scale)))'])
title({['Com : ' com_lab ' / Size-class : ' sz2_  ], 'OOPE'},'fontsize',6,'fontweight','bold')

grid minor;
set(gca,'fontsize',8,'fontweight','bold')

%print(['Figures/compare_forced_coupled_time_series_oope_' com],'-dpng')

%% spectres de taille


mig_last_year = mig_2_oope(:,end-365+180);
mes_last_year = mes_2_oope(:,end-365+180);
epi_last_year = epi_2_oope(:,end-365+180);

figure;
plot(log10(L),log10(mig_last_year),'linewidth',2)
hold on
plot(log10(L),log10(mes_last_year),'linewidth',2)
plot(log10(L),log10(epi_last_year),'linewidth',2)

legend('Mig','Res','Epi')
xlabel('log_{10}(L)')
ylabel('log_{10}(B)')
