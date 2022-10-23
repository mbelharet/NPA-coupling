clearvars;


Version = 'HOT';%'KERFIX';%'BATS';%'one-way';%'two-ways';%'apecosm-forced';%%
station = 'HOT';%'KERFIX';%'BATS';
variable = 'OOPE'; %'FORAGE';% 

suffix = '';%'_2w';%'_mod';%
pathway = ['/home/belharet/Bureau/APECOSM_COUPLED/mokrane-coupling/' Version '/'];

% parameter to convert joules into kg
psi = 4; %j/mg
convert_joule_to_molC = 474600; % joule/molC

load('/media/belharet/HD_belharet/Optimization_admb/outputs/MALASPINA/global/olevel');
id = find(olevel<=1000);

delt_z = ncread('/media/belharet/HD_belharet/Milestone/corrected_mesh_mask_eORCA1_v2.2.nc','e3t_1d');

% les classes de taille
length_ = ncread('/media/belharet/HD_belharet/Milestone/orca1_REA_REF_OOPE_Y1958D364.nc','length');
L = length_(:,1);

weight_ = ncread('/media/belharet/HD_belharet/Milestone/orca1_REA_REF_OOPE_Y1958D364.nc','weight');
weight = weight_(:,1);



list_oope = dir([pathway 'output_apecosm' suffix '/APECOSM_' station '_FORCED_OOPE_Y*.nc']); 

if(isempty(list_oope))
    list_oope = dir([pathway 'output_apecosm' suffix '/APECOSM_' station '_FORCED_OOPE_Y*.nc.*']);
end


filename = [pathway 'output_apecosm' suffix '/' list_oope(end).name];
v = ncread(filename,'OOPE');
var_oope = squeeze(v(:,:,2,2,:)) ;%(weight,com,time) % joule/kg/m2
 

epi_oope= squeeze(var_oope(:,1,:));% (weight,time)
mig_oope= squeeze(var_oope(:,2,:));% (weight,time)
mes_oope= squeeze(var_oope(:,3,:));% (weight,time)

epi_oope_m = mean(epi_oope,2); % joule/kg/m2
mig_oope_m = mean(mig_oope,2); % joule/kg/m2
mes_oope_m = mean(mes_oope,2); % joule/kg/m2

%% Plankton
filename_pisces = [pathway 'output_pisces' suffix '/' station '_1d_20000101_20191231_ptrc_T.nc']; 

% l'unité des variables pisces est molC/L
phy_1 = ncread(filename_pisces,'PHY'); phy_1 = squeeze(phy_1(2,2,:,end-365:end)) *1e-3 * convert_joule_to_molC ; % molC/m3 --> joule/kg/m2
phy_2 = ncread(filename_pisces,'PHY2'); phy_2 = squeeze(phy_2(2,2,:,end-365:end)) *1e-3 * convert_joule_to_molC ;% molC/m3 --> joule/kg/m2
zoo_1 = ncread(filename_pisces,'ZOO'); zoo_1 = squeeze(zoo_1(2,2,:,end-365:end)) *1e-3 * convert_joule_to_molC ;% molC/m3 --> joule/kg/m2
zoo_2 = ncread(filename_pisces,'ZOO2'); zoo_2 = squeeze(zoo_2(2,2,:,end-365:end)) *1e-3 * convert_joule_to_molC ;% molC/m3 --> joule/kg/m2

% delta_z en 2d (depth * time)
delt_z_2d = repmat(delt_z,1,size(phy_1,2));

% intégration par z et moyenne temporelle (dernière année de simulation)
phy_1_t = mean(sum(phy_1 .* delt_z_2d));% flagellates in molC/m2
phy_2_t = mean(sum(phy_2 .* delt_z_2d));% diatoms in molC/m2
zoo_1_t = mean(sum(zoo_1 .* delt_z_2d));% microzoopk in molC/m2
zoo_2_t = mean(sum(zoo_2 .* delt_z_2d));% mesozoopk in molC/m2

% size spectrum for plankton communities
phy_2_length_min = 10.e-6; % diatoms
phy_1_length_min = 1.e-6;
zoo_2_length_min = 200.e-6;
zoo_1_length_min = 20.e-6 ;
phy_2_length_max = 100.e-6;
phy_1_length_max = 10.e-6;
zoo_2_length_max = 2000.e-6;
zoo_1_length_max = 200.e-6;

allom_coef = 15;

pk_com = {'phy_1','phy_2','zoo_1','zoo_2'};

figure;
for i=1:length(pk_com)
    
    eval([pk_com{i} '_weight_min = allom_coef * (' pk_com{i} '_length_min ^3);'])
    eval([pk_com{i} '_weight_max = allom_coef * (' pk_com{i} '_length_max ^3);'])


    eval(['p_' pk_com{i} '=' pk_com{i} '_t / (log(' pk_com{i} '_weight_max) - log(' pk_com{i} '_weight_min));'])
    eval(['B_' pk_com{i} '= p_' pk_com{i} '* [1/' pk_com{i} '_weight_min 1/' pk_com{i} '_weight_max];'])
    
    %plot plankton communities
    eval(['x = [log10(' pk_com{i} '_weight_min) log10(' pk_com{i} '_weight_max)];'])
    eval(['y = log10(B_' pk_com{i} ');'])
    plot(x,y,'linewidth',2)
    hold on
end

%id_l = find(L<=1);
id_l = 1:length(L);


%plot fish communities

plot(log10(weight(id_l)),log10(epi_oope_m(id_l)),'linewidth',2)
plot(log10(weight(id_l)),log10(mig_oope_m(id_l)),'linewidth',2)
plot(log10(weight(id_l)),log10(mes_oope_m(id_l)),'linewidth',2)

legend('Flag','Diat','microzoo','mesozoo','Epi','Mig','Res','location','southwest')

%% Réponse fonctionnelle
plot_reponse_fonctionnelle = 0;

if(plot_reponse_fonctionnelle)
    list_repfct = dir([pathway 'output_apecosm/APECOSM_DYFAMED_FORCED_repfonct_day_Y*.nc']);
    filename = list_repfct(end).name;
    var = ncread([pathway 'output_apecosm/' filename],'repfonct_day');
    f = squeeze(var(:,:,2,2,:));

end