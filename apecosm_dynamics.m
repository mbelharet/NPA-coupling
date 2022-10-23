Version = 'one-way'; %'two-ways' ; %

suffix = '';%'_mod';

pathway_NPA = ['/home/belharet/Bureau/APECOSM_COUPLED/mokrane-coupling/' Version '/'];

% les classes de taille
length_ = ncread('/media/belharet/HD_belharet/Milestone/orca1_REA_REF_OOPE_Y1958D364.nc','length');
len = length_(:,1);

run('/home/belharet/Bureau/APECOSM_COUPLED/mokrane-coupling/apecosm_1d.m')



%weight_ = ncread('/media/belharet/HD_belharet/Milestone/orca1_REA_REF_OOPE_Y1958D364.nc','weight');
%weight = weight_(:,1);

delt_z = ncread('/media/belharet/HD_belharet/Milestone/corrected_mesh_mask_eORCA1_v2.2.nc','e3t_1d');

%% weight_step
omega_var_MIN = 1.;
omega_var_STEP = 1.;
Lmin = [1.e-5; 1.e-5; 1.e-5];
Lmax = [2.; 2.; 2.];
ALPHA = [1.04; 1.04; 1.04];

omega_var_MAX = omega_var_MIN + (nw - 1) * omega_var_STEP;
BETA = (Lmin - Lmax .* (ALPHA .^ (omega_var_MIN - omega_var_MAX - 1))) ./ (1. - (ALPHA .^ (omega_var_MIN - omega_var_MAX - 1)));
GAM = log((Lmin - Lmax) ./  ((ALPHA .^ omega_var_MIN) - (ALPHA .^ (omega_var_MAX + 1)))) ./ log(ALPHA);

omega_var = omega_var_MIN + [0:nw] * omega_var_STEP; %(1,nw)
tmp_length = repmat(ALPHA,1,nw+1) .^ (repmat(omega_var,ng,1) + repmat(GAM,1,nw+1)) + repmat(BETA,1,nw+1);
tmp_weight = repmat(ALPHA,1,nw+1) .* (tmp_length .^ 3.); %(ng,nw)
weight_step = tmp_weight(:,2:end) - tmp_weight(:,1:end-1) ;

weight = (tmp_weight(:,1:end-1) + tmp_weight(:,2:end)) / 2.;%(ng,nw)
%%


% dn = 2;
% ng = 3;
%nw = nL;
%nz = 75;


% for (size_t g = 0; g < NGROUP; ++g) {
%         food[g] = (repfonct_day[cell][g][0] * C_FONCT[g] / w_pow_c[g][0]) / (1. - repfonct_day[cell][g][0]);
% 
%         for (size_t l = 0; l < N_WEIGHT_CLASS; ++l) {
%             tmp_fonct[g][l] = food[g] / (C_FONCT[g] / w_pow_c[g][l] + food[g]);
%         }
% }
% 
% 
% void compute_oope3d(int cell, size_t dn, ma3d &oope3d_tmp) {
%     int ce = indexDFT2DIF[cell];
%     size_t k_bott = get_bottom(cell);

%% oope

oope = ncread([pathway_NPA 'output_apecosm' suffix '/APECOSM_DYFAMED_FORCED_OOPE_Y2019D004.nc'],'OOPE');
OOPE = squeeze(oope(:,:,2,2,end)); OOPE = permute(OOPE,[2 1]);

delta_z_4d = repmat(delt_z(1:nz),1,dn,ng,nw); delta_z_4d = permute(delta_z_4d,[2 3 4 1]);

profile = permute(profile_norm,[1 2 4 3]);

 OOPE_ = repmat(OOPE,1,1,nz,dn); OOPE_ = permute(OOPE_,[4 1 2 3]) ; % OOPE[g][w] , OOPE_[dn][g][w][z]
 oope3d_tmp = profile .* OOPE_ ./ delta_z_4d; % en joule/kg/m3 % oope3d_tmp[dn][g][w][z]
%oope3d_tmp[g][w][z] = profile[dn][cell][g][w][z] * OOPE[ce][g][w];


%%

CRISTAL_SLOPE = 3;
CRISTAL_CRIT = 100;


% void compute_predation(int cell, double dayfrac, size_t dn) {
% 
% // Schooled biomass of OOPE and low trophic levels
%     double tmp_diat = pow((CRISTAL_CRIT_DIAT * (3. * log(LmaxDIAT / LminDIAT))), CRISTAL_SLOPE);
%     double tmp_micro = pow((CRISTAL_CRIT_MICR * (3. * log(LmaxMICR / LminMICR))), CRISTAL_SLOPE);
%     double tmp_meso = pow((CRISTAL_CRIT_MES * (3. * log(LmaxMES / LminMES))), CRISTAL_SLOPE);
%     double tmp_POM = pow((CRISTAL_CRIT_POM * (3. * log(LmaxPOM / LminPOM))), CRISTAL_SLOPE);
%     double tmp_flg = pow((CRISTAL_CRIT_FLAG * (3. * log(LmaxFLAG / LminFLAG))), CRISTAL_SLOPE);

    % computation of schooling between apecosm species
    %kmax_tmp = min(k_max, k_bott); % k_max est un vecteur de taille ng : k_max(gPrey)
    W = repmat(weight,1,1,dn,nz); W = permute(W,[3 1 2 4]);
    school_tmp = (oope3d_tmp .* W).^ CRISTAL_SLOPE; % oope3d_temp[gPrey][wPrey][z], weight[gPrey][wPrey]

    school = school_tmp ./ (school_tmp + CRISTAL_CRIT .^ CRISTAL_SLOPE); %school[z][gPrey][wPrey]

    %%
    convert_molC_to_joule = 474600; % joule/molC
    CRISTAL_CRIT_FLAG = 1e3;
    CRISTAL_CRIT_DIAT = 1e3;
    CRISTAL_CRIT_MICR = 1e3;
    CRISTAL_CRIT_MES = 1e3;
    LminFLAG = 1.e-6;
    LmaxFLAG = 10.e-6;
    LminDIAT = 10.e-6;
    LmaxDIAT = 100.e-6;
    LminMICR = 20.e-6 ;
    LmaxMICR = 200.e-6;
    LminMES = 200.e-6;
    LmaxMES = 2000.e-6;


%%
    pk_var_name = {'flag','diat','zoo','mes'};
    pk_var_name_upper = upper(pk_var_name); pk_var_name_upper{3} = 'MICR';
    pk_var_name_lower = lower(pk_var_name_upper);
    pk_var_name_long = {'flag','diatom','microzoo','mesozoo'};

    
%%
    % feeding on plankton
    filename_pisces = [pathway_NPA 'output_pisces/DYFAMED_1d_20000101_20191231_ptrc_T.nc'];
    phy_1 = ncread(filename_pisces,'PHY'); flag = squeeze(phy_1(2,2,1:nz,1)) *1e-3 * convert_molC_to_joule ; % (molC/m3)*(joule/molC)  --> joule/m3
    phy_2 = ncread(filename_pisces,'PHY2'); diat = squeeze(phy_2(2,2,1:nz,1)) *1e-3 * convert_molC_to_joule ;% (molC/m3)*(joule/molC)  --> joule/m3
    zoo_1 = ncread(filename_pisces,'ZOO'); zoo = squeeze(zoo_1(2,2,1:nz,1)) *1e-3 * convert_molC_to_joule;% (molC/m3)*(joule/molC)  --> joule/m3
    zoo_2 = ncread(filename_pisces,'ZOO2'); mes = squeeze(zoo_2(2,2,1:nz,1)) *1e-3 * convert_molC_to_joule;% (molC/m3)*(joule/molC)  --> joule/m3
  %%  
    for i=1:length(pk_var_name)
        eval(['tmp =' pk_var_name{i} '.^ CRISTAL_SLOPE;'])
        eval(['tmp_' pk_var_name{i} '= (CRISTAL_CRIT_' pk_var_name_upper{i} '* (3. * log(Lmax' pk_var_name_upper{i} '/ Lmin' pk_var_name_upper{i} '))) .^ CRISTAL_SLOPE;'])
        eval([pk_var_name_long{i} '_school = (tmp ./ (tmp + tmp_' pk_var_name{i} '));'])

        % dl
        eval(['dl_' pk_var_name_lower{i} '= (Lmax' pk_var_name_upper{i} '- Lmin' pk_var_name_upper{i} ')/nw;'])
        eval(['len_prey_' pk_var_name_lower{i} '= Lmin' pk_var_name_upper{i} '+ [0.5:nw] * dl_' pk_var_name_lower{i} ';'])
    end

    
%%
    diat_ = repmat(diatom_school .* diat , 1, ng, nw); diat_ = permute(diat_,[2 3 1]);  % diatom_school[z] .* diat[z]
    microzoo_ = repmat(microzoo_school .* zoo , 1, ng, nw); microzoo_ = permute(microzoo_,[2 3 1]); % microzoo_school[z] * zoo[z]
    mesozoo_ = repmat(mesozoo_school .* mes, 1, ng, nw); mesozoo_ = permute(mesozoo_,[2 3 1]); % mesozoo_school[z] * meso[z]
    flag_ = repmat(flag_school .* flag, 1, ng, nw); flag_ = permute(flag_,[2 3 1]); % flag_school[z] * flg[z]


    %%

    a1 = 5;
    a2 = 0.2;
    rho1 = 2.5;
    rho2 = 10;
     
    for i=1:length(pk_var_name)
        eval(['select_' pk_var_name_long{i} '_ = selectivity(len,len_prey_' pk_var_name_lower{i} ',a1,a2,rho1,rho2);']) %(Npred,Nprey)
        eval(['len_prey_2d = repmat(len_prey_' pk_var_name_lower{i} ',nw,1);'])
        eval(['select_' pk_var_name_long{i} '= (sum(select_' pk_var_name_long{i} '_ .* dl_' pk_var_name_lower{i} './len_prey_2d ,2)) / (log(Lmax' pk_var_name_upper{i} ') - log(Lmin' pk_var_name_upper{i} '));'])  % (Npred,1) 
        
        eval(['select_' pk_var_name_long{i} '= repmat(select_' pk_var_name_long{i} ''',ng,1);']) % (ng,nw)
    end

    tmp1_pk_ = repmat(select_diatom,1,1,nz) .* diat_ + repmat(select_microzoo,1,1,nz) .* mesozoo_ + repmat(select_mesozoo,1,1,nz) .* mesozoo_ + ...
        repmat(select_flag,1,1,nz) .* flag_; % select_diatom[gPred][wPred] --> tmp1_pk_[ng][nw][nz]

    tmp1_pk = repmat(tmp1_pk_,1,1,1,dn); tmp1_pk = permute(tmp1_pk,[4 1 2 3]); %tmp1_pk[dn][ng][nw][nz]
    %%
    % feeding on fish

    selectiv = selectivity(len,len,a1,a2,rho1,rho2); % selectiv[wPred][wPrey]
    selectiv_ = repmat(selectiv,1,1,dn,ng,ng,nz); selectiv_ = permute(selectiv_,[3 4 5 1 2 6]);% selectiv_[dn][gPred][gPrey][wPred][wPrey][z] 
    school_ = repmat(school,1,1,1,1,ng,nw); school_ = permute(school_,[1 2 5 3 6 4]); %school[z][gPrey][wPrey]
    oope3d_temp_ = repmat(oope3d_tmp,1,1,1,1,ng,nw); oope3d_temp_ = permute(oope3d_temp_,[1 2 5 3 6 4]);%oope3d_temp[gPrey][wPrey][z]
    weight_step_ = repmat(weight_step,1,1,dn,ng,nw,nz); weight_step_ = permute(weight_step_,[3 1 4 2 5 6]);
    
    
    %%
    % tmp1[ng][nw][nz]
    tmp1_  = squeeze(sum(sum(selectiv_ .* school_ .* oope3d_temp_ .* weight_step_,3),5)) ; % on somme sur gprey, puis sur wprey 
    tmp1  = tmp1_pk + tmp1_;
    %tmp1 += selectiv[gPred][gPrey][wPred][wPrey] * school[z][gPrey][wPrey] * oope3d_temp[gPrey][wPrey][z] * weight_step[gPrey][wPrey];
    %%
    C_FONCT = [1e-5; 1e-5; 1e-5];
    Ta_LIM = 200000.;
    T_inf = [268.15 ; 268.15 ; 268.15]; T_inf_ = repmat(T_inf,1,dn,nw,nz); T_inf_ = permute(T_inf_,[2 1 3 4]); %(dn,ng,nw,nz)
    T_sup = [280.15 ; 280.15 ; 288.15]; T_sup_ = repmat(T_inf,1,dn,nw,nz); T_sup_ = permute(T_sup_,[2 1 3 4]);
    temper_ = repmat(temper,1,dn,ng,nw); temper_ = permute(temper_,[2 3 4 1]);
    
    C_FONCT_ = repmat(C_FONCT,1,nw,nz,dn); C_FONCT_ = permute(C_FONCT_,[4 1 2 3]); % FONCT[gPred] , C_FONCT_[dn][gPred][wPred][z]
    light_pred = permute(light_pref,[1 2 4 3]); % light_pred [dn][gPred][wPred][z] 
    tcor_tmp = repmat(Tcor',1,dn,ng,nw); tcor_tmp = permute(tcor_tmp,[2,3,4,1]); % tcor[z] , tcor_tmp[dn][gPred][wPred][z]
    %tcor_tmp = 1;

    tlim = 1;
    %tlim = (1. / (1. + exp(Ta_LIM ./ temper_ - Ta_LIM ./ T_inf_))) ./ (1. + exp(Ta_LIM ./ T_sup_ - Ta_LIM ./ temper_));%(dn,ng,nw,nz)
    %tlim(2,:,:,:) = 1;
    
    %%
    C_FONCT_W_DEP = [0.333; 0.333; 0.333];
    w_pow_c = weight .^ repmat(C_FONCT_W_DEP,1,nw);% (ng,nw)
    w_pow_c_ = repmat(w_pow_c,1,1,nz,dn); w_pow_c_ = permute(w_pow_c_,[4 1 2 3]); % w_pow_c[gPred][wPred] , w_pow_c_[dn][gPred][wPred][z]
    %%

    tmp4 = C_FONCT_ ./(light_pred .* tcor_tmp .* tlim .* w_pow_c_) + tmp1 ; %tmp4[dn][gPred][wPred][z]
    %tmp4[gPred][wPred][z] = (C_FONCT[gPred] / (light_pred[dn][cell][gPred][z] * tcor_tmp * tlim[gPred][z] * w_pow_c[gPred][wPred])) + tmp1;
    %%
    %Calculate functional response with Holling II size-dependence
%     tmp1_ = repmat(tmp,1,1,1,dn); tmp1_ = permute(tmp,[4 1 2 3])
     repfonct3D = tmp1 ./ tmp4; % repfonct3D[dn][gPred][wPred][z]

    %Ratio of ingested energy / Ec[t][f][i]
%     community_omega_ = repmat(community_omega,1,1,nz,dn); community_omega_ = permute(community_omega_,[4 1 2 3]);
%     w_pow_m1_3_ = repmat(w_pow_m1_3,1,1,nz,dn); w_pow_m1_3_ = permute(w_pow_m1_3_,[4 1 2 3]);
%     ingest3D = community_omega_ .* tcor_tmp .* tlim_ .* repfonct3D .* w_pow_m1_3_ .* dayfrac;
    
    
    
    %ingest3D[z] = community_omega[cell][gPred][wPred] * tcor_tmp * tlim[gPred][z] * repfonct3D[z] * w_pow_m1_3[gPred][wPred] * dayfrac;









%     // Loop over the groups, as predators
%     for (size_t gPred = 0; gPred < NGROUP; ++gPred) {
%         // The group has been disabled by parameter ENABLED_GROUPS
% 
% 
%         for (size_t wPred = 0; wPred < N_WEIGHT_CLASS; ++wPred) {
%             double profiletot = 0.;
% 
%             // ********* calculation of the repfonct3D, ingest3D and habitat3D arrays.
%             // Check whether the group is feeding at this time of the day
%             if (feed[dn][gPred]) {
%                 // Loop over the water column at current cell
%                 for (size_t z = 0; z < kmax_pred; ++z) {
%                     // Group gPred is feeding
%                     double tmp1 = 0.;
%                     tmp1 += select_diatom[gPred][wPred] * diatom_school[z] * diat[z];
%                     tmp1 += select_microzoo[gPred][wPred] * microzoo_school[z] * zoo[z];
%                     tmp1 += select_mesozoo[gPred][wPred] * mesozoo_school[z] * meso[z];
%                     tmp1 += select_POM[gPred][wPred] * POM_school[z] * goc[z];
%                     tmp1 += select_flag[gPred][wPred] * flag_school[z] * flg[z];
%                     double tcor_tmp = tcor[z];
%                     // Loop over the groups, as preys
%                     for (size_t gPrey = 0; gPrey < NGROUP; ++gPrey) {
%                         if (interact[dn][gPred][gPrey] > 0.) {
%                             // Lop over the size class of prey group gPrey
%                             for (size_t wPrey = min_selectiv[gPred][gPrey][wPred]; wPrey <= max_selectiv[gPred][gPrey][wPred]; ++wPrey) {
%                                 tmp1 += selectiv[gPred][gPrey][wPred][wPrey] * school[z][gPrey][wPrey] * oope3d_temp[gPrey][wPrey][z] * weight_step[gPrey][wPrey];
%                             } // end of wprey
%                         }     // end of interact
%                     }   
% 
%                     tmp4[gPred][wPred][z] = (C_FONCT[gPred] / (light_pred[dn][cell][gPred][z] * tcor_tmp * tlim[gPred][z] * w_pow_c[gPred][wPred])) + tmp1;
%                     // Calculate functional response with Holling II size-dependence
%                     repfonct3D[z] = tmp1 / tmp4[gPred][wPred][z];
%                     // Ratio of ingested energy / Ec[t][f][i]
%                     ingest3D[z] = community_omega[cell][gPred][wPred] * tcor_tmp * tlim[gPred][z] * repfonct3D[z] * w_pow_m1_3[gPred][wPred] * dayfrac;
%                     // Calculate the vertical habitat
%                     habitat3D[z] = habitat_env[dn][cell][gPred][z] * repfonct3D[z];
%                 }  // end of z
%             }      // end of if feed






