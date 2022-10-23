%clear all

station = 'DYFAMED';
data_type = 'ctd';
pathway = ['/media/belharet/HD_belharet/Coupling_1D/' station '/'];


winter = {'Dec','Jan','Feb'};
spring = {'Mar','Apr','May'};
summer = {'Jun','Jul','Aug'};
autumn = {'Sep','Oct','Nov'};
%%
filename = [pathway data_type '_data.mat'];

load(filename)

yy = DYFAMED.Year;

ctd_vars = {'temperature', 'salinity', 'oxy'};

   for i_var=1:length(ctd_vars)
         eval(['ctd_data_' ctd_vars{i_var} '= [];'])
   end
   
   dates_ctd_month = {};
   season_ctd = [];

for i = 1:length(yy)
    
    date_ctd{i} = DYFAMED.Data{i}.Date;
    Column = DYFAMED.Data{i}.Column;
    d = DYFAMED.Data{i}.TSO2;
    
    ddd = date_ctd{i};
    for j=1:length(ddd)
        v_ = d{j};
       pressure = v_(:,2); % n * 5 => (datename, pressure, temperature, salinity, oxy[µmol/kg])
       
       for i_var=1:length(ctd_vars)
            var_ = nan(2500,1);
            var_(pressure) = v_(:,2+i_var);
            eval(['ctd_data_' ctd_vars{i_var} '= [ctd_data_' ctd_vars{i_var} ',var_];'])
       end
       
       ddd_ = ddd{j};
       k = strfind(ddd_,'-');
       dates_ctd_month = [dates_ctd_month, ddd_(k(1)+1:k(2)-1)];
       
       if(ismember(dates_ctd_month{end},winter))
            season_ctd = [season_ctd, 1];
       elseif(ismember(dates_ctd_month{end},spring))
            season_ctd = [season_ctd, 2];
       elseif(ismember(dates_ctd_month{end},summer))
            season_ctd = [season_ctd, 3]; 
       elseif(ismember(dates_ctd_month{end},autumn))
            season_ctd = [season_ctd , 4];
       end
       
    end
    
end


%% bottle
data_type = 'bottle';
filename = [pathway data_type '_data.mat'];

load(filename)

toto = DYFAMED_BTLE.date;
dates = cellstr(toto);

[dd, id] = unique(dates);
[~,id_id] = sort(id);
dates_unq = dd(id_id);

% months of the date

for i=1:length(dates_unq)
    dt = dates_unq{i};
    k = strfind(dt,'-');
    dates_month{i} = dt(k(1)+1:k(2)-1);
end

% depth
dep = DYFAMED_BTLE.CTDpres;

bio_var = {'oxygen','nitrat','nitrit','phosphat','silcat','alkali','tcarbn'};



% depth levels from Nemo
load('/media/belharet/HD_belharet/Optimization_admb/outputs/MALASPINA/global/olevel');

% allocate variables
for j = 1:length(bio_var)
    eval(['Bottle_nemo_' bio_var{j} '= [];'])
end


for i=1:length(dates_unq)
    d = dates_unq{i};
    id = find(strcmp(dates,d)==1);
    
    Bottle.depth{i} = dep(id);
    
    %depth following the nemo resolution
    dep_ = repmat(dep(id),1,length(olevel));
    olevel_ = repmat(olevel',length(id),1);
    dif = dep_ - olevel_;
    [~,id_dep_diff] = min(abs(dif),[],2);
    
     Bottle.depth_nemo{i} = id_dep_diff;
     
    for j = 1:length(bio_var)
        
        eval(['y = DYFAMED_BTLE.' bio_var{j} ';'])
        eval(['Bottle.' bio_var{j} '{i} = y(id);'])
        
        % variable with nemo vertical grid
        X = nan(size(olevel));
        X(id_dep_diff) = y(id);
        %eval(['Bottle_nemo.' bio_var{j} '{i} = X;'])
        eval(['Bottle_nemo_' bio_var{j} '= [Bottle_nemo_' bio_var{j} ',X];'])
    end
    
    % month index : 1=winter, 2=spring, 3=summer, 4=autumn
    if(ismember(dates_month{i},winter))
        season(i) = 1;
    elseif(ismember(dates_month{i},spring))
        season(i) = 2;
    elseif(ismember(dates_month{i},summer))
        season(i) = 3; 
    elseif(ismember(dates_month{i},autumn))
        season(i) = 4;
    end
    
end

%% mean profiles by season for bottles

for i=1:4
   id = find(season == i);
   for j = 1:length(bio_var)
        eval([bio_var{j} '_season(:,i) = nanmean(Bottle_nemo_' bio_var{j} '(:,id),2);'])
        eval([bio_var{j} '_season_std(:,i) = nanstd(Bottle_nemo_' bio_var{j} '(:,id),1,2);'])
        eval([bio_var{j} '_season_q1(:,i) = quantile(Bottle_nemo_' bio_var{j} '(:,id),0.25,2);'])
        eval([bio_var{j} '_season_q3(:,i) = quantile(Bottle_nemo_' bio_var{j} '(:,id),0.75,2);'])
   end
end

% mean profiles by month for bottles
months = {'Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec'};
for i=1:12
    mth = months{i};
   id = find(strcmp(dates_month,mth) == 1);
   for j = 1:length(bio_var)
        eval([bio_var{j} '_month(:,i) = nanmean(Bottle_nemo_' bio_var{j} '(:,id),2);'])
        eval([bio_var{j} '_month_std(:,i) = nanstd(Bottle_nemo_' bio_var{j} '(:,id),1,2);'])
        eval([bio_var{j} '_month_q1(:,i) = quantile(Bottle_nemo_' bio_var{j} '(:,id),0.25,2);'])
        eval([bio_var{j} '_month_q3(:,i) = quantile(Bottle_nemo_' bio_var{j} '(:,id),0.75,2);'])
   end
   
   
end

%% mean profiles by season for ctd

for i=1:4
   id = find(season_ctd == i);
   for j = 1:length(ctd_vars)
        eval([ctd_vars{j} '_ctd_season(:,i) = nanmean(ctd_data_' ctd_vars{j} '(:,id),2);'])
        eval([ctd_vars{j} '_ctd_season_std(:,i) = nanstd(ctd_data_' ctd_vars{j} '(:,id),1,2);'])
        eval([ctd_vars{j} '_ctd_season_q1(:,i) = quantile(ctd_data_' ctd_vars{j} '(:,id),0.25,2);'])
        eval([ctd_vars{j} '_ctd_season_q3(:,i) = quantile(ctd_data_' ctd_vars{j} '(:,id),0.75,2);'])
   end
end


%% time-serie of averaged chl-a from satellite (30/12/1999 -> 31/12/2020) -- source: copernicus

chl = ncread('/media/belharet/HD_belharet/Coupling_1D/DYFAMED/dataset-oc-glo-bio-multi-l4-chl_interpolated_4km_daily-rep_1653494106403.nc','CHL');
chl_m = squeeze(nanmean(nanmean(chl),2));

t = ncread('/media/belharet/HD_belharet/Coupling_1D/DYFAMED/dataset-oc-glo-bio-multi-l4-chl_interpolated_4km_daily-rep_1653494106403.nc','time');
toto = datevec(double(t));
toto = toto(:,1:3);

toto_unq = unique(toto(:,1));
 
CHL =[];
for i =1:length(toto_unq)
    chl_an = nan(12,31);
    
    id = find(toto(:,1)==toto_unq(i));
    
    chl_bis = chl_m(id);
    
    toto_bis = toto(id,:);
    
    for j=1:12
       id_ = find(toto_bis(:,2)== j);
       id_day = toto_bis(id_,3);
       chl_an(j,id_day) = chl_bis(id_);
    end
    
    CHL = cat(3,CHL,chl_an);
    
end

CHL_m = nanmean(CHL,3);

CHL_m_ = reshape(CHL_m',1,numel(CHL_m));

CHL_m_(isnan(CHL_m_)==1) = [];

