clear all

run('dyfamed_data_processing.m')
run('biogeochemical_process.m')

import mlreportgen.report.*
R = Report("biogeochemical_validation","pdf");
open("R")

import mlreportgen.dom.*


message1 = Text("This report shows the vertical distribution of biogeochemical tracers simulated using the fully coupled version of Nemo-Pisces-Apecosm and confronted with field observed data " + ...
    "collected in Dyfamed station.");


add(R,message1);

plot_no3 = 1 ;
plot_NH4 = 1;
plot_po4 = 1;
plot_alca = 0;
plot_oxy = 0;
plot_si = 0;
plot_tem = 0;
plot_sal = 0;
plot_chl = 0;

%clearvars id;

id = find(olevel<=1000);%2500
depth = olevel(id);

depth_ctd = 1:2500;
y_fill_ctd = [depth_ctd,depth_ctd(end:-1:1)];
%% nitrate

y_fill = [-depth',-depth(end:-1:1)'];

var_mod_name = {'NO3','NH4','PO4','Alkalini','O2','Si','tem','sal'};
var_obs_name = {'nitrat','nitrit','phosphat','alkali','oxygen','silcat','temperature_ctd','salinity'};
var_cond = [plot_no3, plot_NH4, plot_po4, plot_alca, plot_oxy, plot_si, plot_tem, plot_sal, plot_chl];
units = {'µmol l^{-1}','µmol l^{-1}','µmol l^{-1}','µmol l^{-1}','µmol l^{-1}','µmol l^{-1}','°C',''};

x_lim_min = [0, 0, 0 , 2500 , 150, 0, 10];%, ,0, , , 5];
x_lim_max = [12, 0.7, 0.6, 2700, 300, 10, 25];%, ,0.55, , , 25];
seasons = {'Winter','Spring','Summer','Autumn'};


compteur =0;
figure;


for i_var = 1:length(var_mod_name)

if(var_cond(i_var))
    compteur = compteur + 1;
   message2 = var_mod_name{i_var};
   T1 = Text(strcat(num2str(compteur),") ", message2));
   T1.Bold = true;
   T1.FontSize = "16pt";
   add(R,T1)

% obs_q1 = nitrat_season_q1(id,1);
% obs_q3 = nitrat_season_q3(id,1);
% x_fill = [obs_q1',obs_q3(end:-1:1)'];
% 
% fill(x_fill,y_fill,col(1,:),'EdgeColor','none','FaceAlpha',0.3)

for i=1:length(seasons)

    message3 = seasons{i};

    if(mod(i,2)==1)
        
        %
        add(R,Figure)

        

        figure;
       

        col = get(gca,'ColorOrder');
        subplot(1,2,1)
    else
        subplot(1,2,2)
    end
    
    
    eval(['var_mod =' var_mod_name{i_var} '_season(id,i);'])
    eval(['var_obs =' var_obs_name{i_var} '_season(id,i);'])
    eval(['var_obs_std =' var_obs_name{i_var} '_season_std(id,i);'])
    
    % 
    eval(['var_q1 = fillmissing(' var_obs_name{i_var} '_season_q1(id,i),''linear'');'])
    eval(['var_q3 = fillmissing(' var_obs_name{i_var} '_season_q3(id,i),''linear'');'])

    fill([var_q1(:,1)',var_q3(end:-1:1,1)'],[-depth',-depth(end:-1:1)'],'b','facealpha',0.2,'edgecolor','none')
    hold on
    plot(var_mod,-depth,'k','linewidth',2) % Model
    plot(var_obs,-depth,'.','markersize',10)
%     hold on
%     errorbar(var_obs,-depth,var_obs_std,'horizontal','.','markersize',6)% obs


    grid minor
    xlabel([var_mod_name{i_var} ' (' units{i_var} ')'],'fontweight','bold','fontsize',10)
    title(seasons{i},'fontweight','bold','fontsize',10)
    xlim([x_lim_min(i_var) x_lim_max(i_var)])% [0 12])
    set(gca,'fontsize',8,'fontweight','bold')

end


end


end



%% Chlorophyll

if(plot_chl)
    
    id = find(olevel<=50);
    depth_chl = olevel(id);
    
    %model
    tchl = nanmean(TCHL(id,end-365:end));
    % obs
    chl_obs = CHL_m_;
    
    fig = figure;
    plot(tchl,'linewidth',2)
    hold on
    plot(chl_obs,'.','markersize',10)
    ylabel('Chlorophyll-a')
    xlabel ('Time (days)')
    
    legend('Mod','Obs')
    
    %print('Figures/validate_biogeoch_CHL','-dpng')
end

add(R,Figure)

close(R)
