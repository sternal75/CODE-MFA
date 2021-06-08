% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% Number of reactions whose direction of net flux is uniquely inferred 
% across CODE-MFA iterations
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 

clear all;
close all;


load('../mat_files/net_fluxes.mat', 'net_fluxes');
load('../mat_files/directionalities.mat', 'directionalities');
load('../mat_files/model_thermodynamics.mat','model_thermodynamics');    
load('../mat_files/model_net_fluxes.mat','model_net_fluxes');
load('../mat_files/model.mat','model');
load('../mat_files/EMU.mat','EMU');
load('../mat_files/idv.mat','idv');
load('../mat_files/WC_known_metabolites_idv.mat','WC_known_metabolites_idv');
load('../mat_files/WC_known_metabolites_concentration.mat','WC_known_metabolites_concentration');
load('../mat_files/MILP_results.mat','MILP_results');
load('../mat_files/iterations_result.mat','iterations_result');
load('../mat_files/sensitiviy_analysis_dG.mat', 'sensitiviy_analysis_dG');


addpath('../functions/emu') 
addpath('../functions/general') 
addpath('../') 
run ../load_constants;

[d_known_net_flux_directions,t_known_net_flux_directions] = xlsread('../xls_input_files/known_directionalities_per_iteration.xlsx');

known_net_flux_directiony_per_iteration = sum(~isnan(d_known_net_flux_directions));
known_net_flux_directiony_per_iteration(2:end+1) = known_net_flux_directiony_per_iteration;
known_net_flux_directiony_per_iteration(1) = sum(~model_net_fluxes.is_net_flux);

% non biological reactions are the ones that were added to the model just
% for mathematical reasons
NON_BIOLOGICAL_REACTIONS = 17;
known_net_flux_directiony_per_iteration = known_net_flux_directiony_per_iteration-NON_BIOLOGICAL_REACTIONS;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% all reactions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure('Position',[1 1 1100 950]);
hold on;
plot([0:length(known_net_flux_directiony_per_iteration)-1], (known_net_flux_directiony_per_iteration),'-*','MarkerSize',14,'LineWidth',2);
hold off;
xlabel('Number of algorithm iterations');
ylabel('Number of reactions with unique direction');
grid on
set(gcf,'color','w');
set(gca, 'FontSize', 29);
ax=gca;
ax.Box='on';

total_number_of_reaction = size(d_known_net_flux_directions,1)-NON_BIOLOGICAL_REACTIONS;
ylim([0 total_number_of_reaction+6]);
line([0 size(d_known_net_flux_directions,2)], [total_number_of_reaction total_number_of_reaction], 'color',RED_COLOR,'LineStyle',':','linewidth',2);
text(1,total_number_of_reaction+2.5,'Reactions in the model','color',RED_COLOR, 'FontSize', 32);

set(gcf,'PaperSize',[50 30]);
s = sprintf('./output_images/3a.pdf');
saveas(gcf, s);




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% only reactions with unknown net flux directionality
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
known_net_flux_directiony_bidirectional_fluxes_per_iteration = known_net_flux_directiony_per_iteration-known_net_flux_directiony_per_iteration(1);
figure('Position',[1 1 1100 950]);
hold on;
plot([0:length(known_net_flux_directiony_bidirectional_fluxes_per_iteration)-1], (known_net_flux_directiony_bidirectional_fluxes_per_iteration),'-*','MarkerSize',14,'LineWidth',2);
hold off;
xlabel('Number of algorithm iterations');
ylabel('Number of reactions with unique direction');
grid on
set(gcf,'color','w');
set(gca, 'FontSize', 29);
ax=gca;
ax.Box='on';
total_number_of_reaction = size(d_known_net_flux_directions,1)-NON_BIOLOGICAL_REACTIONS-known_net_flux_directiony_per_iteration(1);
ylim([0 total_number_of_reaction+3]);
%xlim([0 str2num(ax.XTickLabel{end})]);
%ax.XTickLabel{end}=['>' ax.XTickLabel{end}]
line([0 size(d_known_net_flux_directions,2)], [total_number_of_reaction total_number_of_reaction], 'color',RED_COLOR,'LineStyle',':','linewidth',2);
text(0.1,total_number_of_reaction+1,'Reactions with unknown directionality in the model','color',RED_COLOR, 'FontSize', 28);
