% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% The fit been simulated total cellular metabolite concentration 
% (i.e. convolution of simulated mitochondrial and cytosolic labeling; 
% x-axis) and measurements (y-axis) 
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
load('../mat_files/sensitiviy_analysis_concentration.mat', 'sensitiviy_analysis_concentration');


addpath('../functions/emu') 
addpath('../functions/general') 
addpath('../') 
run ../load_constants;

%     find the index of the best score among all directionalities
best_score = min(directionalities.errors);
index_best_score = find(directionalities.errors==min(directionalities.errors));
best_score_predicted_concentrations     = directionalities.predicted_concentrations(:,index_best_score);    

measured_concentration_lb = [];
measured_concentration_ub = [];
WC_convoluted_concentration = [];
output_for_excel = cell(0);

% metabolites exist in one compartment
metabolite_exists_in_one_compartment = {'6phosphogluconate_CY','Glc6P_CY','Ribose5P_CY','G3P','13BPG','3PG'};
for(i=1:length(metabolite_exists_in_one_compartment))
    index_one_compartment_concentration = find(ismember(model_thermodynamics.mets,metabolite_exists_in_one_compartment{i}));
    measure_lb = exp(model_thermodynamics.mets_lb(index_one_compartment_concentration));
    measure_ub = exp(model_thermodynamics.mets_ub(index_one_compartment_concentration));
    % if measured metabolite has too high confidence intervals
    if((measure_ub/measure_lb) > 10000) && ((measure_ub-measure_lb)>0.1)
        continue;
    end
    measured_concentration_lb(end+1) = measure_lb;
    measured_concentration_ub(end+1) = measure_ub;
    
    WC_convoluted_concentration(end+1)  = exp(best_score_predicted_concentrations(index_one_compartment_concentration));
    
    output_for_excel{end+1} = metabolite_exists_in_one_compartment{i};
end


% metabolites exist in both compartment
for i=1:length(WC_known_metabolites_concentration)        
    x_cy_for_concentrations = WC_known_metabolites_concentration{i}.index_CY;
    x_mt_for_concentrations = WC_known_metabolites_concentration{i}.index_MT;

    measured_concentration_lb(end+1) = WC_known_metabolites_concentration{i}.concentration-2*WC_known_metabolites_concentration{i}.concentrations_STD;
    if(measured_concentration_lb(end)<1e-5)
        measured_concentration_lb(end) = 0.8*WC_known_metabolites_concentration{i}.concentration;
    end
    measured_concentration_ub(end+1) = WC_known_metabolites_concentration{i}.concentration+2*WC_known_metabolites_concentration{i}.concentrations_STD;
    simulated_cy_concentration  = exp(best_score_predicted_concentrations(x_cy_for_concentrations));
    simulated_mt_concentration  = exp(best_score_predicted_concentrations(x_mt_for_concentrations));
    WC_convoluted_concentration(end+1) = CY_WC_VOLUME*simulated_cy_concentration+(1-CY_WC_VOLUME)*simulated_mt_concentration;
    output_for_excel{end+1} = WC_known_metabolites_concentration{i}.met_name;
end

figure('Position',[1 1 1100 950]);
hold on;
plot(log10([WC_convoluted_concentration;WC_convoluted_concentration]), log10([measured_concentration_lb;measured_concentration_ub]),'Color',BLUE_COLOR,'LineWidth',3);
hold off;
xlabel('Simulated concentrations [log10(mM)]');
ylabel('Measured concentrations [log10(mM)]');
grid on
set(gcf,'color','w');
set(gca, 'FontSize', 28);
ax=gca;
ax.Box='on';
line([-4 2], [-4 2], 'color',RED_COLOR,'LineStyle',':','linewidth',2);
% ylim([-3 1.5]);
% xlim([-3 1.5]);
    

    
    
% print gln labeling to excel
xlswrite('temp.xlsx',{'met name'},'temp','A1');  
xlswrite('temp.xlsx',{'measured lb'},'temp','b1');  
xlswrite('temp.xlsx',{'measured ub'},'temp','c1');  
xlswrite('temp.xlsx',{'CODE MFA - best fit concentration'},'temp','d1');  
xlswrite('temp.xlsx',output_for_excel','temp','A2');  
xlswrite('temp.xlsx',measured_concentration_lb','temp','b2');  
xlswrite('temp.xlsx',measured_concentration_ub','temp','c2');  
xlswrite('temp.xlsx',WC_convoluted_concentration','temp','d2');  

set(gcf,'PaperSize',[50 30]);
s = sprintf('./output_images/3b.pdf');
saveas(gcf, s);
    
    
