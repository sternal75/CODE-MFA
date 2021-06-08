% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% Cumulative distribution of reaction metabolite concentration confidence 
% interval sizes inferred by CODE-MFA versus with strictly 
% thermodynamic analysis 
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

known_net_flux_directions=nan(length(model_thermodynamics.rxns),1);
MILP_bounds_results = milp_find_bounds(model_net_fluxes, model_thermodynamics, known_net_flux_directions);


for(i=1:length(model_thermodynamics.mets))
    % low/high concentration (log10)
    low_concentration_code_mfa(i)   = log10(exp(sensitiviy_analysis_concentration{i}.low_concentration));
    high_concentration_code_mfa(i)  = log10(exp(sensitiviy_analysis_concentration{i}.high_concentration));   
    
    low_concentration_based_on_measurements(i)   = log10(exp(model_thermodynamics.mets_lb(i)));
    high_concentration_based_on_measurements(i)  = log10(exp(model_thermodynamics.mets_ub(i)));       
    
    low_concentration_based_on_thermodynamics(i)   = log10(exp(MILP_bounds_results.ln_C.min(i)));
    high_concentration_based_on_thermodynamics(i)  = log10(exp(MILP_bounds_results.ln_C.max(i)));
end 


% cumulative distribution - all metabolites
% MAX_confidence_interval_for_axes_1000 = 1000;
our_method_confidence_intervals_size = high_concentration_code_mfa-low_concentration_code_mfa;
% our_method_confidence_intervals_size(our_method_confidence_intervals_size>MAX_confidence_interval_for_axes_1000)=MAX_confidence_interval_for_axes_1000;
our_method_confidence_intervals_size = our_method_confidence_intervals_size(model_thermodynamics.skip_sensitivity_analysis_metabolite_indices~=1);
%mets = model_thermodynamics.mets(model_thermodynamics.skip_sensitivity_analysis_metabolite_indices~=1);

confidence_intervals_size_based_on_measurements = high_concentration_based_on_measurements-low_concentration_based_on_measurements;
confidence_intervals_size_based_on_measurements = confidence_intervals_size_based_on_measurements(model_thermodynamics.skip_sensitivity_analysis_metabolite_indices~=1);

confidence_intervals_size_based_on_thermodynamics = high_concentration_based_on_thermodynamics-low_concentration_based_on_thermodynamics;
confidence_intervals_size_based_on_thermodynamics = confidence_intervals_size_based_on_thermodynamics(model_thermodynamics.skip_sensitivity_analysis_metabolite_indices~=1);



x_our_method = unique(sort((our_method_confidence_intervals_size)));
x_based_on_measurements = unique(sort((confidence_intervals_size_based_on_measurements)));
x_based_on_thermodynamics = unique(sort((confidence_intervals_size_based_on_thermodynamics)));
our_method_num_of_sensitivity_analysis_range_smaller_than_x=[];
based_on_measurements_num_of_sensitivity_analysis_range_smaller_than_x=[];
based_on_thermodynamics_num_of_sensitivity_analysis_range_smaller_than_x=[];
for(i=1:length(x_our_method))
    our_method_num_of_sensitivity_analysis_range_smaller_than_x(i)  =  sum((our_method_confidence_intervals_size) <= x_our_method(i));
end
for(i=1:length(x_based_on_measurements))
    based_on_measurements_num_of_sensitivity_analysis_range_smaller_than_x(i)  =  sum((confidence_intervals_size_based_on_measurements) <= x_based_on_measurements(i));
end
for(i=1:length(x_based_on_thermodynamics))
    based_on_thermodynamics_num_of_sensitivity_analysis_range_smaller_than_x(i)  =  sum((confidence_intervals_size_based_on_thermodynamics) <= x_based_on_thermodynamics(i));
end

figure('Position',[1 1 1100 950]);
hold on;
% plot(x_our_method, our_method_num_of_sensitivity_analysis_range_smaller_than_x/length(our_method_confidence_intervals_size),'-x','MarkerSize',14,'LineWidth',2);
plot(x_our_method, our_method_num_of_sensitivity_analysis_range_smaller_than_x/length(our_method_confidence_intervals_size),'-*','MarkerSize',14,'LineWidth',2);
% plot(x_based_on_measurements, based_on_measurements_num_of_sensitivity_analysis_range_smaller_than_x/length(confidence_intervals_size_based_on_measurements),'-*','MarkerSize',14,'LineWidth',2);
plot(x_based_on_thermodynamics, based_on_thermodynamics_num_of_sensitivity_analysis_range_smaller_than_x/length(confidence_intervals_size_based_on_thermodynamics),'-*','MarkerSize',14,'LineWidth',2);
% plot(x_our_method, our_method_num_of_sensitivity_analysis_range_smaller_than_x/length(our_method_confidence_intervals_size),'.-','MarkerSize',26,'LineWidth',2);
hold off;
xlabel('Metabolite concentration confidence intervals [orders of magnitude]');
ylabel('Fraction of metabolites');
grid on
set(gcf,'color','w');
set(gca, 'FontSize', 22);
ax=gca;
ax.Box='on';
ylim([0 1.02]);
xlim([0 str2num(ax.XTickLabel{end})]);
% ax.XTickLabel{end}=['>' ax.XTickLabel{end}]
% legend('CoDe-MFA','Measurements', 'Thermodynamics and measurements', 'FontSize',18);
legend('CoDe-MFA','Thermodynamics and measurements', 'FontSize',18, 'Location','northwest');


% calculate Wilcoxon rank sum test
p = ranksum(confidence_intervals_size_based_on_thermodynamics,our_method_confidence_intervals_size)
fprintf('Wilcoxon rank sum test for method based on thermodynamics compared to CoDe-MFA. P-value=%e\n', p);

set(gcf,'PaperSize',[50 30]);
s = sprintf('./output_images/3h.pdf');
saveas(gcf, s);
