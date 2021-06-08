% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% Cumulative distribution of reaction Gibbs energy confidence interval 
% sizes inferred by CODE-MFA versus with strictly thermodynamic analysis
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

known_net_flux_directions=nan(length(model_thermodynamics.rxns),1);
MILP_bounds_results = milp_find_bounds(model_net_fluxes, model_thermodynamics, known_net_flux_directions);

% find min/max Gibbs free energy based on measurements
for(i=1:model_net_fluxes.rxn_num)
    f=zeros(1,length(MILP_bounds_results.dG));
    f(model_thermodynamics.product_indexes{i})  = RT;
    f(model_thermodynamics.reactant_indexes{i}) = -RT;
    
    [optimization_values,fval,exitflag,output] = linprog(f, [], [], [], [], model_thermodynamics.mets_lb, model_thermodynamics.mets_ub);
    if(isempty(optimization_values))
        fprintf('intlinprog is empty\n');
    end
    if ((exitflag ~= 1)&&(exitflag ~= 2)&&(exitflag ~= -3))
        fprintf('Error in intlinprog (dG) - %d\n',i);      
    else
        dG_based_on_measurements.min(i) = fval+model_thermodynamics.delta_G0(i);
    end                    
    
    
    f=-f;
    [optimization_values,fval,exitflag,output] = linprog(f, [], [], [], [], model_thermodynamics.mets_lb, model_thermodynamics.mets_ub);
    if(isempty(optimization_values))
        fprintf('intlinprog is empty\n');
    end
    if ((exitflag ~= 1)&&(exitflag ~= 2)&&(exitflag ~= -3))
        fprintf('Error in intlinprog (dG) - %d\n',i);
    else
        dG_based_on_measurements.max(i) = -fval+model_thermodynamics.delta_G0(i);    
    end
    
end





for(i=1:length(sensitiviy_analysis_dG))
    low_dG_code_mfa(i) = sensitiviy_analysis_dG{i}.low_dG;
    high_dG_code_mfa(i) = sensitiviy_analysis_dG{i}.high_dG;
        
    low_dG_based_on_measurements(i)       = dG_based_on_measurements.min(i);
    high_dG_based_on_measurements(i)      = dG_based_on_measurements.max(i);
    
    low_dG_based_on_thermodynamics(i)       = MILP_bounds_results.dG.min(i);
    high_dG_based_on_thermodynamics(i)      = MILP_bounds_results.dG.max(i);
    
end 


confidence_intervals_size_code_mfa = abs(high_dG_code_mfa-low_dG_code_mfa);
confidence_intervals_size_code_mfa = confidence_intervals_size_code_mfa(model_thermodynamics.thermodynamics_of_reaction_defined==1);
confidence_intervals_size_based_on_measurements = abs(high_dG_based_on_measurements-low_dG_based_on_measurements);
confidence_intervals_size_based_on_measurements = confidence_intervals_size_based_on_measurements(model_thermodynamics.thermodynamics_of_reaction_defined==1);
confidence_intervals_size_based_on_thermodynamics = abs(high_dG_based_on_thermodynamics-low_dG_based_on_thermodynamics);
confidence_intervals_size_based_on_thermodynamics = confidence_intervals_size_based_on_thermodynamics(model_thermodynamics.thermodynamics_of_reaction_defined==1);



x_code_mfa = unique(sort(log10(confidence_intervals_size_code_mfa)));
x_based_on_measurements = unique(sort(log10(confidence_intervals_size_based_on_measurements)));
x_based_on_thermodynamics = unique(sort(log10(confidence_intervals_size_based_on_thermodynamics)));

code_mfa_num_of_sensitivity_analysis_range_smaller_than_x=[];
based_on_measurements_num_of_sensitivity_analysis_range_smaller_than_x=[];
based_on_thermodynamics_num_of_sensitivity_analysis_range_smaller_than_x=[];
for(i=1:length(x_code_mfa))
    code_mfa_num_of_sensitivity_analysis_range_smaller_than_x(i)  =  sum(log10(confidence_intervals_size_code_mfa) <= x_code_mfa(i));
end
for(i=1:length(x_based_on_measurements))
    based_on_measurements_num_of_sensitivity_analysis_range_smaller_than_x(i)  =  sum(log10(confidence_intervals_size_based_on_measurements) <= x_based_on_measurements(i));
end
for(i=1:length(x_based_on_thermodynamics))
    based_on_thermodynamics_num_of_sensitivity_analysis_range_smaller_than_x(i)  =  sum(log10(confidence_intervals_size_based_on_thermodynamics) <= x_based_on_thermodynamics(i));
end
figure('Position',[1 1 1100 950]);
hold on;
plot(x_code_mfa, code_mfa_num_of_sensitivity_analysis_range_smaller_than_x/length(confidence_intervals_size_code_mfa),'-*','MarkerSize',14,'LineWidth',2);
%plot(x_based_on_measurements, based_on_measurements_num_of_sensitivity_analysis_range_smaller_than_x/length(confidence_intervals_size_based_on_measurements),'-*','MarkerSize',14,'LineWidth',2);
plot(x_based_on_thermodynamics, based_on_thermodynamics_num_of_sensitivity_analysis_range_smaller_than_x/length(confidence_intervals_size_based_on_thermodynamics),'-*','MarkerSize',14,'LineWidth',2);
hold off;
xlabel('dG confidence interval size [log10(kJ/mol)]');
ylabel('Fraction of reactions');
legend('CoDe-MFA','MFA' );
grid on
set(gcf,'color','w');
set(gca, 'FontSize', 28);
ax=gca;
ax.Box='on';
ylim([0 1.02]);
xlim([0 str2num(ax.XTickLabel{end})]);
%ax.XTickLabel{end}=['>' ax.XTickLabel{end}]
%legend('CoDe-MFA','Measurements', 'Thermodynamics and measurements', 'FontSize',18);
legend('CoDe-MFA','Thermodynamics and measurements', 'FontSize',18, 'Location','northwest');


% calculate Wilcoxon rank sum test
p = ranksum(confidence_intervals_size_based_on_thermodynamics,confidence_intervals_size_code_mfa)
fprintf('Wilcoxon rank sum test for method based on thermodynamics compared to CoDe-MFA. P-value=%e\n', p);

set(gcf,'PaperSize',[50 30]);
s = sprintf('./output_images/3i.pdf');
saveas(gcf, s);
