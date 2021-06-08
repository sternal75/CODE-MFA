% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% Calculate confidence intervals for metabolite concentrations
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
load('mat_files/model.mat','model');
load('mat_files/model_net_fluxes.mat','model_net_fluxes');
load('mat_files/model_thermodynamics.mat','model_thermodynamics');    
load('mat_files/EMU.mat','EMU');
load('mat_files/idv.mat','idv');
load('mat_files/directionalities.mat', 'directionalities');    
load('mat_files/WC_known_metabolites_idv.mat','WC_known_metabolites_idv');
load('mat_files/WC_known_metabolites_concentration.mat','WC_known_metabolites_concentration');
load('mat_files/MILP_results.mat','MILP_results');
load('mat_files/iterations_result.mat','iterations_result');    


addpath('./functions/emu') 
addpath('./functions/general') 

load_constants



% sort all directionalities per their error, in order to start
% sensitivity analysis with the order of the errors. This is supposed
% to be faster as the convergence is faster
[val ind]=sort(directionalities.errors);
directionalities.errors=directionalities.errors(ind);
directionalities.directionality_matrix=directionalities.directionality_matrix(:,ind);
directionalities.predicted_net_fluxes = directionalities.predicted_net_fluxes(:,ind);
directionalities.predicted_fb_fluxes = directionalities.predicted_fb_fluxes(:,ind);
directionalities.predicted_concentrations = directionalities.predicted_concentrations(:,ind);


SENSITIVITY_ANALYSIS_JUMPS_FOR_CONCENTRATION = [1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20];            


%     find the index of the best score among all directionalities
best_score = min(directionalities.errors);
index_best_score = find(directionalities.errors==min(directionalities.errors));

model_positive_lower_bound=model_net_fluxes.positive_direction_lb;
model_positive_upper_bound=model_net_fluxes.positive_direction_ub;
model_negative_lower_bound=model_net_fluxes.negative_direction_lb;
model_negative_upper_bound=model_net_fluxes.negative_direction_ub;    
model_lower_bound = -model_negative_upper_bound;
model_lower_bound(isnan(model_lower_bound))=model_positive_lower_bound(isnan(model_lower_bound));
model_upper_bound = model_positive_upper_bound;

% predicted values of best score over all directionalities
best_score_predicted_fluxes_fb  = directionalities.predicted_fb_fluxes(:,index_best_score);    
best_score_predicted_net_fluxes = directionalities.predicted_net_fluxes(:,index_best_score);
best_score_predicted_concentrations     = directionalities.predicted_concentrations(:,index_best_score);    
best_score_directionalities             = directionalities.directionality_matrix(:,index_best_score);    


sensitivity_analysis_high_fluxes            = cell(length(model_thermodynamics.mets),1);
sensitivity_analysis_high_fluxes_fb         = cell(length(model_thermodynamics.mets),1);
sensitivity_analysis_high_concentrations    = cell(length(model_thermodynamics.mets),1);
sensitivity_analysis_low_fluxes             = cell(length(model_thermodynamics.mets),1);
sensitivity_analysis_low_fluxes_fb          = cell(length(model_thermodynamics.mets),1);
sensitivity_analysis_low_concentrations     = cell(length(model_thermodynamics.mets),1);


all_indices_within_confidence_intervals = find(directionalities.errors < (min(directionalities.errors)+CONSTANT_VALUE_FOR_CONFIDENCE_INTERVAL));
concentration_sensitivity_analysis_min_max = cell(length(model_thermodynamics.mets),1);
res = cell(length(model_thermodynamics.mets),1);
for(i=1:length(all_indices_within_confidence_intervals))
    current_directionalities = directionalities.directionality_matrix(:,all_indices_within_confidence_intervals(i)); 
    MILP_bounds_results_all_indices_within_confidence_intervals{i} = milp_find_bounds(model_net_fluxes, model_thermodynamics, current_directionalities);
    for(j=1:length(model_thermodynamics.mets))
        concentration_sensitivity_analysis_min_max{j}(end+1,1)  = MILP_bounds_results_all_indices_within_confidence_intervals{i}.ln_C.min(j);
        concentration_sensitivity_analysis_min_max{j}(end,2)    = MILP_bounds_results_all_indices_within_confidence_intervals{i}.ln_C.max(j);    
        
        if(model_thermodynamics.skip_sensitivity_analysis_metabolite_indices(j)==1)
            res{j}.skip_sensitivity_analysis_concentration=1;
        else
            res{j}.skip_sensitivity_analysis_concentration=0;
        end    
        
        sensitivity_analysis_high_fluxes{j}         = best_score_predicted_net_fluxes;
        sensitivity_analysis_low_fluxes{j}          = best_score_predicted_net_fluxes;
        sensitivity_analysis_high_fluxes_fb{j}      = best_score_predicted_fluxes_fb;        
        sensitivity_analysis_low_fluxes_fb{j}       = best_score_predicted_fluxes_fb;    
        sensitivity_analysis_high_concentrations{j} = best_score_predicted_concentrations;
        sensitivity_analysis_low_concentrations{j}  = best_score_predicted_concentrations;                        
    end
end




% go to the right
last_sensitivity_analysis_right             = cell(length(model_thermodynamics.mets),1);
parfor(i=1:length(model_thermodynamics.mets))    
    sensitivity_analysis_index = 1;   
    used_fva_indexes=zeros(length(all_indices_within_confidence_intervals),1);    
    consecutive_counter = 0;
    sensitivity_analysis_concentration = best_score_predicted_concentrations(i)+SENSITIVITY_ANALYSIS_JUMPS_FOR_CONCENTRATION(1);
    while(true)      
        consecutive_counter = consecutive_counter+1;
        current_indexes = find((sensitivity_analysis_concentration>(concentration_sensitivity_analysis_min_max{i}(:,1)+1e-5))&(sensitivity_analysis_concentration<(concentration_sensitivity_analysis_min_max{i}(:,2)-1e-5))&(used_fva_indexes~=1));
        if(isempty(current_indexes))
            break;
        end
        

        if(sum(used_fva_indexes)==length(all_indices_within_confidence_intervals))
            if(consecutive_counter==length(all_indices_within_confidence_intervals))
                break;  %finish
            else
                used_fva_indexes=zeros(length(all_indices_within_confidence_intervals),1);    
                consecutive_counter = 0;
            end
        end
        used_fva_indexes(current_indexes(1))=1;        

        while (sensitivity_analysis_index <= length(SENSITIVITY_ANALYSIS_JUMPS_FOR_CONCENTRATION))
            fprintf('+++ met con index=%d, sa index=%d, best score con=%2f, sa high con=%2f\n', i, sensitivity_analysis_index, best_score_predicted_concentrations(i), sensitivity_analysis_high_concentrations{i}(i));   

            current_predicted_fluxes_fb         = directionalities.predicted_fb_fluxes(:,all_indices_within_confidence_intervals(current_indexes(1)));
            current_predicted_fluxes_net        = directionalities.predicted_net_fluxes(:,all_indices_within_confidence_intervals(current_indexes(1)));
            current_predicted_concentrations    = directionalities.predicted_concentrations(:,all_indices_within_confidence_intervals(current_indexes(1)));            

            current_directionalities = directionalities.directionality_matrix(:,all_indices_within_confidence_intervals(current_indexes(1)));             

            initial_fluxes_net          = current_predicted_fluxes_net;
            initial_fluxes_fb           = current_predicted_fluxes_fb;
            initial_concentrations      = current_predicted_concentrations;
            
            sensitivity_analysis_concentration = best_score_predicted_concentrations(i)+SENSITIVITY_ANALYSIS_JUMPS_FOR_CONCENTRATION(sensitivity_analysis_index);
            
            max_concentration = concentration_sensitivity_analysis_min_max{i}(current_indexes(1),2);
            min_concentration = concentration_sensitivity_analysis_min_max{i}(current_indexes(1),1);
            if ((sensitivity_analysis_concentration > max_concentration-EPSILON_VALUE) || (sensitivity_analysis_concentration < min_concentration+EPSILON_VALUE))
                break;
            end

            
            initial_concentrations(i) = sensitivity_analysis_concentration;

            directionality_vector = current_directionalities;
            MILP_bounds_results = MILP_bounds_results_all_indices_within_confidence_intervals{current_indexes(1)};
            current_MILP_bounds_results = MILP_bounds_results;
            current_MILP_bounds_results.ln_C.min(i) = sensitivity_analysis_concentration;
            current_MILP_bounds_results.ln_C.max(i) = sensitivity_analysis_concentration;      

            mymodel = model;

            tic;
            [exitflag error error_match_labeling error_match_concentrations predicted_net_flux predicted_fb_flux concentrations number_of_independent_variables number_of_fitted_elements] = ComputeEMUOptFlux(current_directionalities, mymodel, model_net_fluxes, model_thermodynamics, EMU, idv, WC_known_metabolites_idv, WC_known_metabolites_concentration, initial_fluxes_fb, initial_fluxes_net, initial_concentrations, current_MILP_bounds_results,best_score+CONSTANT_VALUE_FOR_CONFIDENCE_INTERVAL);
            elapsed_time=toc;

            if(error > best_score+CONSTANT_VALUE_FOR_CONFIDENCE_INTERVAL)
                break;
            end

            if ((exitflag ~= 1)&&(exitflag ~= 2)&&(exitflag ~= -3)&&(exitflag ~= 0))
                fprintf('Error in fmincon\n');
                break;
            end

            sensitivity_analysis_high_fluxes{i}            = predicted_net_flux;
            sensitivity_analysis_high_fluxes_fb{i}         = predicted_fb_flux;
            sensitivity_analysis_high_concentrations{i}    = concentrations;

            sensitivity_analysis_index = sensitivity_analysis_index+1;
            consecutive_counter = 0;
        end    
    end
    last_sensitivity_analysis_right{i} = sensitivity_analysis_concentration;    
    if(last_sensitivity_analysis_right{i} > max(concentration_sensitivity_analysis_min_max{i}(:,2))-EPSILON_VALUE)
        last_sensitivity_analysis_right{i} = max(concentration_sensitivity_analysis_min_max{i}(:,2));
    end
        
end


% go to the left
last_sensitivity_analysis_left              = cell(length(model_thermodynamics.mets),1);
parfor(i=1:length(model_thermodynamics.mets))
    
    sensitivity_analysis_index = 1;    
    used_fva_indexes=zeros(length(all_indices_within_confidence_intervals),1);    
    consecutive_counter = 0;
    sensitivity_analysis_concentration = best_score_predicted_concentrations(i)-SENSITIVITY_ANALYSIS_JUMPS_FOR_CONCENTRATION(1);
    while(true)                

        consecutive_counter = consecutive_counter+1;
        current_indexes = find((sensitivity_analysis_concentration>(concentration_sensitivity_analysis_min_max{i}(:,1)+1e-5))&(sensitivity_analysis_concentration<(concentration_sensitivity_analysis_min_max{i}(:,2)-1e-5))&(used_fva_indexes~=1));
        if(isempty(current_indexes))
            break;
        end

        if(sum(used_fva_indexes)==length(all_indices_within_confidence_intervals))
            if(consecutive_counter==length(all_indices_within_confidence_intervals))
                break;  %finish
            else
                used_fva_indexes=zeros(length(all_indices_within_confidence_intervals),1);    
                consecutive_counter = 0;
            end
        end        
        used_fva_indexes(current_indexes(1))=1;

        while (sensitivity_analysis_index <= length(SENSITIVITY_ANALYSIS_JUMPS_FOR_CONCENTRATION))

            fprintf('--- met con index=%d, sa index=%d, sa low con=%2f, best score con=%2f, sa high con=%2f\n', i, sensitivity_analysis_index, sensitivity_analysis_low_concentrations{i}(i), best_score_predicted_concentrations(i), sensitivity_analysis_high_concentrations{i}(i));   

            current_predicted_fluxes_fb = directionalities.predicted_fb_fluxes(:,all_indices_within_confidence_intervals(current_indexes(1)));
            current_predicted_fluxes_net = directionalities.predicted_net_fluxes(:,all_indices_within_confidence_intervals(current_indexes(1)));
            current_predicted_concentrations    = directionalities.predicted_concentrations(:,all_indices_within_confidence_intervals(current_indexes(1)));            

            current_directionalities = directionalities.directionality_matrix(:,all_indices_within_confidence_intervals(current_indexes(1))); 
            initial_fluxes_net          = current_predicted_fluxes_net;
            initial_fluxes_fb           = current_predicted_fluxes_fb;
            initial_concentrations      = current_predicted_concentrations;
            
            sensitivity_analysis_concentration = best_score_predicted_concentrations(i)-SENSITIVITY_ANALYSIS_JUMPS_FOR_CONCENTRATION(sensitivity_analysis_index);

            
            max_concentration = concentration_sensitivity_analysis_min_max{i}(current_indexes(1),2);
            min_concentration = concentration_sensitivity_analysis_min_max{i}(current_indexes(1),1);
            if ((sensitivity_analysis_concentration > max_concentration-EPSILON_VALUE) || (sensitivity_analysis_concentration < min_concentration+EPSILON_VALUE))
                break;
            end
                        

            initial_concentrations(i) = sensitivity_analysis_concentration;

            directionality_vector = current_directionalities;
            MILP_bounds_results = MILP_bounds_results_all_indices_within_confidence_intervals{current_indexes(1)};
            current_MILP_bounds_results = MILP_bounds_results;
            current_MILP_bounds_results.ln_C.min(i) = sensitivity_analysis_concentration;
            current_MILP_bounds_results.ln_C.max(i) = sensitivity_analysis_concentration;      
            
            mymodel = model;

            tic;
            [exitflag error error_match_labeling error_match_concentrations predicted_net_flux predicted_fb_flux concentrations number_of_independent_variables number_of_fitted_elements] = ComputeEMUOptFlux(current_directionalities, mymodel, model_net_fluxes, model_thermodynamics, EMU, idv, WC_known_metabolites_idv, WC_known_metabolites_concentration, initial_fluxes_fb, initial_fluxes_net, initial_concentrations, current_MILP_bounds_results,best_score+CONSTANT_VALUE_FOR_CONFIDENCE_INTERVAL);
            elapsed_time=toc;

            if(error > best_score+CONSTANT_VALUE_FOR_CONFIDENCE_INTERVAL)
                break;
            end

            if ((exitflag ~= 1)&&(exitflag ~= 2)&&(exitflag ~= -3)&&(exitflag ~= 0))
                fprintf('Error in fmincon\n');
                break;
            end

            sensitivity_analysis_low_fluxes{i}         = predicted_net_flux;
            sensitivity_analysis_low_fluxes_fb{i}      = predicted_fb_flux;      
            sensitivity_analysis_low_concentrations{i} = concentrations;
            
            sensitivity_analysis_index = sensitivity_analysis_index+1;
            consecutive_counter = 0;
        end    
    end
    
    last_sensitivity_analysis_left{i} = sensitivity_analysis_concentration;
    if(last_sensitivity_analysis_left{i} < min(concentration_sensitivity_analysis_min_max{i}(:,1))+EPSILON_VALUE)
        last_sensitivity_analysis_left{i} = min(concentration_sensitivity_analysis_min_max{i}(:,1));
    end            
end


for(i=1:length(model_thermodynamics.mets))
    res{i}.metabolite_index    = i;
    res{i}.metabolite_name     = model_thermodynamics.mets{i};
    res{i}.low_concentration        = last_sensitivity_analysis_left{i};
    res{i}.high_concentration       = last_sensitivity_analysis_right{i};
    res{i}.net_fluxes.low  = sensitivity_analysis_low_fluxes{i};
    res{i}.net_fluxes.high = sensitivity_analysis_high_fluxes{i};
    res{i}.fb_fluxes.low   = sensitivity_analysis_low_fluxes_fb{i};
    res{i}.fb_fluxes.high  = sensitivity_analysis_high_fluxes_fb{i};
    res{i}.concentrations.low   = sensitivity_analysis_low_concentrations{i};
    res{i}.concentrations.high  = sensitivity_analysis_high_concentrations{i};
end


sensitiviy_analysis_concentration=res;
save('mat_files/sensitiviy_analysis_concentration.mat','sensitiviy_analysis_concentration');



