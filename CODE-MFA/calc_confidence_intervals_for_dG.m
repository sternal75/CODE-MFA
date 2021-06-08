% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% Calculate confidence intervals for gibbs free energies
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


SENSITIVITY_ANALYSIS_JUMPS_FOR_dG = [2 4 6 8 10 12 14 16 18 20 25 30 35 40 45 50 55 60 65 70 80 90 100 120 140 160 200];            %fixed size jumps in case the flux is too small to avoid too sensitivity analysis checked values


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
for(i=1:length(model_net_fluxes.rxns))
    numerator=0;
    denominator=0;
    numerator(model_thermodynamics.product_indexes{i})      = (best_score_predicted_concentrations(model_thermodynamics.product_indexes{i}));
    numerator=sum(numerator);
    denominator(model_thermodynamics.reactant_indexes{i})    = (best_score_predicted_concentrations(model_thermodynamics.reactant_indexes{i}));            
    denominator=sum(denominator);
    best_score_predicted_dG(i)         = model_thermodynamics.delta_G0(i)+RT*(numerator-denominator);        
end    


sensitivity_analysis_high_fluxes            = cell(length(model_net_fluxes.rxns),1);
sensitivity_analysis_high_fluxes_fb         = cell(length(model_net_fluxes.rxns),1);
sensitivity_analysis_high_concentrations    = cell(length(model_net_fluxes.rxns),1);
sensitivity_analysis_low_fluxes             = cell(length(model_net_fluxes.rxns),1);
sensitivity_analysis_low_fluxes_fb          = cell(length(model_net_fluxes.rxns),1);
sensitivity_analysis_low_concentrations     = cell(length(model_net_fluxes.rxns),1);
sensitivity_analysis_high_dG                = cell(length(model_net_fluxes.rxns),1);
sensitivity_analysis_low_dG                 = cell(length(model_net_fluxes.rxns),1);


all_indices_within_confidence_intervals = find(directionalities.errors < (min(directionalities.errors)+CONSTANT_VALUE_FOR_CONFIDENCE_INTERVAL));
dG_sensitivity_analysis_min_max = cell(length(model_net_fluxes.rxns),1);
res = cell(length(model_net_fluxes.rxns),1);
for(i=1:length(all_indices_within_confidence_intervals))
    current_directionalities = directionalities.directionality_matrix(:,all_indices_within_confidence_intervals(i)); 
    MILP_bounds_results_all_indices_within_confidence_intervals{i} = milp_find_bounds(model_net_fluxes, model_thermodynamics, current_directionalities);
    for(j=1:length(model_net_fluxes.rxns))
        dG_sensitivity_analysis_min_max{j}(end+1,1) = MILP_bounds_results_all_indices_within_confidence_intervals{i}.dG.min(j);
        dG_sensitivity_analysis_min_max{j}(end,2) = MILP_bounds_results_all_indices_within_confidence_intervals{i}.dG.max(j);    
        if(model_net_fluxes.skip_sensitivity_analysis_reaction_indices(j)==1)
            res{j}.skip_sensitivity_analysis_reaction=1;
        else
            res{j}.skip_sensitivity_analysis_reaction=0;
        end            

        sensitivity_analysis_high_fluxes{j}         = best_score_predicted_net_fluxes;
        sensitivity_analysis_low_fluxes{j}          = best_score_predicted_net_fluxes;
        sensitivity_analysis_high_fluxes_fb{j}      = best_score_predicted_fluxes_fb;        
        sensitivity_analysis_low_fluxes_fb{j}       = best_score_predicted_fluxes_fb;    
        sensitivity_analysis_high_concentrations{j} = best_score_predicted_concentrations;
        sensitivity_analysis_low_concentrations{j}  = best_score_predicted_concentrations; 
        sensitivity_analysis_high_dG{j}(:,1)             = best_score_predicted_dG; 
        sensitivity_analysis_low_dG{j}(:,1)              = best_score_predicted_dG;
    end        
end




% go to the right
last_sensitivity_analysis_right             = cell(length(model_net_fluxes.rxns),1);
parfor(i=1:length(model_net_fluxes.rxns))    
    sensitivity_analysis_index = 1;   
    used_fva_indexes=zeros(length(all_indices_within_confidence_intervals),1);    
    consecutive_counter = 0;
    sensitivity_analysis_dG = best_score_predicted_dG(i)+SENSITIVITY_ANALYSIS_JUMPS_FOR_dG(1);
    while(true)      
        consecutive_counter = consecutive_counter+1;
        current_indexes = find((sensitivity_analysis_dG>(dG_sensitivity_analysis_min_max{i}(:,1)+EPSILON_VALUE))&(sensitivity_analysis_dG<(dG_sensitivity_analysis_min_max{i}(:,2)-EPSILON_VALUE))&(used_fva_indexes~=1));
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

        while (sensitivity_analysis_index <= length(SENSITIVITY_ANALYSIS_JUMPS_FOR_dG))
            fprintf('+++ dG index=%d, sa index=%d, best score dG=%2f, sa high dG=%2f\n', i, sensitivity_analysis_index, best_score_predicted_dG(i), sensitivity_analysis_dG);   

            current_predicted_fluxes_fb = directionalities.predicted_fb_fluxes(:,all_indices_within_confidence_intervals(current_indexes(1)));
            current_predicted_fluxes_net = directionalities.predicted_net_fluxes(:,all_indices_within_confidence_intervals(current_indexes(1)));

            current_directionalities = directionalities.directionality_matrix(:,all_indices_within_confidence_intervals(current_indexes(1)));             

            initial_fluxes_net = current_predicted_fluxes_net;
            initial_fluxes_fb = current_predicted_fluxes_fb;
            sensitivity_analysis_dG = best_score_predicted_dG(i)+SENSITIVITY_ANALYSIS_JUMPS_FOR_dG(sensitivity_analysis_index);
            % when trying delta G with value between the zone that is too
            % close to zero, choose the closest permited one (closest to
            % zero)
            if(sensitivity_analysis_dG < dG_MIN_DISTANCE_FROM_ZERO)&&(sensitivity_analysis_dG >- dG_MIN_DISTANCE_FROM_ZERO)
                if(sensitivity_analysis_dG>=0)
                    sensitivity_analysis_dG = dG_MIN_DISTANCE_FROM_ZERO+2*EPSILON_VALUE; 
                else
                    sensitivity_analysis_dG = -dG_MIN_DISTANCE_FROM_ZERO-2*EPSILON_VALUE;
                end
            end

            max_dG = dG_sensitivity_analysis_min_max{i}(current_indexes(1),2);
            min_dG = dG_sensitivity_analysis_min_max{i}(current_indexes(1),1);
            if ((sensitivity_analysis_dG > max_dG-EPSILON_VALUE) || (sensitivity_analysis_dG < min_dG+EPSILON_VALUE))
                break;
            end

            initial_concentrations = directionalities.predicted_concentrations(:,all_indices_within_confidence_intervals(current_indexes(1)));

            directionality_vector = current_directionalities;
            MILP_bounds_results = MILP_bounds_results_all_indices_within_confidence_intervals{current_indexes(1)};
            current_MILP_bounds_results = MILP_bounds_results;
            current_MILP_bounds_results.dG.min(i) = sensitivity_analysis_dG;
            current_MILP_bounds_results.dG.max(i) = sensitivity_analysis_dG;   

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
            sensitivity_analysis_high_dG{i} = [];
            for(j=1:length(model_net_fluxes.rxns))
                numerator=0;
                denominator=0;
                numerator(model_thermodynamics.product_indexes{j})      = (sensitivity_analysis_high_concentrations{i}(model_thermodynamics.product_indexes{j}));
                numerator=sum(numerator);
                denominator(model_thermodynamics.reactant_indexes{j})    = (sensitivity_analysis_high_concentrations{i}(model_thermodynamics.reactant_indexes{j}));            
                denominator=sum(denominator);
                sensitivity_analysis_high_dG{i}(end+1,1)         = model_thermodynamics.delta_G0(j)+RT*(numerator-denominator);    
    
            end    

            sensitivity_analysis_index = sensitivity_analysis_index+1;
            consecutive_counter = 0;
        end    
    end

    last_sensitivity_analysis_right{i} = sensitivity_analysis_dG;
    if(last_sensitivity_analysis_right{i} > max(dG_sensitivity_analysis_min_max{i}(:,2))-EPSILON_VALUE)
        last_sensitivity_analysis_right{i} = max(dG_sensitivity_analysis_min_max{i}(:,2));
    end    
end

% go to the left
last_sensitivity_analysis_reft             = cell(length(model_net_fluxes.rxns),1);
parfor(i=1:length(model_net_fluxes.rxns))    
    sensitivity_analysis_index = 1;   
    used_fva_indexes=zeros(length(all_indices_within_confidence_intervals),1);    
    consecutive_counter = 0;
    sensitivity_analysis_dG = best_score_predicted_dG(i)-SENSITIVITY_ANALYSIS_JUMPS_FOR_dG(1);
    while(true)                
        consecutive_counter = consecutive_counter+1;
        current_indexes = find((sensitivity_analysis_dG>(dG_sensitivity_analysis_min_max{i}(:,1)+EPSILON_VALUE))&(sensitivity_analysis_dG<(dG_sensitivity_analysis_min_max{i}(:,2)-EPSILON_VALUE))&(used_fva_indexes~=1));
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

        while (sensitivity_analysis_index <= length(SENSITIVITY_ANALYSIS_JUMPS_FOR_dG))
            fprintf('--- dG index=%d, sa index=%d, sa low dG=%2f, best score dG=%2f, sa high dG=%2f\n', i, sensitivity_analysis_index, sensitivity_analysis_dG, best_score_predicted_dG(i), last_sensitivity_analysis_right{i});   

            current_predicted_fluxes_fb = directionalities.predicted_fb_fluxes(:,all_indices_within_confidence_intervals(current_indexes(1)));
            current_predicted_fluxes_net = directionalities.predicted_net_fluxes(:,all_indices_within_confidence_intervals(current_indexes(1)));

            current_directionalities = directionalities.directionality_matrix(:,all_indices_within_confidence_intervals(current_indexes(1))); 
            initial_fluxes_net = current_predicted_fluxes_net;
            initial_fluxes_fb = current_predicted_fluxes_fb;
            sensitivity_analysis_dG = best_score_predicted_dG(i)-SENSITIVITY_ANALYSIS_JUMPS_FOR_dG(sensitivity_analysis_index);
            % when trying delta G with value between the zone that is too
            % close to zero, choose the closest permited one (closest to
            % zero)
            if(sensitivity_analysis_dG < dG_MIN_DISTANCE_FROM_ZERO)&&(sensitivity_analysis_dG >- dG_MIN_DISTANCE_FROM_ZERO)
                if(sensitivity_analysis_dG>=0)
                    sensitivity_analysis_dG = dG_MIN_DISTANCE_FROM_ZERO+2*EPSILON_VALUE; 
                else
                    sensitivity_analysis_dG = -dG_MIN_DISTANCE_FROM_ZERO-2*EPSILON_VALUE;
                end
            end
            

            max_dG = dG_sensitivity_analysis_min_max{i}(current_indexes(1),2);
            min_dG = dG_sensitivity_analysis_min_max{i}(current_indexes(1),1);
            if ((sensitivity_analysis_dG > max_dG-EPSILON_VALUE) || (sensitivity_analysis_dG < min_dG+EPSILON_VALUE))
                break;
            end

            initial_concentrations = directionalities.predicted_concentrations(:,all_indices_within_confidence_intervals(current_indexes(1)));            

            directionality_vector = current_directionalities;
            MILP_bounds_results = MILP_bounds_results_all_indices_within_confidence_intervals{current_indexes(1)};
            current_MILP_bounds_results = MILP_bounds_results;
            current_MILP_bounds_results.dG.min(i) = sensitivity_analysis_dG;
            current_MILP_bounds_results.dG.max(i) = sensitivity_analysis_dG;

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
            sensitivity_analysis_low_dG{i} = [];
            for(j=1:length(model_net_fluxes.rxns))
                numerator=0;
                denominator=0;
                numerator(model_thermodynamics.product_indexes{j})      = (sensitivity_analysis_low_concentrations{i}(model_thermodynamics.product_indexes{j}));
                numerator=sum(numerator);
                denominator(model_thermodynamics.reactant_indexes{j})    = (sensitivity_analysis_low_concentrations{i}(model_thermodynamics.reactant_indexes{j}));            
                denominator=sum(denominator);
                sensitivity_analysis_low_dG{i}(end+1,1)         = model_thermodynamics.delta_G0(j)+RT*(numerator-denominator);        

            end                

            sensitivity_analysis_index = sensitivity_analysis_index+1;
            consecutive_counter = 0;
        end    
    end

    last_sensitivity_analysis_left{i} = sensitivity_analysis_dG;
    if(last_sensitivity_analysis_left{i} < min(dG_sensitivity_analysis_min_max{i}(:,1))+EPSILON_VALUE)
        last_sensitivity_analysis_left{i} = min(dG_sensitivity_analysis_min_max{i}(:,1));
    end        
end


for(i=1:length(model_net_fluxes.rxns))
    res{i}.flux_index=i;
    res{i}.low_dG          = last_sensitivity_analysis_left{i};
    res{i}.high_dG         = last_sensitivity_analysis_right{i};
    res{i}.net_fluxes.low  = sensitivity_analysis_low_fluxes{i};
    res{i}.net_fluxes.high = sensitivity_analysis_high_fluxes{i};
    res{i}.fb_fluxes.low   = sensitivity_analysis_low_fluxes_fb{i};
    res{i}.fb_fluxes.high  = sensitivity_analysis_high_fluxes_fb{i};    
    res{i}.concentrations.low   = sensitivity_analysis_low_concentrations{i};
    res{i}.concentrations.high  = sensitivity_analysis_high_concentrations{i};  
    res{i}.dG.low               = sensitivity_analysis_low_dG{i};
    res{i}.dG.high              = sensitivity_analysis_high_dG{i};  
end

sensitiviy_analysis_dG=res;
save('mat_files/sensitiviy_analysis_dG.mat','sensitiviy_analysis_dG');




