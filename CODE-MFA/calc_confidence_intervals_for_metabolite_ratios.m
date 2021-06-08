% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% Calculate confidence intervals for metabolite ratios
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
clear all;
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


co_factors{1}   = 'NAD_CY/NADH_CY';
co_factors{2}   = 'NADP_CY/NADPH_CY';
co_factors{3}   = 'NAD_MT/NADH_MT';
co_factors{4}   = 'NADP_MT/NADPH_MT';
co_factors{5}   = 'GDP_CY/GTP_CY';
co_factors{6}   = 'GDP_MT/GTP_MT';
co_factors{7}   = 'ATP_CY/ADP_CY';
co_factors{8}   = 'ATP_MT/ADP_MT';    
co_factors{9}   = 'ADP_CY/ADP_MT';    
co_factors{10}  = 'ATP_CY/ATP_MT';    
co_factors{11}  = 'NAD_CY/NAD_MT';    
co_factors{12}  = 'NADH_CY/NADH_MT';    
co_factors{13}  = 'NADP_CY/NADP_MT';    
co_factors{14}  = 'NADPH_CY/NADPH_MT';    
co_factors{15}  = 'GDP_CY/GDP_MT';
co_factors{16}  = 'GTP_CY/GTP_MT';        
co_factors{17}  = 'Glutamine_CY/Glutamine_MT';    
co_factors{18}  = 'Glutamate_CY/Glutamate_MT';    
co_factors{19}  = 'Orthophosphate_CY/Orthophosphate_MT';    
co_factors{20}  = 'Acetyl_CoA_CY/Acetyl_CoA_MT';    
co_factors{21}  = 'OAA_CY/OAA_MT';
co_factors{22}  = 'Citrate_CY/Citrate_MT';
co_factors{23}  = 'CoA_CY/CoA_MT';
co_factors{24}  = 'Fumarate_CY/Fumarate_MT';
co_factors{25}  = 'Malate_CY/Malate_MT'; 
co_factors{26}  = 'Pyruvate_CY/Pyruvate_MT';
co_factors{27}  = 'PEP_CY/PEP_MT';
co_factors{28}  = 'AKG_CY/AKG_MT';
co_factors{29}  = 'Aspartate_CY/Aspartate_MT';
co_factors{30}  = 'Glutathione_CY/Glutathione_MT';
co_factors{31}  = 'Alanine_CY/Alanine_MT';        


% sort all directionalities per their error, in order to start
% sensitivity analysis with the order of the errors. This is supposed
% to be faster as the convergence is faster
[val ind]=sort(directionalities.errors);
directionalities.errors=directionalities.errors(ind);
directionalities.directionality_matrix=directionalities.directionality_matrix(:,ind);
directionalities.predicted_net_fluxes = directionalities.predicted_net_fluxes(:,ind);
directionalities.predicted_fb_fluxes = directionalities.predicted_fb_fluxes(:,ind);
directionalities.predicted_concentrations = directionalities.predicted_concentrations(:,ind);


SENSITIVITY_ANALYSIS_JUMPS_FOR_COFACTOR_RATIO = [0.5 1 1.5 2 2.5 3 3.5 4 4.5 5 5.5 6 6.5 7 7.5 8 8.5 9 9.5 10];            


% find the index of the best score among all directionalities
best_score = min(directionalities.errors);
index_best_score = find(directionalities.errors==min(directionalities.errors));


% predicted values of best score over all directionalities
best_score_predicted_fluxes_fb  = directionalities.predicted_fb_fluxes(:,index_best_score);    
best_score_predicted_net_fluxes = directionalities.predicted_net_fluxes(:,index_best_score);
best_score_predicted_concentrations     = directionalities.predicted_concentrations(:,index_best_score);    
best_score_directionalities             = directionalities.directionality_matrix(:,index_best_score);    


add_co_factors_ratio_to_optimization     = cell(length(co_factors),1);
best_score_predicted_cofactor_ratio      = cell(length(co_factors),1);
sensitivity_analysis_high_fluxes         = cell(length(co_factors),1);
sensitivity_analysis_low_fluxes          = cell(length(co_factors),1);
sensitivity_analysis_high_fluxes_fb      = cell(length(co_factors),1);
sensitivity_analysis_low_fluxes_fb       = cell(length(co_factors),1);
sensitivity_analysis_high_concentrations = cell(length(co_factors),1);
sensitivity_analysis_low_concentrations  = cell(length(co_factors),1);



all_indices_within_confidence_intervals = find(directionalities.errors < (min(directionalities.errors)+CONSTANT_VALUE_FOR_CONFIDENCE_INTERVAL));
% change confidence interval value to be 10% of the best score rather than
% chi square one degree of freedom (as this is a ratio and not one value)
CONSTANT_VALUE_FOR_CONFIDENCE_INTERVAL = best_score/10;

all_indices_within_confidence_intervals_predicted_cofactor_ratio = cell(length(co_factors),length(all_indices_within_confidence_intervals));

cofactor_ratio_sensitivity_analysis_min_max = cell(length(co_factors),1);
res = cell(length(co_factors),1);

for(i=1:length(all_indices_within_confidence_intervals))    
    current_directionalities = directionalities.directionality_matrix(:,all_indices_within_confidence_intervals(i)); 
    MILP_bounds_results_all_indices_within_confidence_intervals{i} = milp_find_bounds(model_net_fluxes, model_thermodynamics, current_directionalities);

    for(j=1:length(co_factors))
        metabolites = strsplit(co_factors{j},'/');
        index_metabolite_1=strcmp(model_thermodynamics.mets,metabolites{1});
        index_metabolite_1=find(index_metabolite_1==1);
        index_metabolite_2=strcmp(model_thermodynamics.mets,metabolites{2});
        index_metabolite_2=find(index_metabolite_2==1);
        add_co_factors_ratio_to_optimization{j}.indices=[index_metabolite_1 index_metabolite_2]; 
        best_score_predicted_cofactor_ratio{j} = best_score_predicted_concentrations(add_co_factors_ratio_to_optimization{j}.indices(1))-best_score_predicted_concentrations(add_co_factors_ratio_to_optimization{j}.indices(2));
        all_indices_within_confidence_intervals_predicted_cofactor_ratio{j,i} = directionalities.predicted_concentrations(add_co_factors_ratio_to_optimization{j}.indices(1),all_indices_within_confidence_intervals(i)) - directionalities.predicted_concentrations(add_co_factors_ratio_to_optimization{j}.indices(2),all_indices_within_confidence_intervals(i));
    
        cofactor_ratio_sensitivity_analysis_min_max{j}(end+1,1)  = MILP_bounds_results_all_indices_within_confidence_intervals{i}.ln_C.min(index_metabolite_1)-MILP_bounds_results_all_indices_within_confidence_intervals{i}.ln_C.max(index_metabolite_2);
        cofactor_ratio_sensitivity_analysis_min_max{j}(end,2)    = MILP_bounds_results_all_indices_within_confidence_intervals{i}.ln_C.max(index_metabolite_1)-MILP_bounds_results_all_indices_within_confidence_intervals{i}.ln_C.min(index_metabolite_2);                
        % check if we have lower and upper bound for these ratio as an
        % input
        for(k=1:length((model_thermodynamics.co_factors)))
            if(model_thermodynamics.co_factors{k}.indices == [index_metabolite_1 index_metabolite_2])
                cofactor_ratio_sensitivity_analysis_min_max{j}(end,1) = log(model_thermodynamics.co_factors{k}.ratio_lb);
                cofactor_ratio_sensitivity_analysis_min_max{j}(end,2) = log(model_thermodynamics.co_factors{k}.ratio_ub);
            elseif(model_thermodynamics.co_factors{k}.indices == [index_metabolite_2 index_metabolite_1])
                cofactor_ratio_sensitivity_analysis_min_max{j}(end,1) = log(1/model_thermodynamics.co_factors{k}.ratio_ub);
                cofactor_ratio_sensitivity_analysis_min_max{j}(end,2) = log(1/model_thermodynamics.co_factors{k}.ratio_lb);                
            end
        end
        index_cofactor=strcmp(model_thermodynamics.co_factors,co_factors{j});
        index_cofactor=find(index_cofactor==1);        
        model_thermodynamics.co_factors{index_cofactor};

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
parfor(i=1:length(co_factors))    
    sensitivity_analysis_index = 1;   
    used_fva_indexes=zeros(length(all_indices_within_confidence_intervals),1);    
    consecutive_counter = 0;
    sensitivity_analysis_cofactor_ratio = best_score_predicted_cofactor_ratio{i}+SENSITIVITY_ANALYSIS_JUMPS_FOR_COFACTOR_RATIO(1);
    while(true)      
        consecutive_counter = consecutive_counter+1;
        current_indexes = find((sensitivity_analysis_cofactor_ratio>(cofactor_ratio_sensitivity_analysis_min_max{i}(:,1)+EPSILON_VALUE))&(sensitivity_analysis_cofactor_ratio<(cofactor_ratio_sensitivity_analysis_min_max{i}(:,2)-EPSILON_VALUE))&(used_fva_indexes~=1));
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

        while (sensitivity_analysis_index <= length(SENSITIVITY_ANALYSIS_JUMPS_FOR_COFACTOR_RATIO))

            fprintf('cofactor index=%d +++++++ %d  +++++++ %2f %2f\n', i, sensitivity_analysis_index, best_score_predicted_cofactor_ratio{i}, sensitivity_analysis_cofactor_ratio);   

            my_model_thermodynamics = model_thermodynamics;

            current_predicted_fluxes_fb         = directionalities.predicted_fb_fluxes(:,all_indices_within_confidence_intervals(current_indexes(1)));
            current_predicted_fluxes_net        = directionalities.predicted_net_fluxes(:,all_indices_within_confidence_intervals(current_indexes(1)));
            current_predicted_concentrations    = directionalities.predicted_concentrations(:,all_indices_within_confidence_intervals(current_indexes(1)));

            current_directionalities = directionalities.directionality_matrix(:,all_indices_within_confidence_intervals(current_indexes(1)));             

            initial_fluxes_net = current_predicted_fluxes_net;
            initial_fluxes_fb = current_predicted_fluxes_fb;
            initial_concentrations = current_predicted_concentrations;

            sensitivity_analysis_cofactor_ratio = best_score_predicted_cofactor_ratio{i}+SENSITIVITY_ANALYSIS_JUMPS_FOR_COFACTOR_RATIO(sensitivity_analysis_index);                       

            max_cofactor_ratio = cofactor_ratio_sensitivity_analysis_min_max{i}(current_indexes(1),2);
            min_cofactor_ratio = cofactor_ratio_sensitivity_analysis_min_max{i}(current_indexes(1),1);
            if ((sensitivity_analysis_cofactor_ratio > max_cofactor_ratio-EPSILON_VALUE) || (sensitivity_analysis_cofactor_ratio < min_cofactor_ratio+EPSILON_VALUE))
                break;
            end


            directionality_vector = current_directionalities;
            MILP_bounds_results = MILP_bounds_results_all_indices_within_confidence_intervals{current_indexes(1)};
            current_MILP_bounds_results = MILP_bounds_results;

            add_co_factors_ratio_to_optimization{i}.ratio_lb = exp(sensitivity_analysis_cofactor_ratio);
            add_co_factors_ratio_to_optimization{i}.ratio_ub = exp(sensitivity_analysis_cofactor_ratio);
            my_model_thermodynamics.co_factors{end+1}=add_co_factors_ratio_to_optimization{i};


            mymodel = model;
            try
            tic;                
            [exitflag error error_match_labeling error_match_concentrations predicted_net_flux predicted_fb_flux concentrations number_of_independent_variables number_of_fitted_elements] = ComputeEMUOptFlux(current_directionalities, mymodel, model_net_fluxes, my_model_thermodynamics, EMU, idv, WC_known_metabolites_idv, WC_known_metabolites_concentration, initial_fluxes_fb, initial_fluxes_net, initial_concentrations, current_MILP_bounds_results,best_score+CONSTANT_VALUE_FOR_CONFIDENCE_INTERVAL);
            elapsed_time=toc;            
            catch
                fprintf('ComputeEMUOptFlux failure\n');
                break;
            end

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
    last_sensitivity_analysis_right{i} = sensitivity_analysis_cofactor_ratio;    
    if(last_sensitivity_analysis_right{i} > max(cofactor_ratio_sensitivity_analysis_min_max{i}(:,2))-EPSILON_VALUE)
        last_sensitivity_analysis_right{i} = max(cofactor_ratio_sensitivity_analysis_min_max{i}(:,2));
    end
end

% go to the left
last_sensitivity_analysis_left              = cell(length(model_thermodynamics.mets),1);
parfor(i=1:length(co_factors))
    sensitivity_analysis_index = 1;    
    used_fva_indexes=zeros(length(all_indices_within_confidence_intervals),1);    
    consecutive_counter = 0;
    sensitivity_analysis_cofactor_ratio = best_score_predicted_cofactor_ratio{i}-SENSITIVITY_ANALYSIS_JUMPS_FOR_COFACTOR_RATIO(1);
    while(true)                

        consecutive_counter = consecutive_counter+1;
        current_indexes = find((sensitivity_analysis_cofactor_ratio>(cofactor_ratio_sensitivity_analysis_min_max{i}(:,1)+EPSILON_VALUE))&(sensitivity_analysis_cofactor_ratio<(cofactor_ratio_sensitivity_analysis_min_max{i}(:,2)-EPSILON_VALUE))&(used_fva_indexes~=1));
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

        while (sensitivity_analysis_index <= length(SENSITIVITY_ANALYSIS_JUMPS_FOR_COFACTOR_RATIO))
            fprintf('cofactor index=%d ------- %d  ------- %2f %2f %2f\n', i, sensitivity_analysis_index, sensitivity_analysis_cofactor_ratio, best_score_predicted_cofactor_ratio{i}, last_sensitivity_analysis_right{i});   

            my_model_thermodynamics = model_thermodynamics;

            current_predicted_fluxes_fb         = directionalities.predicted_fb_fluxes(:,all_indices_within_confidence_intervals(current_indexes(1)));
            current_predicted_fluxes_net        = directionalities.predicted_net_fluxes(:,all_indices_within_confidence_intervals(current_indexes(1)));            
            current_predicted_concentrations    = directionalities.predicted_concentrations(:,all_indices_within_confidence_intervals(current_indexes(1)));

            current_directionalities = directionalities.directionality_matrix(:,all_indices_within_confidence_intervals(current_indexes(1))); 
            initial_fluxes_net = current_predicted_fluxes_net;
            initial_fluxes_fb = current_predicted_fluxes_fb;
            initial_concentrations = current_predicted_concentrations;

            sensitivity_analysis_cofactor_ratio = best_score_predicted_cofactor_ratio{i}-SENSITIVITY_ANALYSIS_JUMPS_FOR_COFACTOR_RATIO(sensitivity_analysis_index);           

            max_cofactor_ratio = cofactor_ratio_sensitivity_analysis_min_max{i}(current_indexes(1),2);
            min_cofactor_ratio = cofactor_ratio_sensitivity_analysis_min_max{i}(current_indexes(1),1);
            if ((sensitivity_analysis_cofactor_ratio > max_cofactor_ratio-1e-5) || (sensitivity_analysis_cofactor_ratio < min_cofactor_ratio+1e-5))
                break;
            end



            directionality_vector = current_directionalities;
            MILP_bounds_results = MILP_bounds_results_all_indices_within_confidence_intervals{current_indexes(1)};
            current_MILP_bounds_results = MILP_bounds_results;

            add_co_factors_ratio_to_optimization{i}.ratio_lb = exp(sensitivity_analysis_cofactor_ratio);
            add_co_factors_ratio_to_optimization{i}.ratio_ub = exp(sensitivity_analysis_cofactor_ratio);
            my_model_thermodynamics.co_factors{end+1}=add_co_factors_ratio_to_optimization{i};



            mymodel = model;
%             mymodel.positive_direction_lb = model.positive_direction_lb*0;
%             mymodel.positive_direction_ub = model.positive_direction_ub*1000;            

            try
            tic;
            [exitflag error error_match_labeling error_match_concentrations predicted_net_flux predicted_fb_flux concentrations number_of_independent_variables number_of_fitted_elements] = ComputeEMUOptFlux(current_directionalities, mymodel, model_net_fluxes, my_model_thermodynamics, EMU, idv, WC_known_metabolites_idv, WC_known_metabolites_concentration, initial_fluxes_fb, initial_fluxes_net, initial_concentrations, current_MILP_bounds_results,best_score+CONSTANT_VALUE_FOR_CONFIDENCE_INTERVAL);
            elapsed_time=toc;
            catch
                fprintf('ComputeEMUOptFlux failure\n');
                break;
            end
            

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
    last_sensitivity_analysis_left{i} = sensitivity_analysis_cofactor_ratio;
    if(last_sensitivity_analysis_left{i} < min(cofactor_ratio_sensitivity_analysis_min_max{i}(:,1))+EPSILON_VALUE)
        last_sensitivity_analysis_left{i} = min(cofactor_ratio_sensitivity_analysis_min_max{i}(:,1));
    end        
end

for(i=1:length(co_factors))
    res{i}.co_factors_name           = co_factors{i};
    res{i}.low_cofactor_ratio        = exp(last_sensitivity_analysis_left{i});
    res{i}.high_cofactor_ratio       = exp(last_sensitivity_analysis_right{i});
    res{i}.net_fluxes.low  = sensitivity_analysis_low_fluxes{i};
    res{i}.net_fluxes.high = sensitivity_analysis_high_fluxes{i};
    res{i}.fb_fluxes.low   = sensitivity_analysis_low_fluxes_fb{i};
    res{i}.fb_fluxes.high  = sensitivity_analysis_high_fluxes_fb{i};    
    res{i}.concentrations.low   = sensitivity_analysis_low_concentrations{i};
    res{i}.concentrations.high  = sensitivity_analysis_high_concentrations{i};        
end


sensitiviy_analysis_cofactors_ratio=res;
save('mat_files/sensitiviy_analysis_cofactors_ratio.mat','sensitiviy_analysis_cofactors_ratio');
