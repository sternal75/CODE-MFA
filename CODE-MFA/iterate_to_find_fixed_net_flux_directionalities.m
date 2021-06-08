% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% Step I and II - see figure 2 in the paper
% Step I - Determine bounds on the direction of new flux based on thermodynamic analysis
% Step II - Determine bounds on the direction of new flux based on thermodynamic analysis
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
function [iterations_result] = iterate_to_find_fixed_net_flux_directionalities(model, model_net_fluxes, model_thermodynamics, EMU, idv, WC_known_metabolites_idv, WC_known_metabolites_concentration)
load_constants
outputExcelFileName = ['output/results_MILP_MFA_iterations ' strrep(datestr(datetime),':','-') ['.xlsx']];
delete(outputExcelFileName);
excel_column={'B','E','H','K','N','Q','T','W','Z','AC','AF','AI','AL','AO','AR','AU','AX','BA','BD','BG','BJ','BM','BP','BS','BV','BY','CB','CE','CH','CK','CN','CQ','CT','CW','CZ'};
epsilon = 0.0001;

SENSITIVITY_ANALYSIS_VALUES = [-1000 -500 -400 -300 -200 -100 -50 -20 -5 -1 -0.1 0.1 1 5 20 50 100 200 300 400 500 1000];

CONSTANT_VALUE_ABOVE_MIN_FOR_NON_ACCEPTABLE_SCORE   = 6;
CONSTANT_VALUE_FOR_NON_ACCEPTABLE_SCORE   = 130;
MAXIMUM_OPTIMIZATION_FIRST_STEP_THRESHOLD = 1000;
CONSTANT_VALUE_FOR_BEST_SCORE_ITERATIONS  = 2;


dG_same_direction_array       = zeros(length(model_thermodynamics.rxns),1);
net_flux_same_direction_array = zeros(length(model_thermodynamics.rxns),1);
known_net_flux_directions=nan(length(model_thermodynamics.rxns),1);


MILP_bounds_results_cells=cell(0);
directionality_array=[];
MILP_bounds_results_cells_dG=cell(0);
directionality_array_dG=[];
num_of_iterations=1;

rxns_fb=cell(0);
for(i=1:length(model_net_fluxes.is_net_flux))
    if(model_net_fluxes.is_net_flux(i))
        rxns_fb(end+1) = model_thermodynamics.full_rxns(i);
        rxns_fb{end} = [rxns_fb{end},' (f)'];
        rxns_fb(end+1) = model_thermodynamics.full_rxns(i);
        rxns_fb{end} = [rxns_fb{end},' (b)'];
    else
        rxns_fb(end+1) = model_thermodynamics.full_rxns(i);
    end
end
xlswrite(outputExcelFileName,model_thermodynamics.full_rxns','vf minus vb','A2');
xlswrite(outputExcelFileName,model_thermodynamics.mets,'ln C','A2');  
xlswrite(outputExcelFileName,model_thermodynamics.full_rxns','dG','A2');  
xlswrite(outputExcelFileName,model_thermodynamics.full_rxns','directionalities','A2');  
xlswrite(outputExcelFileName,model_thermodynamics.full_rxns','vf minus vb - dG','A2');  
xlswrite(outputExcelFileName,model_thermodynamics.mets,'ln C - dG','A2');  
xlswrite(outputExcelFileName,model_thermodynamics.full_rxns','dG - dG','A2');  
xlswrite(outputExcelFileName,model_thermodynamics.full_rxns','directionalities - dG','A2');  
xlswrite(outputExcelFileName,rxns_fb','best fluxes fb','A2');  
xlswrite(outputExcelFileName,model_thermodynamics.full_rxns','best fluxes net','A2');  
xlswrite(outputExcelFileName,model_thermodynamics.mets,'best concentrations','A2');  


model_net_fluxes_iterator = model_net_fluxes;
while(1)
    try
    % Find bounds for fluxes, concentrations, and Gibbs free energies
    MILP_bounds_results = milp_find_bounds(model_net_fluxes_iterator, model_thermodynamics, known_net_flux_directions);
    catch
        fprintf('ERROR - milp_find_bounds');
    end
    for(i=1:length(model_thermodynamics.rxns))
        dG_min_max_same_direction = sign(MILP_bounds_results.dG.min(i) * MILP_bounds_results.dG.max(i));
        net_flux_min_max_same_direction = sign(MILP_bounds_results.vf_minus_vb.min(i) * MILP_bounds_results.vf_minus_vb.max(i));
        dG_min_max_same_direction(dG_min_max_same_direction==0)=1;
        dG_min_max_same_direction(dG_min_max_same_direction==-1)=0;
        net_flux_min_max_same_direction(net_flux_min_max_same_direction==0)=1;
        net_flux_min_max_same_direction(net_flux_min_max_same_direction==-1)=0;
        if(model_net_fluxes_iterator.is_net_flux(i))
            dG_same_direction_array(i)          = dG_min_max_same_direction;
            net_flux_same_direction_array(i)    = net_flux_min_max_same_direction;
        else
            if((~dG_min_max_same_direction)||(~net_flux_min_max_same_direction))
                fprintf('ERROR %d\n',i);
            end
            dG_same_direction_array(i)          = dG_min_max_same_direction;
            net_flux_same_direction_array(i)    = net_flux_min_max_same_direction;        
        end
    end
    if(sum(abs(dG_same_direction_array-net_flux_same_direction_array))~=0)
        fprintf('ERROR %d\n',i);        
    end
   
    
    % update directionality vector to find know directions (net flux min and max > 0 or
    % net flux min and max < 0. Directionality vector will be 1 of positive net
    % flux, 0 if negative net flux, and nan if the direction is unknown. 
    directionality_vector=nan(length(model_thermodynamics.rxns),1);
    directionality_vector(net_flux_same_direction_array==1)=sign((MILP_bounds_results.vf_minus_vb.min(net_flux_same_direction_array==1))+(MILP_bounds_results.vf_minus_vb.max(net_flux_same_direction_array==1)));
    directionality_vector(directionality_vector==-1)=0;
    MILP_bounds_results = milp_find_bounds(model_net_fluxes_iterator, model_thermodynamics, directionality_vector);
    MILP_bounds_results_cells{end+1}=MILP_bounds_results;
    directionality_array=[directionality_array directionality_vector];

    
    % output vf_minus_vb to excel file
    arrayForExcelOutput=[MILP_bounds_results.vf_minus_vb.min' MILP_bounds_results.vf_minus_vb.max'];
    xlswrite(outputExcelFileName,{'min','max'},'vf minus vb',excel_column{num_of_iterations}); 
    xlswrite(outputExcelFileName,arrayForExcelOutput,'vf minus vb',strcat(excel_column{num_of_iterations},'2'));     
    % output ln_C to excel file
    arrayForExcelOutput=[MILP_bounds_results.ln_C.min' MILP_bounds_results.ln_C.max'];
    xlswrite(outputExcelFileName,{'min','max'},'ln C',excel_column{num_of_iterations}); 
    xlswrite(outputExcelFileName,arrayForExcelOutput,'ln C',strcat(excel_column{num_of_iterations},'2')); 
    % output dG to excel file
    arrayForExcelOutput=[MILP_bounds_results.dG.min' MILP_bounds_results.dG.max'];
    xlswrite(outputExcelFileName,{'min','max'},'dG',excel_column{num_of_iterations}); 
    xlswrite(outputExcelFileName,arrayForExcelOutput,'dG',strcat(excel_column{num_of_iterations},'2')); 
    % output directionalities
    xlswrite(outputExcelFileName,{'directionalities'},'directionalities',excel_column{num_of_iterations}); 
    xlswrite(outputExcelFileName,directionality_vector,'directionalities',strcat(excel_column{num_of_iterations},'2')); 

    
    % If its not the first iteration
    if(num_of_iterations~=1)
        is_no_more_directionalities = isequalwithequalnans(directionality_array(:,end),directionality_array(:,end-1));
         if(is_no_more_directionalities)
            iterations_result.MILP_bounds_results=MILP_bounds_results;
            iterations_result.directionalities=directionality_array(:,end);
            return;
        end
    end
    

    
    current_MILP_bounds_results = MILP_bounds_results;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % find best score result
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    NUM_OF_TRIALS=100;
    error_matrix=[];
    error_labeling_matrix=cell(1,NUM_OF_TRIALS);
    error_concentrations_matrix=cell(1,NUM_OF_TRIALS);
    exitflag_matrix=[];
    predicted_flux_matrix=[];
    predicted_net_flux_matrix=[];
    concentration_matrix=[];
      

    tic
    % Start from random fluxes and concentrations to find global minimim
    % for the non convex optimization
    parfor(i=1:NUM_OF_TRIALS)
        i

        initial_fluxes_fb = rand(length(model.rxns),1).*(model.positive_direction_ub-model.positive_direction_lb)+model.positive_direction_lb;  % Initial flux vector for non-convex optimization
        initial_fluxes = rand(length(model_net_fluxes_iterator.rxns),1).*(current_MILP_bounds_results.vf_minus_vb.max'-current_MILP_bounds_results.vf_minus_vb.min')+current_MILP_bounds_results.vf_minus_vb.min';  % Initial flux vector for non-convex optimization
        initial_concentrations = rand(length(model_thermodynamics.mets),1).*(current_MILP_bounds_results.ln_C.max'-current_MILP_bounds_results.ln_C.min')+current_MILP_bounds_results.ln_C.min';  % Initial metabolites concentrations
    
        [exitflag error error_match_labeling error_match_concentrations predicted_net_flux predicted_flux concentrations number_of_independent_variables number_of_fitted_elements] = ComputeEMUOptFlux(directionality_vector, model, model_net_fluxes_iterator, model_thermodynamics, EMU, idv, WC_known_metabolites_idv, WC_known_metabolites_concentration, initial_fluxes_fb, initial_fluxes, initial_concentrations, current_MILP_bounds_results);
        fprintf('%d - error=%d\n',i,error);
        error_matrix(i) = error;
        error_labeling_matrix{i} = error_match_labeling;
        error_concentrations_matrix{i} = error_match_concentrations;
        exitflag_matrix(i) = exitflag;
        predicted_flux_matrix(:,i)  = predicted_flux;
        predicted_net_flux_matrix(:,i)  = predicted_net_flux;
        concentration_matrix(:,i)   = concentrations;
        
    end
    toc
    MFA_results.error_matrix=error_matrix;
    MFA_results.exitflag_matrix=exitflag_matrix;
    MFA_results.predicted_flux_matrix=predicted_flux_matrix;
    MFA_results.predicted_net_flux_matrix=predicted_net_flux_matrix;
    MFA_results.concentration_matrix=concentration_matrix;


    best_error=min(error_matrix);
    best_error_index = find(error_matrix==best_error);
    MAX_ACCEPTABLE_SCORE=max(CONSTANT_VALUE_FOR_NON_ACCEPTABLE_SCORE,best_error+CONSTANT_VALUE_ABOVE_MIN_FOR_NON_ACCEPTABLE_SCORE);  

    % save the fluxes with the best score among all optimizations
    initial_fluxes          = predicted_net_flux_matrix(:,best_error_index);

    
    % output best fb flux predictions to excel file
    arrayForExcelOutput=predicted_flux_matrix(:,best_error_index);
    xlswrite(outputExcelFileName,{'best fluxes fb'},'best fluxes fb',excel_column{num_of_iterations}); 
    xlswrite(outputExcelFileName,arrayForExcelOutput,'best fluxes fb',strcat(excel_column{num_of_iterations},'2')); 
    % output best net flux predictions to excel file
    arrayForExcelOutput=predicted_net_flux_matrix(:,best_error_index);
    xlswrite(outputExcelFileName,{'best fluxes net'},'best fluxes net',excel_column{num_of_iterations}); 
    xlswrite(outputExcelFileName,arrayForExcelOutput,'best fluxes net',strcat(excel_column{num_of_iterations},'2')); 
    % output best ln concentrations predictions to excel file
    arrayForExcelOutput=concentration_matrix(:,best_error_index);
    xlswrite(outputExcelFileName,{'best concentrations'},'best concentrations',excel_column{num_of_iterations}); 
    xlswrite(outputExcelFileName,arrayForExcelOutput,'best concentrations',strcat(excel_column{num_of_iterations},'2')); 
    % output best error
    arrayForExcelOutput=[best_error exitflag_matrix(best_error_index)];
    xlswrite(outputExcelFileName,{'best error', 'best exitflag'},'best error',excel_column{num_of_iterations}); 
    xlswrite(outputExcelFileName,arrayForExcelOutput,'best error',strcat(excel_column{num_of_iterations},'2'));
    % output errors of random iterations
    arrayForExcelOutput=[error_matrix' exitflag_matrix'];
    xlswrite(outputExcelFileName,{'error', 'exitflag'},'error',excel_column{num_of_iterations}); 
    xlswrite(outputExcelFileName,arrayForExcelOutput,'error',strcat(excel_column{num_of_iterations},'2')); 
       
    
    % take sll results with accepted error
    [sorted_error sorted_error_indicess] = sort(error_matrix);
    accepted_score_indices          = sorted_error_indicess(sorted_error < MAX_ACCEPTABLE_SCORE)
    accepted_score_fluxes_fb        = predicted_flux_matrix(:,accepted_score_indices);
    accepted_score_fluxes           = predicted_net_flux_matrix(:,accepted_score_indices);
    accepted_score_concentrations   = concentration_matrix(:,accepted_score_indices);
    accepted_score_error            = error_matrix(accepted_score_indices);
    accepted_score_error_labeling   = error_labeling_matrix(accepted_score_indices);  
    
       
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % sensitivity analysis - drop large error results based on fluxes
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
    num_of_net_fluxes = length(directionality_vector);
    sensitivity_analysis_last_positive_value_vector = MILP_bounds_results.vf_minus_vb.max;
    sensitivity_analysis_last_negative_value_vector(num_of_net_fluxes+1:2*num_of_net_fluxes) = MILP_bounds_results.vf_minus_vb.min;
    parfor(i=1:2*num_of_net_fluxes)
        if(i<=length(directionality_vector))
            if(model_net_fluxes.skip_sensitivity_analysis_reaction_indices(i)==1)
                continue;
            end
            % if the directionality is not determined yet, perform sensitivity
            % analysis to the positive and negative directions
            first_sensitivity_analysis_value_to_the_right = max(accepted_score_fluxes(i,:)); 
            j=1;
            sensitiviy_analysis_values_to_the_right = SENSITIVITY_ANALYSIS_VALUES(first_sensitivity_analysis_value_to_the_right < SENSITIVITY_ANALYSIS_VALUES);
            sensitivity_analysis_last_positive_value_vector(i) = first_sensitivity_analysis_value_to_the_right;
            last_flux_sensitivity_analysis_value_positive = first_sensitivity_analysis_value_to_the_right;
            for(k=1:length(accepted_score_indices))
                current_initial_concentrations  = accepted_score_concentrations(:,k);
                current_initial_fluxes          = accepted_score_fluxes(:,k);
                current_initial_fluxes_fb       = accepted_score_fluxes_fb(:,k);
                while (j < length(sensitiviy_analysis_values_to_the_right))
                    % try to go to the negativ direction            
                    current_MILP_bounds_results = MILP_bounds_results;
                    current_flux_sensitivity_analysis_value=sensitiviy_analysis_values_to_the_right(j);
                    sensitivity_analysis_last_positive_value_vector(i) = current_flux_sensitivity_analysis_value;
                    if((current_flux_sensitivity_analysis_value) >= MILP_bounds_results.vf_minus_vb.max(i)-epsilon)
                        break;
                    end

                    current_initial_fluxes(i)=(current_flux_sensitivity_analysis_value);
                    model_net_fluxes_iterator_current = model_net_fluxes_iterator;
                    current_directionality_vector = directionality_vector;
                    if(current_flux_sensitivity_analysis_value>=0)
                        current_directionality_vector(i)=1;
                        model_net_fluxes_iterator_current.positive_direction_lb(i) = current_flux_sensitivity_analysis_value-epsilon;
                        model_net_fluxes_iterator_current.positive_direction_ub(i) = current_flux_sensitivity_analysis_value+epsilon;
                    else
                        current_directionality_vector(i)=0;
                        model_net_fluxes_iterator_current.negative_direction_lb(i) = abs(current_flux_sensitivity_analysis_value)-epsilon;
                        model_net_fluxes_iterator_current.negative_direction_ub(i) = abs(current_flux_sensitivity_analysis_value)+epsilon; 
                    end

                    % find all bounds as it may change due to known directionalities and fixed
                    % values asigned by the sensitivity analysis
                    i
                    current_flux_sensitivity_analysis_value
                    current_MILP_bounds_results = milp_find_bounds(model_net_fluxes_iterator_current, model_thermodynamics, current_directionality_vector, current_MILP_bounds_results);
                    net_flux_same_direction_array = zeros(length(model_thermodynamics.rxns),1);
                    for(ind=1:length(model_thermodynamics.rxns))
                        net_flux_min_max_same_direction = sign(current_MILP_bounds_results.vf_minus_vb.min(ind) * current_MILP_bounds_results.vf_minus_vb.max(ind));
                        net_flux_min_max_same_direction(net_flux_min_max_same_direction==0)=1;
                        net_flux_min_max_same_direction(net_flux_min_max_same_direction==-1)=0;
                        if(model_net_fluxes_iterator_current.is_net_flux(ind))
                            net_flux_same_direction_array(ind)    = net_flux_min_max_same_direction;
                        else
                            net_flux_same_direction_array(ind)    = net_flux_min_max_same_direction;        
                        end
                    end
                    current_directionality_vector=nan(length(model_thermodynamics.rxns),1);
                    current_directionality_vector(net_flux_same_direction_array==1)=sign((current_MILP_bounds_results.vf_minus_vb.min(net_flux_same_direction_array==1))+(current_MILP_bounds_results.vf_minus_vb.max(net_flux_same_direction_array==1)));
                    current_directionality_vector(current_directionality_vector==-1)=0;

                    [exitflag error error_match_labeling error_match_concentrations predicted_net_flux predicted_flux concentrations number_of_independent_variables number_of_fitted_elements] = ComputeEMUOptFlux(current_directionality_vector, model, model_net_fluxes_iterator_current, model_thermodynamics, EMU, idv, WC_known_metabolites_idv, WC_known_metabolites_concentration, current_initial_fluxes_fb, current_initial_fluxes, current_initial_concentrations, current_MILP_bounds_results, MAX_ACCEPTABLE_SCORE);
                    if((error > MAX_ACCEPTABLE_SCORE) || ((exitflag ~= 1)&&(exitflag ~= 2)&&(exitflag ~= -3)&&(exitflag ~= 0)))
                        fprintf('exitflag=%d error=%d best_error=%d\n', exitflag, error, best_error);  
                        if(current_flux_sensitivity_analysis_value==0.1)
                            sensitivity_analysis_last_positive_value_vector(i)=last_flux_sensitivity_analysis_value_positive;
                        end
                        break;
                    end
                    fprintf('net flux +++++++ %d:%d:%d  +++++++ %2f %2f\n', i, j, k, initial_fluxes(i), current_flux_sensitivity_analysis_value);   
                    j=j+1;  
                    last_flux_sensitivity_analysis_value_positive = current_flux_sensitivity_analysis_value;
                end                    
            end
        else
            if(model_net_fluxes.skip_sensitivity_analysis_reaction_indices(i-num_of_net_fluxes)==1)
                continue;
            end
            % if the directionality is not determined yet, perform sensitivity
            % analysis to the positive and negative directions

            j=1;
            first_sensitivity_analysis_value_to_the_left = min(accepted_score_fluxes(i-num_of_net_fluxes,:));
            sensitiviy_analysis_values_to_the_left = SENSITIVITY_ANALYSIS_VALUES(first_sensitivity_analysis_value_to_the_left > SENSITIVITY_ANALYSIS_VALUES);
            sensitiviy_analysis_values_to_the_left = flipud(sensitiviy_analysis_values_to_the_left')';
            sensitivity_analysis_last_negative_value_vector(i) = first_sensitivity_analysis_value_to_the_left;
            last_flux_sensitivity_analysis_value_negative = first_sensitivity_analysis_value_to_the_left;
            for(k=1:length(accepted_score_indices))        
                current_initial_concentrations  = accepted_score_concentrations(:,k);
                current_initial_fluxes          = accepted_score_fluxes(:,k);
                current_initial_fluxes_fb       = accepted_score_fluxes_fb(:,k);                    
                while (j <= length(sensitiviy_analysis_values_to_the_left))

                    % try to go to the negativ direction            
                    current_MILP_bounds_results = MILP_bounds_results;
                    current_flux_sensitivity_analysis_value = sensitiviy_analysis_values_to_the_left(j);

                    sensitivity_analysis_last_negative_value_vector(i) = current_flux_sensitivity_analysis_value;
                    if((current_flux_sensitivity_analysis_value) <= (MILP_bounds_results.vf_minus_vb.min(i-num_of_net_fluxes)+epsilon))
                        break;
                    end
                    current_initial_fluxes(i-num_of_net_fluxes)=(current_flux_sensitivity_analysis_value);
                    model_net_fluxes_iterator_current = model_net_fluxes_iterator;
                    current_directionality_vector = directionality_vector;

                    if(current_flux_sensitivity_analysis_value>=0)
                        current_directionality_vector(i-num_of_net_fluxes)=1;
                        model_net_fluxes_iterator_current.positive_direction_lb(i-num_of_net_fluxes) = current_flux_sensitivity_analysis_value-epsilon;
                        model_net_fluxes_iterator_current.positive_direction_ub(i-num_of_net_fluxes) = current_flux_sensitivity_analysis_value+epsilon;                           
                    else
                        current_directionality_vector(i-num_of_net_fluxes)=0;
                        model_net_fluxes_iterator_current.negative_direction_lb(i-num_of_net_fluxes) = abs(current_flux_sensitivity_analysis_value)-epsilon;
                        model_net_fluxes_iterator_current.negative_direction_ub(i-num_of_net_fluxes) = abs(current_flux_sensitivity_analysis_value)+epsilon; 
                    end

                    % find all bounds as it may change due to known directionalities and fixed
                    % values asigned by the sensitivity analysis
                    current_MILP_bounds_results = milp_find_bounds(model_net_fluxes_iterator_current, model_thermodynamics, current_directionality_vector, current_MILP_bounds_results);
                    net_flux_same_direction_array = zeros(length(model_thermodynamics.rxns),1);
                    for(ind=1:length(model_thermodynamics.rxns))
                        net_flux_min_max_same_direction = sign(current_MILP_bounds_results.vf_minus_vb.min(ind) * current_MILP_bounds_results.vf_minus_vb.max(ind));
                        net_flux_min_max_same_direction(net_flux_min_max_same_direction==0)=1;
                        net_flux_min_max_same_direction(net_flux_min_max_same_direction==-1)=0;
                        if(model_net_fluxes_iterator_current.is_net_flux(ind))
                            net_flux_same_direction_array(ind)    = net_flux_min_max_same_direction;
                        else
                            net_flux_same_direction_array(ind)    = net_flux_min_max_same_direction;        
                        end
                    end
                    current_directionality_vector=nan(length(model_thermodynamics.rxns),1);
                    current_directionality_vector(net_flux_same_direction_array==1)=sign((current_MILP_bounds_results.vf_minus_vb.min(net_flux_same_direction_array==1))+(current_MILP_bounds_results.vf_minus_vb.max(net_flux_same_direction_array==1)));
                    current_directionality_vector(current_directionality_vector==-1)=0;

                    [exitflag error error_match_labeling error_match_concentrations predicted_net_flux predicted_flux concentrations number_of_independent_variables number_of_fitted_elements] = ComputeEMUOptFlux(current_directionality_vector, model, model_net_fluxes_iterator_current, model_thermodynamics, EMU, idv, WC_known_metabolites_idv, WC_known_metabolites_concentration, current_initial_fluxes_fb, current_initial_fluxes, current_initial_concentrations, current_MILP_bounds_results, MAX_ACCEPTABLE_SCORE);
                    if((error > MAX_ACCEPTABLE_SCORE) || ((exitflag ~= 1)&&(exitflag ~= 2)&&(exitflag ~= -3)&&(exitflag ~= 0)))                    
                        fprintf('exitflag=%d error=%d best_error=%d\n', exitflag, error, best_error);                      
                        if(current_flux_sensitivity_analysis_value==-0.1)
                            sensitivity_analysis_last_negative_value_vector(i)=last_flux_sensitivity_analysis_value_negative;
                        end                            
                        break;
                    end
                    fprintf('net flux ------- %d:%d:%d  ------- %2f %2f\n', i-num_of_net_fluxes, j, k, initial_fluxes(i-num_of_net_fluxes), current_flux_sensitivity_analysis_value);   
                    j=j+1;
                    last_flux_sensitivity_analysis_value_negative = current_flux_sensitivity_analysis_value;
                end
            end
        end
    end
    
    last_right_predicted_flux = sensitivity_analysis_last_positive_value_vector';
    last_right_predicted_flux(last_right_predicted_flux>MILP_bounds_results.vf_minus_vb.max')=MILP_bounds_results.vf_minus_vb.max(last_right_predicted_flux>MILP_bounds_results.vf_minus_vb.max');
    model_net_fluxes_iterator.positive_direction_ub(last_right_predicted_flux>0)=last_right_predicted_flux(last_right_predicted_flux>0);
    model_net_fluxes_iterator.negative_direction_lb(last_right_predicted_flux<0)=abs(last_right_predicted_flux(last_right_predicted_flux<0));
    known_net_flux_directions(last_right_predicted_flux<0) = 0;
        
    last_left_predicted_flux = sensitivity_analysis_last_negative_value_vector(num_of_net_fluxes+1:2*num_of_net_fluxes)';
    last_left_predicted_flux(last_left_predicted_flux<MILP_bounds_results.vf_minus_vb.min')=MILP_bounds_results.vf_minus_vb.min(last_left_predicted_flux<MILP_bounds_results.vf_minus_vb.min');
    model_net_fluxes_iterator.positive_direction_lb(last_left_predicted_flux>0)=last_left_predicted_flux(last_left_predicted_flux>0);
    model_net_fluxes_iterator.negative_direction_ub(last_left_predicted_flux<0)=abs(last_left_predicted_flux(last_left_predicted_flux<0));
    known_net_flux_directions(last_left_predicted_flux>0) = 1;
    
    num_of_iterations=num_of_iterations+1;
end

