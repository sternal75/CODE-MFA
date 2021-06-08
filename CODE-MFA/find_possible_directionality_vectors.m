% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% Calculate confidence intervals for metabolite ratios
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
load('mat_files/model.mat','model');
load('mat_files/model_net_fluxes.mat','model_net_fluxes');
load('mat_files/model_thermodynamics.mat','model_thermodynamics');    
load('mat_files/EMU.mat','EMU');
load('mat_files/idv.mat','idv');
load('mat_files/WC_known_metabolites_idv.mat','WC_known_metabolites_idv');
load('mat_files/WC_known_metabolites_concentration.mat','WC_known_metabolites_concentration');
load('mat_files/MILP_results.mat','MILP_results');
load('mat_files/iterations_result.mat','iterations_result');


addpath('./functions/emu') 
addpath('./functions/general') 

% Overcome potential local minima of the non-convex optimization
NUM_OF_TRIALS=100;

% seed for random - fluxes and concentrations
ts=clock;
rand('seed',(ts(6)*10000));

current_MILP_results = MILP_results;
mymodel_net_fluxes = model_net_fluxes;
mymodel = model;
mymodel_net_fluxes = model_net_fluxes;
mymodel_thermodynamics = model_thermodynamics;
myEMU = EMU;
myidv = idv;
myWC_known_metabolites_idv = WC_known_metabolites_idv;
myWC_known_metabolites_concentration = WC_known_metabolites_concentration;

directionality_vector_all=[];
predicted_net_flux_matrix_all=[];
predicted_fb_flux_matrix_all=[];
concentration_matrix_all=[];
error_matrix_all=[];
exitflag_matrix_all=[];

number_of_directionality_vectors = size(current_MILP_results.all_directionalities,2);
parfor(i=1:number_of_directionality_vectors)

    directionality_vector_index = i;

    directionality_vector = current_MILP_results.all_directionalities(:,directionality_vector_index);

    MILP_bounds_results = milp_find_bounds(mymodel_net_fluxes, mymodel_thermodynamics, directionality_vector, iterations_result.MILP_bounds_results);
    current_MILP_bounds_results = MILP_bounds_results;

    error_matrix=[];
    exitflag_matrix=[];
    predicted_net_flux_matrix=[];
    predicted_fb_flux_matrix=[];
    concentration_matrix=[];
    for(j=1:NUM_OF_TRIALS)
        fprintf('**** %d %d\n', i, j)

        initial_fluxes_fb = rand(length(mymodel.rxns),1).*(mymodel.positive_direction_ub-mymodel.positive_direction_lb)+mymodel.positive_direction_lb;  % Initial flux vector for non-convex optimization
        initial_fluxes_net = rand(length(mymodel_net_fluxes.rxns),1).*(current_MILP_bounds_results.vf_minus_vb.max'-current_MILP_bounds_results.vf_minus_vb.min')+current_MILP_bounds_results.vf_minus_vb.min';  % Initial flux vector for non-convex optimization
        initial_concentrations = rand(length(mymodel_thermodynamics.mets),1).*(current_MILP_bounds_results.ln_C.max'-current_MILP_bounds_results.ln_C.min')+current_MILP_bounds_results.ln_C.min';  % Initial metabolites concentrations        

        [exitflag error error_match_labeling error_match_concentrations predicted_net_flux predicted_fb_flux concentrations number_of_independent_variables number_of_fitted_elements] = ComputeEMUOptFlux(directionality_vector, mymodel, mymodel_net_fluxes, mymodel_thermodynamics, myEMU, myidv, myWC_known_metabolites_idv, myWC_known_metabolites_concentration, initial_fluxes_fb, initial_fluxes_net, initial_concentrations, current_MILP_bounds_results);

        error_matrix(j,1) = error;
        exitflag_matrix(j,1) = exitflag;
        predicted_net_flux_matrix(:,j)  = predicted_net_flux;
        predicted_fb_flux_matrix(:,j)  = predicted_fb_flux;
        concentration_matrix(:,j)   = concentrations;

    end
    directionality_vector_all(:,i) = directionality_vector;
    predicted_net_flux_matrix_all(:,:,i)=predicted_net_flux_matrix;    
    predicted_fb_flux_matrix_all(:,:,i)=predicted_fb_flux_matrix;    
    concentration_matrix_all(:,:,i)=concentration_matrix;    
    error_matrix_all(:,i)=error_matrix;
    exitflag_matrix_all(:,i)=exitflag_matrix;




end


results = cell(length(number_of_directionality_vectors),1);
for(i=1:number_of_directionality_vectors)
    results{i}.directionality_vector = directionality_vector_all(:,i);
    results{i}.predicted_net_flux_matrix=predicted_net_flux_matrix_all(:,:,i);    
    results{i}.predicted_fb_flux_matrix=predicted_fb_flux_matrix_all(:,:,i);   
    results{i}.concentration_matrix=concentration_matrix_all(:,:,i);    
    results{i}.error_matrix=error_matrix_all(:,i);
    results{i}.exitflag_matrix=exitflag_matrix_all(:,i);

    [sorted_error sorted_error_indicess] = sort(results{i}.error_matrix);        

    results{i}.best_predicted_net_flux         = results{i}.predicted_net_flux_matrix(:,sorted_error_indicess(1));   
    results{i}.best_predicted_fb_flux          = results{i}.predicted_fb_flux_matrix(:,sorted_error_indicess(1));   
    results{i}.best_predicted_concentration    = results{i}.concentration_matrix(:,sorted_error_indicess(1));   
    results{i}.best_error                      = results{i}.error_matrix(sorted_error_indicess(1));   
    results{i}.best_exitflag                   = results{i}.exitflag_matrix(sorted_error_indicess(1));   
end
directionalities.directionality_matrix=[];
directionalities.predicted_net_fluxes=[];
directionalities.predicted_fb_fluxes=[];
directionalities.errors=[];
directionalities.predicted_concentrations=[];
for(i=1:number_of_directionality_vectors)    
    try
        if((results{i}.best_exitflag~=1)&(results{i}.best_exitflag~=2)&(results{i}.best_exitflag~=-3)&(results{i}.best_exitflag~=0))
            continue;
        end
    catch
        fprintf('job %d did not end\n', i);
        continue;
    end
    directionalities.directionality_matrix      = [directionalities.directionality_matrix results{i}.directionality_vector];
    directionalities.predicted_net_fluxes       = [directionalities.predicted_net_fluxes results{i}.best_predicted_net_flux];
    directionalities.predicted_fb_fluxes        = [directionalities.predicted_fb_fluxes results{i}.best_predicted_fb_flux];
    directionalities.predicted_concentrations   = [directionalities.predicted_concentrations results{i}.best_predicted_concentration];
    directionalities.errors                     = [directionalities.errors results{i}.best_error];
end 
save('mat_files/directionalities.mat','directionalities');


