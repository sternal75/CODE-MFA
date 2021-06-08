% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% Infer fluxes and concentrations by:
% 1. fitting convoluted simulated labeling to whole cell measured labeling
% 2. fitting convoluted simulated metabolite concentrations to whole cell
% measured concentrations
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
function [exitflag error error_match_labeling error_match_concentrations fluxes_net one_direction_flux concentrations number_of_independent_variables number_of_fitted_elements] = ComputeEMUOptFlux(flux_directionalities, model, model_net_fluxes, model_thermodynamics, EMU, idv, WC_known_metabolites_idv, WC_known_metabolites_concentration, initial_fluxes_fb, initial_flux, initial_concentrations, MILP_bounds_results, ObjectiveLimit)
global g_constants  g_EMU g_idv g_idv_known_arr g_idv_known_mat g_iteration_num g_model g_model_net_fluxes g_model_thermodynamics g_WC_known_metabolites_idv g_WC_known_metabolites_concentration g_MILP_bounds_results g_flux_directionalities g_num_of_experiments;

load_constants
g_constants.CY_WC_VOLUME    = CY_WC_VOLUME;
g_constants.MT_WC_VOLUME    = MT_WC_VOLUME;
g_constants.RT              = RT;
g_constants.EPSILON_VALUE = EPSILON_VALUE;
g_constants.WC_CONCENTRATION_STD_FACTOR = WC_CONCENTRATION_STD_FACTOR;



% change matrices according to the directionalities
g_EMU = EMU;
g_idv = idv;
% g_idv_known = idv_known;
% g_idv_known_variance = idv_known_variance;
g_WC_known_metabolites_idv = WC_known_metabolites_idv;
g_WC_known_metabolites_concentration = WC_known_metabolites_concentration;

g_iteration_num = 0;
g_model = model;
g_model_net_fluxes = model_net_fluxes;
g_model_thermodynamics = model_thermodynamics;
g_MILP_bounds_results = MILP_bounds_results;
g_flux_directionalities = flux_directionalities;
g_num_of_experiments = length(g_WC_known_metabolites_idv);



g_idv_known_arr = cell(g_num_of_experiments,1);
g_idv_known_mat = cell(g_num_of_experiments,1);
for (i=1:g_num_of_experiments)
    g_idv_known_arr{i} = zeros(0,1);
    g_idv_known_mat{i}=[];
    for j=1:length(g_WC_known_metabolites_idv{i})
        if(isnan(g_WC_known_metabolites_idv{i}{j}.index_MT))
            alon=1;
        end
        g_idv_known_arr{i} = [g_idv_known_arr{i};g_WC_known_metabolites_idv{i}{j}.index_CY;g_WC_known_metabolites_idv{i}{j}.index_MT];
        g_idv_known_mat{i} = [g_idv_known_mat{i};g_WC_known_metabolites_idv{i}{j}.index_CY g_WC_known_metabolites_idv{i}{j}.index_MT];
    end    
end




% for linear inequality constraints
[met_num, rxn_num] = size(model.S);
Aeq = model_net_fluxes.S(model_net_fluxes.met_extra == 0, :);
% Aeq = [Aeq;model_net_fluxes.equality_constraints];
Beq = zeros(size(Aeq,1),1);
Aeq = [Aeq;model_net_fluxes.non_equality_constraint_fluxes];
Aeq = [Aeq zeros(size(Aeq,1),length(initial_concentrations)+length(initial_fluxes_fb))];

% Aeq = [Aeq zeros(size(Aeq,1),length(g_WC_known_metabolites))];
Beq = [Beq;model_net_fluxes.non_equality_constraint_values];
one_direction_flux_ind = 1;
for(i=1:length(model_thermodynamics.rxns))
    Aeq(end+1,i)  = 1;
    Aeq(end,length(initial_flux)+length(initial_concentrations)+one_direction_flux_ind) = -1;
    Beq(end+1) = 0;                  
    if(model_net_fluxes.is_net_flux(i))
        Aeq(end,length(initial_flux)+length(initial_concentrations)+one_direction_flux_ind+1) = 1;
        one_direction_flux_ind = one_direction_flux_ind+2;
    else
        one_direction_flux_ind=one_direction_flux_ind+1;
    end
    
end

fb_equality_constraints = [zeros(size(model.equality_constraints,1),length(initial_flux)+length(initial_concentrations)) model.equality_constraints];
Aeq = [Aeq;fb_equality_constraints];
Beq = [Beq;zeros(size(fb_equality_constraints,1),1)];
% for linear inequality constraints
A1 = zeros(length(model_thermodynamics.rxns), length(model_thermodynamics.mets));
b1 = zeros(length(model_thermodynamics.rxns),1);
A2 = zeros(length(model_thermodynamics.rxns), length(model_thermodynamics.mets));
b2 = zeros(length(model_thermodynamics.rxns),1);
for(i=1:length(model_thermodynamics.rxns))
    if(model_thermodynamics.thermodynamics_of_reaction_defined(i))
        
        A1(i,model_thermodynamics.product_indexes{i})  = g_constants.RT;
        A1(i,model_thermodynamics.reactant_indexes{i}) = -g_constants.RT;
        b1(i) = -model_thermodynamics.delta_G0(i)+MILP_bounds_results.dG.max(i);         

        A2(i,model_thermodynamics.product_indexes{i})  = -g_constants.RT;
        A2(i,model_thermodynamics.reactant_indexes{i}) = g_constants.RT;
        b2(i) = model_thermodynamics.delta_G0(i)-MILP_bounds_results.dG.min(i);                     

    end
end
A1 = [zeros(size(A1,1),length(initial_flux)) A1 zeros(size(A1,1),length(initial_fluxes_fb))];
A2 = [zeros(size(A2,1),length(initial_flux)) A2 zeros(size(A2,1),length(initial_fluxes_fb))];
A = [A1;A2];
b = [b1;b2];

% % add lb/ub constratins for co-factor ratios
A5 = zeros(length(model_thermodynamics.co_factors), length(model_thermodynamics.mets));
A6 = zeros(length(model_thermodynamics.co_factors), length(model_thermodynamics.mets));
for(i=1:length(model_thermodynamics.co_factors))
    A5(i,model_thermodynamics.co_factors{i}.indices(1)) = 1;
    A5(i,model_thermodynamics.co_factors{i}.indices(2)) = -1;
    b5(i,1)=log(model_thermodynamics.co_factors{i}.ratio_ub);
    A6(i,model_thermodynamics.co_factors{i}.indices(1)) = -1;
    A6(i,model_thermodynamics.co_factors{i}.indices(2)) = 1;    
    b6(i,1)=-log(model_thermodynamics.co_factors{i}.ratio_lb);
end
A5 = [zeros(size(A5,1),length(initial_flux)) A5 zeros(size(A5,1),length(initial_fluxes_fb))];
A6 = [zeros(size(A6,1),length(initial_flux)) A6 zeros(size(A6,1),length(initial_fluxes_fb))];
A=[A;A5;A6];
b=[b;b5;b6];



options = optimset('Algorithm','sqp', 'GradObj','on', 'MaxIter', 1500);
if(nargin==13)
    options = optimset('Algorithm','sqp', 'GradObj','on', 'MaxIter', 400, 'ObjectiveLimit',ObjectiveLimit);
end


%clc;
fprintf('Running fmincon...\n');


lb = [MILP_bounds_results.vf_minus_vb.min';MILP_bounds_results.ln_C.min';model.positive_direction_lb];
ub = [MILP_bounds_results.vf_minus_vb.max';MILP_bounds_results.ln_C.max';model.positive_direction_ub];
f = rand(size([initial_flux;initial_concentrations;initial_fluxes_fb]))-0.5; 
[initial_values_new,fval,exitflag] = linprog(f, A, b, Aeq, Beq, lb, ub);
if(exitflag ~= 1)
    exitflag=nan;
    error=nan;
    predicted_net_flux=nan;
    predicted_fb_flux=nan;
    fluxes_net=nan;
    one_direction_flux=nan;
    concentrations=nan;
    number_of_independent_variables=nan;
    number_of_fitted_elements=nan;
    error_match_labeling=nan;
    error_match_concentrations=nan;
    
    fprintf('ERROR in LINEAR PROGRAMING\n');    
    return;
end

initial_values = [initial_flux;initial_concentrations;initial_fluxes_fb];
[optimization_parameters ,fval,exitflag]= fmincon(@opt_func, initial_values,A,b, Aeq, Beq, lb, ub, [], options);

fluxes_net          = optimization_parameters(1:model_net_fluxes.rxn_num);
concentrations  = optimization_parameters(model_net_fluxes.rxn_num+1:g_model_net_fluxes.rxn_num+length(g_model_thermodynamics.mets));
flux_fb     = optimization_parameters(model_net_fluxes.rxn_num+length(g_model_thermodynamics.mets)+1:end);
one_direction_flux = ComputeOneDirectionFluxes(fluxes_net, concentrations, flux_fb, model_net_fluxes, model_thermodynamics);

number_of_fitted_elements=0;
for (i=1:g_num_of_experiments)
    for j=1:size(g_idv_known_mat{i},1)
        number_of_fitted_elements = number_of_fitted_elements+length(g_WC_known_metabolites_idv{i}{j}.idv(2:end));
    end
end
number_of_independent_variables = rank(Aeq);


fprintf('Finished fmincon\n');

idv_opt = cell(g_num_of_experiments,1);
idv_d   = cell(g_num_of_experiments,1);
fluxes_net_after_force_zero         = cell(g_num_of_experiments,1);
one_direction_flux_after_force_zero = cell(g_num_of_experiments,1);
for(i=1:g_num_of_experiments)
    fluxes_net_after_force_zero{i} = fluxes_net;
    fluxes_net_after_force_zero{i}(g_model_net_fluxes.force_zero_flux{i}==1)=g_model_net_fluxes.positive_direction_lb(g_model_net_fluxes.force_zero_flux{i}==1);
    one_direction_flux_after_force_zero{i} = ComputeOneDirectionFluxes(fluxes_net_after_force_zero{i}, concentrations, flux_fb, g_model_net_fluxes, g_model_thermodynamics);    
    fcn_name = ['ComputeEmuIDV_Opt' int2str(i)];
    [idv_opt{i} idv_d{i}] = feval(fcn_name,g_idv{i},one_direction_flux_after_force_zero{i});
    %[idv_opt{i} idv_d{i} cycle_error] = ComputeEmuIDV(EMU{i}, idv{i}, g_idv_known_arr{i}, one_direction_flux_after_force_zero{i}, fcn_name);
end
% fprintf('\nIDVs:\n');
% for i=1:length(EMU.list)
%     fprintf('\tEMU %s (%d) = %s\n', EMU.name{i}, i, DispIDV(idv_opt{i}));
% end
[e e_match_labeling e_match_concentrations]= ComputeError(g_WC_known_metabolites_idv, idv_opt, g_idv_known_mat, g_WC_known_metabolites_concentration, g_EMU, concentrations, g_num_of_experiments);
error = e;
error_match_labeling = e_match_labeling;
error_match_concentrations = e_match_concentrations;
% fprintf('\n\tIDV Error = %f\t cycle error %f\n', e, cycle_error+0);


% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% Optimization function - called by fmincon
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
function [e d] = opt_func(optimization_parameters)
global g_EMU g_idv g_WC_known_metabolites_idv g_idv_known_arr g_iteration_num g_model_net_fluxes g_model g_model_thermodynamics g_idv_known_mat g_WC_known_metabolites_concentration g_num_of_experiments;
fluxes_net          = optimization_parameters(1:g_model_net_fluxes.rxn_num);
concentrations  = optimization_parameters(g_model_net_fluxes.rxn_num+1:g_model_net_fluxes.rxn_num+length(g_model_thermodynamics.mets));
flux_fb         = optimization_parameters(g_model_net_fluxes.rxn_num+length(g_model_thermodynamics.mets)+1:end);
% cy_mt_ratio     = optimizatio n_parameters(g_model_net_fluxes.rxn_num+length(g_model_thermodynamics.mets)+1:end);
g_iteration_num = g_iteration_num+1;

one_direction_flux = ComputeOneDirectionFluxes(fluxes_net, concentrations, flux_fb, g_model_net_fluxes, g_model_thermodynamics);

idv = cell(g_num_of_experiments,1);
idv_d   = cell(g_num_of_experiments,1);
fluxes_net_after_force_zero         = cell(g_num_of_experiments,1);
one_direction_flux_after_force_zero = cell(g_num_of_experiments,1);
for(i=1:g_num_of_experiments)
    fluxes_net_after_force_zero{i} = fluxes_net;
    fluxes_net_after_force_zero{i}(g_model_net_fluxes.force_zero_flux{i}==1)=g_model_net_fluxes.positive_direction_lb(g_model_net_fluxes.force_zero_flux{i}==1);
    one_direction_flux_after_force_zero{i} = ComputeOneDirectionFluxes(fluxes_net_after_force_zero{i}, concentrations, flux_fb, g_model_net_fluxes, g_model_thermodynamics);    
    gx_matrix{i} = ComputeGxDeriviatives(fluxes_net_after_force_zero{i}, concentrations, flux_fb, g_model_net_fluxes, g_model_thermodynamics, g_model, g_model_net_fluxes.force_zero_flux{i});    

    fcn_name = ['ComputeEmuIDV_Opt' int2str(i)];    
    [idv{i} idv_d{i}] = feval(fcn_name,g_idv{i},one_direction_flux_after_force_zero{i});
    
end


[e e_match_labeling e_match_concentrations] = ComputeError(g_WC_known_metabolites_idv, idv, g_idv_known_mat, g_WC_known_metabolites_concentration, g_EMU, concentrations, g_num_of_experiments);
d = ComputeErrorDeriv(g_WC_known_metabolites_idv, idv, idv_d, g_idv_known_mat, g_WC_known_metabolites_concentration, g_EMU, fluxes_net, one_direction_flux, gx_matrix, concentrations, g_num_of_experiments);


 %e = e+cycle_error
 
%  fprintf('%f\n', e);


% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% Compute forward and backward fluxes
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
function one_direction_flux = ComputeOneDirectionFluxes(fluxes_net, concentrations, fluxes_fb, model_net_fluxes, model_thermodynamics)
global g_constants g_flux_directionalities;
one_direction_flux = [];
one_direction_flux_ind = 1;
for(i=1:length(model_net_fluxes.is_net_flux))
    if(model_net_fluxes.is_net_flux(i))
        if((~isnan(g_flux_directionalities(i)))&&(model_thermodynamics.thermodynamics_of_reaction_defined(i)))
            numerator   = sum(concentrations(model_thermodynamics.product_indexes{i}));
            denominator = sum(concentrations(model_thermodynamics.reactant_indexes{i}));        

            fb_ratio = [((exp(numerator)/exp(denominator))^-1)*exp(-model_thermodynamics.delta_G0(i)/g_constants.RT)];
            if(fb_ratio==1)
                fb_ratio = fb_ratio+g_constants.EPSILON_VALUE;
            end
                % forward flux
                one_direction_flux(one_direction_flux_ind,1)    = fluxes_net(i)*fb_ratio/(fb_ratio-1);
                % backward flux
                one_direction_flux(one_direction_flux_ind+1,1)  = fluxes_net(i)/(fb_ratio-1);                
            if((one_direction_flux(one_direction_flux_ind)==Inf)||(one_direction_flux(one_direction_flux_ind+1)==Inf)||(one_direction_flux(one_direction_flux_ind)==-Inf)||(one_direction_flux(one_direction_flux_ind+1)==-Inf)||(isnan(one_direction_flux(one_direction_flux_ind)))||(isnan(one_direction_flux(one_direction_flux_ind+1))))
                alon=2
                fb_ratio
            end
        else
            one_direction_flux(one_direction_flux_ind,1)    = fluxes_fb(one_direction_flux_ind);
            % backward flux
            one_direction_flux(one_direction_flux_ind+1,1)  = fluxes_fb(one_direction_flux_ind+1);
        end        
        one_direction_flux_ind = one_direction_flux_ind+2;
    else
        one_direction_flux(one_direction_flux_ind,1) = fluxes_net(i);
        one_direction_flux_ind = one_direction_flux_ind+1;
    end
end     



% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% Compute gx matrix to transform gradient based forward backward flux to 
% gradient based net flux and concentrations
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
function gx_matrix = ComputeGxDeriviatives(fluxes_net, concentrations, flux_fb, model_net_fluxes, model_thermodynamics, model, force_flux_to_zero)
gx_matrix = zeros(model.rxn_num,length(fluxes_net)+length(concentrations)+length(flux_fb));
one_direction_flux_ind = 1;
global g_constants g_flux_directionalities;
for(i=1:length(fluxes_net))
    if(model_net_fluxes.is_net_flux(i))
        if((~isnan(g_flux_directionalities(i)))&&(model_thermodynamics.thermodynamics_of_reaction_defined(i)))
            numerator   = sum(concentrations(model_thermodynamics.product_indexes{i}));
            denominator = sum(concentrations(model_thermodynamics.reactant_indexes{i}));        

            fb_ratio = [((exp(numerator)/exp(denominator))^-1)*exp(-model_thermodynamics.delta_G0(i)/g_constants.RT)];
            % forward flux
            gx_matrix(one_direction_flux_ind, i) = fb_ratio/(fb_ratio-1);
            % backward flux
            gx_matrix(one_direction_flux_ind+1, i) = 1/(fb_ratio-1);

            for(j=1:length(model_thermodynamics.product_indexes{i}))
                gx_matrix(one_direction_flux_ind, length(fluxes_net)+model_thermodynamics.product_indexes{i}(j)) = fluxes_net(i)*((fb_ratio-1)^-2)*fb_ratio;
                gx_matrix(one_direction_flux_ind+1, length(fluxes_net)+model_thermodynamics.product_indexes{i}(j)) = fluxes_net(i)*((fb_ratio-1)^-2)*fb_ratio;
            end 
            for(j=1:length(model_thermodynamics.reactant_indexes{i}))
                gx_matrix(one_direction_flux_ind, length(fluxes_net)+model_thermodynamics.reactant_indexes{i}(j)) = -fluxes_net(i)*((fb_ratio-1)^-2)*fb_ratio;
                gx_matrix(one_direction_flux_ind+1, length(fluxes_net)+model_thermodynamics.reactant_indexes{i}(j)) = -fluxes_net(i)*((fb_ratio-1)^-2)*fb_ratio;
            end         
        else
            gx_matrix(one_direction_flux_ind, length(fluxes_net)+length(model_thermodynamics.mets)+one_direction_flux_ind) = 1;
            gx_matrix(one_direction_flux_ind+1, length(fluxes_net)+length(model_thermodynamics.mets)+one_direction_flux_ind+1) = 1;
        end
        one_direction_flux_ind = one_direction_flux_ind+2;
    else
        if((force_flux_to_zero(i)==1))
            gx_matrix(one_direction_flux_ind, i) = 0; %make derivative zero - we do not want it to try to change 
        else
            gx_matrix(one_direction_flux_ind, i) = 1;
        end
        
        one_direction_flux_ind = one_direction_flux_ind+1;
        
    end    
end    


% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% Compute error of of the fit of convoluted simulated fluxes and
% concentrations into measurements
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
function [e e_match_labeling e_match_concentrations]= ComputeError(g_WC_known_metabolites_idv, idv, g_idv_known_mat, g_WC_known_metabolites_concentration, g_EMU, concentrations, num_of_experiments)
global g_constants;
e = 0;
e_match_labeling = cell(num_of_experiments,1);
e_match_concentrations = [];

% Match IDVs
for(i=1:num_of_experiments)
    e_match_labeling{i} = [];
    for j=1:size(g_idv_known_mat{i},1)    
        x_cy = g_idv_known_mat{i}(j,1);
        x_mt = g_idv_known_mat{i}(j,2);

        EMU_indices_cy = find(g_EMU{i}.list(:,1)==x_cy);
        EMU_indices_mt = find(g_EMU{i}.list(:,1)==x_mt);

        if(~isnan(x_cy))
            idv_cy = idv{i}{EMU_indices_cy(1)};
        else
            idv_cy = idv{i}{EMU_indices_mt(1)};
        end
        if(~isnan(x_mt))
            idv_mt = idv{i}{EMU_indices_mt(1)};
        else
            idv_mt = idv{i}{EMU_indices_cy(1)};
        end
        

        x_cy_for_concentrations = g_WC_known_metabolites_idv{i}{j}.index_CY_for_concentrations;
        x_mt_for_concentrations = g_WC_known_metabolites_idv{i}{j}.index_MT_for_concentrations;


        iterator_cy_mt_ratio=g_constants.CY_WC_VOLUME*exp(concentrations(x_cy_for_concentrations))/(g_constants.CY_WC_VOLUME*exp(concentrations(x_cy_for_concentrations))+g_constants.MT_WC_VOLUME*exp(concentrations(x_mt_for_concentrations)));
        % compare known idv to the first EMU idv of this metabolite
        % it must be the first one from all EMU of the same metabolite, as the
        % first one contains all the carbons
        g_WC_known_metabolites_idv{i}{j}.idv_variance(g_WC_known_metabolites_idv{i}{j}.idv_variance<0.0001)=0.0001;
        try
        % if m+0 is contaminated, do normalized based to m+0 and compare
        % normalized measured labeling to normalized computational labeling            
        if(g_WC_known_metabolites_idv{i}{j}.contaminated_m0)
            fixed_idv=g_WC_known_metabolites_idv{i}{j}.idv(2:end)/(1-g_WC_known_metabolites_idv{i}{j}.idv(1));
            fixed_idv_var = g_WC_known_metabolites_idv{i}{j}.idv_variance(2:end)/(1-g_WC_known_metabolites_idv{i}{j}.idv(1));
            fixed_idv_simulated_cy = idv_cy(2:end)/(1-idv_cy(1));
            fixed_idv_simulated_mt = idv_mt(2:end)/(1-idv_mt(1));
            e_match_labeling{i}(end+1,1) = sum(((fixed_idv(2:end)-((iterator_cy_mt_ratio*fixed_idv_simulated_cy(2:end))+((1-iterator_cy_mt_ratio)*fixed_idv_simulated_mt(2:end))) ).^2)./fixed_idv_var(2:end));
        else
            e_match_labeling{i}(end+1,1) = sum(((g_WC_known_metabolites_idv{i}{j}.idv(2:end)-((iterator_cy_mt_ratio*idv_cy(2:end))+((1-iterator_cy_mt_ratio)*idv_mt(2:end))) ).^2)./g_WC_known_metabolites_idv{i}{j}.idv_variance(2:end));
        end
        
        e = e + e_match_labeling{i}(end);
        catch
            error('Error in function: ComputeError');
        end
    end
end
% e=0;

% Match concentrations
for i=1:length(g_WC_known_metabolites_concentration)        
    x_cy = g_WC_known_metabolites_concentration{i}.index_CY;
    x_mt = g_WC_known_metabolites_concentration{i}.index_MT;
    convoluted_calculated_concentration = ((g_constants.CY_WC_VOLUME*exp(concentrations(x_cy)))+(g_constants.MT_WC_VOLUME*exp(concentrations(x_mt))));
    % STD is based on the smaller one out of the measured concentrations or
    % the calculated convoluted concentrations
    if(g_WC_known_metabolites_concentration{i}.concentration <= convoluted_calculated_concentration)
        concentration_std = g_WC_known_metabolites_concentration{i}.concentration*g_constants.WC_CONCENTRATION_STD_FACTOR;
    else
        concentration_std = convoluted_calculated_concentration*g_constants.WC_CONCENTRATION_STD_FACTOR;
        if(concentration_std < 0.001)
            concentration_std = 0.001;
        end
    end
    e_match_concentrations(end+1,1) =  ((g_WC_known_metabolites_concentration{i}.concentration-convoluted_calculated_concentration)/concentration_std)^2;    
    e = e + e_match_concentrations(end);
end

    
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% Compute first order derivatives of the objective function (ObjF) with respect to each of the the optimization parameters
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
function d_fluxes_and_concentrations = ComputeErrorDeriv(g_WC_known_metabolites_idv, idv, idv_d, g_idv_known_mat, g_WC_known_metabolites_concentration, g_EMU, fluxes_net, one_direction_flux, gx_matrix, concentrations, num_of_experiments)
global g_constants;
rxn_num     = length(one_direction_flux);
rxn_net_num = length(fluxes_net);

% Compute derivatives of matched IDVs
d_fluxes_and_concentrations = 0;
for(i=1:num_of_experiments)
    d_one_way_fluxes = 0;
    d_cy_mt_ratios = [];
    iterator_cy_mt_ratio_d = zeros(length(gx_matrix{i}),size(g_idv_known_mat{i},1));
    for j=1:size(g_idv_known_mat{i},1)    
        x_cy_for_concentrations = g_WC_known_metabolites_idv{i}{j}.index_CY_for_concentrations;
        x_mt_for_concentrations = g_WC_known_metabolites_idv{i}{j}.index_MT_for_concentrations;

        iterator_cy_mt_ratio_d(rxn_net_num+x_cy_for_concentrations,j) = g_constants.CY_WC_VOLUME*exp(concentrations(x_cy_for_concentrations))/(g_constants.CY_WC_VOLUME*exp(concentrations(x_cy_for_concentrations))+g_constants.MT_WC_VOLUME*exp(concentrations(x_mt_for_concentrations))) - g_constants.CY_WC_VOLUME*g_constants.CY_WC_VOLUME*exp(concentrations(x_cy_for_concentrations))^2/((g_constants.CY_WC_VOLUME*exp(concentrations(x_cy_for_concentrations))+g_constants.MT_WC_VOLUME*exp(concentrations(x_mt_for_concentrations)))^2);
        iterator_cy_mt_ratio_d(rxn_net_num+x_mt_for_concentrations,j) = -g_constants.CY_WC_VOLUME*exp(concentrations(x_cy_for_concentrations))*g_constants.MT_WC_VOLUME*exp(concentrations(x_mt_for_concentrations))/((g_constants.CY_WC_VOLUME*exp(concentrations(x_cy_for_concentrations))+g_constants.MT_WC_VOLUME*exp(concentrations(x_mt_for_concentrations)))^2);    
    end

    d_fluxes_and_concentrations2=0;
    for j=1:size(g_idv_known_mat{i},1)    
        x_cy = g_idv_known_mat{i}(j,1);
        x_mt = g_idv_known_mat{i}(j,2);

        x_cy_for_concentrations = g_WC_known_metabolites_idv{i}{j}.index_CY_for_concentrations;
        x_mt_for_concentrations = g_WC_known_metabolites_idv{i}{j}.index_MT_for_concentrations;

        EMU_indices_cy = find(g_EMU{i}.list(:,1)==x_cy);
        EMU_indices_mt = find(g_EMU{i}.list(:,1)==x_mt);

        if(~isnan(x_cy))
            idv_cy = idv{i}{EMU_indices_cy(1)};
            idv_cy_d = idv_d{i}{EMU_indices_cy(1)};
        else
            idv_cy = idv{i}{EMU_indices_mt(1)};
            idv_cy_d = idv_d{i}{EMU_indices_mt(1)};
        end
        if(~isnan(x_mt))
            idv_mt = idv{i}{EMU_indices_mt(1)};
            idv_mt_d = idv_d{i}{EMU_indices_mt(1)};
        else
            idv_mt = idv{i}{EMU_indices_cy(1)};
            idv_mt_d = idv_d{i}{EMU_indices_cy(1)};
        end



        iterator_cy_mt_ratio=g_constants.CY_WC_VOLUME*exp(concentrations(x_cy_for_concentrations))/(g_constants.CY_WC_VOLUME*exp(concentrations(x_cy_for_concentrations))+g_constants.MT_WC_VOLUME*exp(concentrations(x_mt_for_concentrations)));

        % compare known idv to the first EMU idv of this metabolite
        % it must be the first one from all EMU of the same metabolite, as the
        % first one contains all the carbons
        g_WC_known_metabolites_idv{i}{j}.idv_variance(g_WC_known_metabolites_idv{i}{j}.idv_variance<0.0001)=0.0001;


        % compare known idv to the first EMU idv of this metabolite
        % it must be the first one from all EMU of the same metabolite, as the
        % first one contains all the carbons    
        
        % if m+0 is contaminated, do normalized based to m+0 and compare
        % normalized measured labeling to normalized computational labeling
        if(g_WC_known_metabolites_idv{i}{j}.contaminated_m0)
            fixed_idv=g_WC_known_metabolites_idv{i}{j}.idv(2:end)/(1-g_WC_known_metabolites_idv{i}{j}.idv(1));
            fixed_idv_var = g_WC_known_metabolites_idv{i}{j}.idv_variance(2:end)/(1-g_WC_known_metabolites_idv{i}{j}.idv(1));
            fixed_idv_simulated_cy = idv_cy(2:end)/(1-idv_cy(1));
            fixed_idv_simulated_mt = idv_mt(2:end)/(1-idv_mt(1));
            fixed_idv_simulated_cy_d = (idv_cy_d(:,2:end)*(1-idv_cy(1))-(-idv_cy_d(:,1))*idv_cy(2:end))/(1-idv_cy(1))^2;
            fixed_idv_simulated_mt_d = (idv_mt_d(:,2:end)*(1-idv_mt(1))-(-idv_mt_d(:,1))*idv_mt(2:end))/(1-idv_mt(1))^2;
            d_one_way_fluxes = d_one_way_fluxes + 2*sum( repmat((((iterator_cy_mt_ratio*fixed_idv_simulated_cy(2:end))+((1-iterator_cy_mt_ratio)*fixed_idv_simulated_mt(2:end)))-fixed_idv(2:end))./fixed_idv_var(2:end), rxn_num,1).*((iterator_cy_mt_ratio*fixed_idv_simulated_cy_d(:,2:end))+((1-iterator_cy_mt_ratio)*fixed_idv_simulated_mt_d(:,2:end))),2 );

            d_fluxes_and_concentrations1 = d_one_way_fluxes'*gx_matrix{i};
            d_fluxes_and_concentrations2 = d_fluxes_and_concentrations2 + 2*sum( repmat((((iterator_cy_mt_ratio*fixed_idv_simulated_cy(2:end))+((1-iterator_cy_mt_ratio)*fixed_idv_simulated_mt(2:end)))-fixed_idv(2:end))./fixed_idv_var(2:end), length(gx_matrix{i}),1).*((iterator_cy_mt_ratio_d(:,j)*fixed_idv_simulated_cy(2:end))+((-iterator_cy_mt_ratio_d(:,j))*fixed_idv_simulated_mt(2:end))),2 );            
        else
            d_one_way_fluxes = d_one_way_fluxes + 2*sum( repmat((((iterator_cy_mt_ratio*idv_cy(2:end))+((1-iterator_cy_mt_ratio)*idv_mt(2:end)))-g_WC_known_metabolites_idv{i}{j}.idv(2:end))./g_WC_known_metabolites_idv{i}{j}.idv_variance(2:end), rxn_num,1).*((iterator_cy_mt_ratio*idv_cy_d(:,2:end))+((1-iterator_cy_mt_ratio)*idv_mt_d(:,2:end))),2 );
            d_fluxes_and_concentrations1 = d_one_way_fluxes'*gx_matrix{i};
            d_fluxes_and_concentrations2 = d_fluxes_and_concentrations2 + 2*sum( repmat((((iterator_cy_mt_ratio*idv_cy(2:end))+((1-iterator_cy_mt_ratio)*idv_mt(2:end)))-g_WC_known_metabolites_idv{i}{j}.idv(2:end))./g_WC_known_metabolites_idv{i}{j}.idv_variance(2:end), length(gx_matrix{i}),1).*((iterator_cy_mt_ratio_d(:,j)*idv_cy(2:end))+((-iterator_cy_mt_ratio_d(:,j))*idv_mt(2:end))),2 );            
        end
    end
    d_fluxes_and_concentrations = d_fluxes_and_concentrations + d_fluxes_and_concentrations1' + d_fluxes_and_concentrations2;
end

% Compute derivatives of matched concentrations
for i=1:length(g_WC_known_metabolites_concentration)        
    x_cy = g_WC_known_metabolites_concentration{i}.index_CY;
    x_mt = g_WC_known_metabolites_concentration{i}.index_MT;
    
    convoluted_calculated_concentration = ((g_constants.CY_WC_VOLUME*exp(concentrations(x_cy)))+(g_constants.MT_WC_VOLUME*exp(concentrations(x_mt))));

    %derivatives based on fixed concentration STD 
    if(g_WC_known_metabolites_concentration{i}.concentration <= convoluted_calculated_concentration)
        concentration_std = g_WC_known_metabolites_concentration{i}.concentration*g_constants.WC_CONCENTRATION_STD_FACTOR;
        d_fluxes_and_concentrations(rxn_net_num+x_cy) = d_fluxes_and_concentrations(rxn_net_num+x_cy) + 2*((convoluted_calculated_concentration-g_WC_known_metabolites_concentration{i}.concentration)/concentration_std) * g_constants.CY_WC_VOLUME*exp(concentrations(x_cy))/concentration_std;
        d_fluxes_and_concentrations(rxn_net_num+x_mt) = d_fluxes_and_concentrations(rxn_net_num+x_mt) + 2*((convoluted_calculated_concentration-g_WC_known_metabolites_concentration{i}.concentration)/concentration_std) * g_constants.MT_WC_VOLUME*exp(concentrations(x_mt))/concentration_std;        
    else
        %derivatives based on calculated convoluted concentrations - in
        %this case the derivatives are based on the calculated
        %concentration => harder to calculate derivatives
        concentration_std = convoluted_calculated_concentration*g_constants.WC_CONCENTRATION_STD_FACTOR;
        %if the calculated STD is too small - in this case the derivatives
        %will be based again on fixed values - easier to calculated derivatives
        if(concentration_std < 0.001)
            concentration_std = 0.001;
            d_fluxes_and_concentrations(rxn_net_num+x_cy) = d_fluxes_and_concentrations(rxn_net_num+x_cy) + 2*((convoluted_calculated_concentration-g_WC_known_metabolites_concentration{i}.concentration)/concentration_std) * g_constants.CY_WC_VOLUME*exp(concentrations(x_cy))/concentration_std;
            d_fluxes_and_concentrations(rxn_net_num+x_mt) = d_fluxes_and_concentrations(rxn_net_num+x_mt) + 2*((convoluted_calculated_concentration-g_WC_known_metabolites_concentration{i}.concentration)/concentration_std) * g_constants.MT_WC_VOLUME*exp(concentrations(x_mt))/concentration_std;                    
        else
            numerator = (convoluted_calculated_concentration-g_WC_known_metabolites_concentration{i}.concentration);
            numerator_cy_d = g_constants.CY_WC_VOLUME*exp(concentrations(x_cy));
            numerator_mt_d = g_constants.MT_WC_VOLUME*exp(concentrations(x_mt));
            denumerator = concentration_std;
            denumerator_cy_d = g_constants.CY_WC_VOLUME*exp(concentrations(x_cy))*g_constants.WC_CONCENTRATION_STD_FACTOR;
            denumerator_mt_d = g_constants.MT_WC_VOLUME*exp(concentrations(x_mt))*g_constants.WC_CONCENTRATION_STD_FACTOR;
            division_cy_d = (numerator_cy_d*denumerator-denumerator_cy_d*numerator)/(denumerator^2);
            division_mt_d = (numerator_mt_d*denumerator-denumerator_mt_d*numerator)/(denumerator^2);

            d_fluxes_and_concentrations(rxn_net_num+x_cy) = d_fluxes_and_concentrations(rxn_net_num+x_cy) + 2*((convoluted_calculated_concentration-g_WC_known_metabolites_concentration{i}.concentration)/concentration_std) * division_cy_d;
            d_fluxes_and_concentrations(rxn_net_num+x_mt) = d_fluxes_and_concentrations(rxn_net_num+x_mt) + 2*((convoluted_calculated_concentration-g_WC_known_metabolites_concentration{i}.concentration)/concentration_std) * division_mt_d;               
        end        
    end
end
    
