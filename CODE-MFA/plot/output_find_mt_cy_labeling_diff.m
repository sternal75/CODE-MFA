% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% Find differences in labeling patterns in mitochondria and cytososl
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
function output_find_mt_cy_labeling_diff()
    load('../mat_files/sensitiviy_analysis_concentration.mat', 'sensitiviy_analysis_concentration');
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
    load('../mat_files/net_fluxes.mat', 'net_fluxes');
    load('../mat_files/sensitiviy_analysis_dG.mat', 'sensitiviy_analysis_dG');    

    addpath('../functions/emu') 
    addpath('../functions/general') 
    addpath('../')     
    run ../load_constants;


    % per concentrations lower/upper bounds
    % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
    for(i=1:length(sensitiviy_analysis_concentration))
        % lower bound
        predicted_fluxes_fb       = sensitiviy_analysis_concentration{i}.fb_fluxes.low;
        predicted_net_fluxes      = sensitiviy_analysis_concentration{i}.net_fluxes.low;
        predicted_concentrations  = sensitiviy_analysis_concentration{i}.concentrations.low;
        [labeling_diff_concentration_low_gln(i,:) labeling_diff_concentration_low_glc(i,:)] = find_mt_cy_labeling_diff(model, WC_known_metabolites_idv, WC_known_metabolites_concentration, model_net_fluxes, model_thermodynamics, EMU, idv, predicted_net_fluxes, predicted_fluxes_fb, predicted_concentrations);
    end
    for(i=1:length(sensitiviy_analysis_concentration))
        % lower bound
        predicted_fluxes_fb       = sensitiviy_analysis_concentration{i}.fb_fluxes.high;
        predicted_net_fluxes      = sensitiviy_analysis_concentration{i}.net_fluxes.high;
        predicted_concentrations  = sensitiviy_analysis_concentration{i}.concentrations.high;
        [labeling_diff_concentration_high_gln(i,:) labeling_diff_concentration_high_glc(i,:)] = find_mt_cy_labeling_diff(model, WC_known_metabolites_idv, WC_known_metabolites_concentration, model_net_fluxes, model_thermodynamics, EMU, idv, predicted_net_fluxes, predicted_fluxes_fb, predicted_concentrations);
    end    
    % per net_fluxes lower/upper bounds
    % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
    for(i=1:length(net_fluxes))
        % lower bound
        predicted_fluxes_fb       = net_fluxes{i}.fb_fluxes.low;
        predicted_net_fluxes      = net_fluxes{i}.net_fluxes.low;
        predicted_concentrations  = net_fluxes{i}.concentrations.low;
        [labeling_diff_net_flux_low_gln(i,:) labeling_diff_net_flux_low_glc(i,:)] = find_mt_cy_labeling_diff(model, WC_known_metabolites_idv, WC_known_metabolites_concentration, model_net_fluxes, model_thermodynamics, EMU, idv, predicted_net_fluxes, predicted_fluxes_fb, predicted_concentrations);
    end
    for(i=1:length(net_fluxes))
        % lower bound
        predicted_fluxes_fb       = net_fluxes{i}.fb_fluxes.high;
        predicted_net_fluxes      = net_fluxes{i}.net_fluxes.high;
        predicted_concentrations  = net_fluxes{i}.concentrations.high;
        [labeling_diff_net_flux_high_gln(i,:) labeling_diff_net_flux_high_glc(i,:)] = find_mt_cy_labeling_diff(model, WC_known_metabolites_idv, WC_known_metabolites_concentration, model_net_fluxes, model_thermodynamics, EMU, idv, predicted_net_fluxes, predicted_fluxes_fb, predicted_concentrations);
    end    
    % per dG lower/upper bounds
    % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
    for(i=1:length(sensitiviy_analysis_dG))
        % lower bound
        predicted_fluxes_fb       = sensitiviy_analysis_dG{i}.fb_fluxes.low;
        predicted_net_fluxes      = sensitiviy_analysis_dG{i}.net_fluxes.low;
        predicted_concentrations  = sensitiviy_analysis_dG{i}.concentrations.low;
        [labeling_diff_dG_low_gln(i,:) labeling_diff_dG_low_glc(i,:)] = find_mt_cy_labeling_diff(model, WC_known_metabolites_idv, WC_known_metabolites_concentration, model_net_fluxes, model_thermodynamics, EMU, idv, predicted_net_fluxes, predicted_fluxes_fb, predicted_concentrations);
    end
    for(i=1:length(sensitiviy_analysis_dG))
        % lower bound
        predicted_fluxes_fb       = sensitiviy_analysis_dG{i}.fb_fluxes.high;
        predicted_net_fluxes      = sensitiviy_analysis_dG{i}.net_fluxes.high;
        predicted_concentrations  = sensitiviy_analysis_dG{i}.concentrations.high;
        [labeling_diff_dG_high_gln(i,:) labeling_diff_dG_high_glc(i,:)] = find_mt_cy_labeling_diff(model, WC_known_metabolites_idv, WC_known_metabolites_concentration, model_net_fluxes, model_thermodynamics, EMU, idv, predicted_net_fluxes, predicted_fluxes_fb, predicted_concentrations);
    end    
    
    sensitivity_analysis_label={'low concentration','high concentration','low flux', 'high flux', 'low dG', 'high dG'};
    
    
    [max_concentration_low_diff  max_concentration_low_location]    = max(labeling_diff_concentration_low_gln);
    [max_concentration_high_diff max_concentration_high_location]   = max(labeling_diff_concentration_high_gln);    
    [max_net_flux_low_diff  max_net_flux_low_location]              = max(labeling_diff_net_flux_low_gln);
    [max_net_flux_high_diff max_net_flux_high_location]             = max(labeling_diff_net_flux_high_gln);        
    [max_dG_low_diff  max_dG_low_location]                          = max(labeling_diff_dG_low_gln);
    [max_dG_high_diff max_dG_high_location]                         = max(labeling_diff_dG_high_gln);            
    
    max_concentration_sensitivity_analysis = [max_concentration_low_location;max_concentration_high_location;max_net_flux_low_location;max_net_flux_high_location;max_dG_low_location;max_dG_high_location];
    
    [max_diff max_location] = max([max_concentration_low_diff;max_concentration_high_diff;max_net_flux_low_diff;max_net_flux_high_diff;max_dG_low_diff;max_dG_high_diff]);
    fprintf('********** max labeling diff from gln **********\n')
    
    for(i=1:size(max_diff,2))
        if(max_diff(i)>1e-5)
            if((max_location(i)==1)||(max_location(i)==2))
                sensitivity_analysis_on_met_or_flux = model_net_fluxes.mets{max_concentration_sensitivity_analysis(max_location(i),i)};
            else
                sensitivity_analysis_on_met_or_flux = model_net_fluxes.rxns{max_concentration_sensitivity_analysis(max_location(i),i)};
            end
            
            fprintf('%s max diff is %.4f, found at %s %s\n',WC_known_metabolites_idv{1}{i}.met_name,max_diff(i), sensitivity_analysis_label{max_location(i)}, sensitivity_analysis_on_met_or_flux);
        end
    end
    fprintf('\n\n\n\n\n');
    [max_concentration_low_diff  max_concentration_low_location]    = max(labeling_diff_concentration_low_glc);
    [max_concentration_high_diff max_concentration_high_location]   = max(labeling_diff_concentration_high_glc);    
    [max_net_flux_low_diff  max_net_flux_low_location]              = max(labeling_diff_net_flux_low_glc);
    [max_net_flux_high_diff max_net_flux_high_location]             = max(labeling_diff_net_flux_high_glc);        
    [max_dG_low_diff  max_dG_low_location]                          = max(labeling_diff_dG_low_glc);
    [max_dG_high_diff max_dG_high_location]                         = max(labeling_diff_dG_high_glc);            
    
    max_concentration_sensitivity_analysis = [max_concentration_low_location;max_concentration_high_location;max_net_flux_low_location;max_net_flux_high_location;max_dG_low_location;max_dG_high_location];
    
    [max_diff max_location] = max([max_concentration_low_diff;max_concentration_high_diff;max_net_flux_low_diff;max_net_flux_high_diff;max_dG_low_diff;max_dG_high_diff]);
    fprintf('********** max labeling diff from glc **********\n')
    
    for(i=1:size(max_diff,2))
        if(max_diff(i)>1e-5)
            if((max_location(i)==1)||(max_location(i)==2))
                sensitivity_analysis_on_met_or_flux = model_net_fluxes.mets{max_concentration_sensitivity_analysis(max_location(i),i)};
            else
                sensitivity_analysis_on_met_or_flux = model_net_fluxes.rxns{max_concentration_sensitivity_analysis(max_location(i),i)};
            end
            
            fprintf('%s max diff is %.4f, found at %s %s\n',WC_known_metabolites_idv{2}{i}.met_name,max_diff(i), sensitivity_analysis_label{max_location(i)}, sensitivity_analysis_on_met_or_flux);
        end
    end
    

    function [labeling_diff_gln labeling_diff_glc] = find_mt_cy_labeling_diff(model, WC_known_metabolites_idv, WC_known_metabolites_concentration, model_net_fluxes, model_thermodynamics, EMU, idv, predicted_net_fluxes, predicted_fluxes_fb, predicted_concentrations);
    num_of_experiments = length(WC_known_metabolites_idv);
    idv_opt = cell(num_of_experiments,1);
    idv_d   = cell(num_of_experiments,1);
    idv_known_arr = cell(num_of_experiments,1);
    idv_known_mat = cell(num_of_experiments,1);
    for (i=1:num_of_experiments)
        idv_known_arr{i} = zeros(0,1);
        idv_known_mat{i}=[];
        for j=1:length(WC_known_metabolites_idv{i})
            idv_known_arr{i} = [idv_known_arr{i};WC_known_metabolites_idv{i}{j}.index_CY;WC_known_metabolites_idv{i}{j}.index_MT];
            idv_known_mat{i} = [idv_known_mat{i};WC_known_metabolites_idv{i}{j}.index_CY WC_known_metabolites_idv{i}{j}.index_MT];
        end    
    end

    fluxes_net_after_force_zero         = cell(num_of_experiments,1);
    one_direction_flux_after_force_zero = cell(num_of_experiments,1);
    for(i=1:num_of_experiments)
        fluxes_net_after_force_zero{i} = predicted_net_fluxes;
        fluxes_net_after_force_zero{i}(model_net_fluxes.force_zero_flux{i}==1)=0.1;
        one_direction_flux_after_force_zero{i} = ComputeOneDirectionFluxes(fluxes_net_after_force_zero{i}, predicted_concentrations, predicted_fluxes_fb, model_net_fluxes, model_thermodynamics);    
        fcn_name = ['ComputeEmuIDV_Opt' int2str(1)];
        [idv_opt{i} idv_d{i}] = feval(fcn_name,idv{i},one_direction_flux_after_force_zero{i});        
    end


    [e e_match_labeling e_match_concentrations] = ComputeError(WC_known_metabolites_idv, idv_opt, idv_known_mat, WC_known_metabolites_concentration, EMU, predicted_concentrations, num_of_experiments);
    labeling_diff_gln = cy_vs_mt_labaling_diff(WC_known_metabolites_idv{1}, WC_known_metabolites_concentration, EMU{1}, idv{1}, 1, one_direction_flux_after_force_zero{1}, predicted_concentrations, e_match_labeling{1}, e_match_concentrations, model, 'Gln');
    labeling_diff_glc = cy_vs_mt_labaling_diff(WC_known_metabolites_idv{2}, WC_known_metabolites_concentration, EMU{2}, idv{2}, 2, one_direction_flux_after_force_zero{2}, predicted_concentrations, e_match_labeling{2}, e_match_concentrations, model, 'Glc');


    
    function one_direction_flux = ComputeOneDirectionFluxes(fluxes_net, concentrations, fluxes_fb, model_net_fluxes, model_thermodynamics)
    flux_directionalities = fluxes_net > 0;
    load_constants
    one_direction_flux = [];
    one_direction_flux_ind = 1;
    for(i=1:length(model_net_fluxes.is_net_flux))
        if(model_net_fluxes.is_net_flux(i))
            if((~isnan(flux_directionalities(i)))&&(model_thermodynamics.thermodynamics_of_reaction_defined(i)))
                numerator   = sum(concentrations(model_thermodynamics.product_indexes{i}));
                denominator = sum(concentrations(model_thermodynamics.reactant_indexes{i}));        

                fb_ratio = [((exp(numerator)/exp(denominator))^-1)*exp(-model_thermodynamics.delta_G0(i)/RT)];
                    one_direction_flux(one_direction_flux_ind,1)    = fluxes_net(i)*fb_ratio/(fb_ratio-1);
                    % backward flux
                    one_direction_flux(one_direction_flux_ind+1,1)  = fluxes_net(i)/(fb_ratio-1);                
                if((one_direction_flux(one_direction_flux_ind)==Inf)||(one_direction_flux(one_direction_flux_ind+1)==Inf)||(one_direction_flux(one_direction_flux_ind)==-Inf)||(one_direction_flux(one_direction_flux_ind+1)==-Inf)||(isnan(one_direction_flux(one_direction_flux_ind)))||(isnan(one_direction_flux(one_direction_flux_ind+1))))
                    alon=1
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


function labeling_diff_result = cy_vs_mt_labaling_diff(WC_known_metabolites_idv, WC_known_metabolites_concentration, EMU, idv, exp_index, final_predicted_flux, final_predicted_concentrations, error_labeling, error_concentrations, model, labeling_str)

    load_constants
    
    idv_known_arr = zeros(0,1);
    idv_known_mat=[];
    for i=1:length(WC_known_metabolites_idv)
        idv_known_arr = [idv_known_arr;WC_known_metabolites_idv{i}.index_CY;WC_known_metabolites_idv{i}.index_MT];
        idv_known_mat = [idv_known_mat;WC_known_metabolites_idv{i}.index_CY WC_known_metabolites_idv{i}.index_MT];
    end

    fcn_name = ['ComputeEmuIDV_Opt' int2str(exp_index)];
    [idv_opt idv_d] = feval(fcn_name,idv,final_predicted_flux);

    for i=1:size(idv_known_mat,1)
        %figure('Position', get(0, 'Screensize'));
%         subplot(3,floor(((size(idv_known_mat,1)+2)/3)),i)
        idv_known_vs_opt=[];
        xLabelForBar=cell(0);
        mass_isotopomer_legend=cell(0);        
        x_CY = idv_known_mat(i,1);
        x_MT = idv_known_mat(i,2);
        EMU_indices_CY = find(EMU.list(:,1)==x_CY);
        EMU_indices_MT = find(EMU.list(:,1)==x_MT);
        for(j=1:length(WC_known_metabolites_idv{i}.idv))
            mass_isotopomer_legend{end+1}=sprintf('m+%s',num2str(j-1));
        end
    
        % calculate cy/mt ratio for metabolie
        x_cy_for_concentrations = WC_known_metabolites_idv{i}.index_CY_for_concentrations;
        x_mt_for_concentrations = WC_known_metabolites_idv{i}.index_MT_for_concentrations;
        cy_mt_ratio=CY_WC_VOLUME*exp(final_predicted_concentrations(x_cy_for_concentrations))/(CY_WC_VOLUME*exp(final_predicted_concentrations(x_cy_for_concentrations))+MT_WC_VOLUME*exp(final_predicted_concentrations(x_mt_for_concentrations)));
        
        idv_known_vs_opt(end+1,:)=WC_known_metabolites_idv{i}.idv;
        xLabelForBar{end+1}='WC (m)'; 
        if((~isnan(x_CY))&&(~isnan(x_MT)))
            idv_known_vs_opt(end+1,1:length(idv_opt{EMU_indices_CY(1)}))=cy_mt_ratio*idv_opt{EMU_indices_CY(1)}+(1-cy_mt_ratio)*idv_opt{EMU_indices_MT(1)};
            xLabelForBar{end+1}='WC (s)';                       
            idv_known_vs_opt(end+1,1:length(idv_opt{EMU_indices_CY(1)}))=idv_opt{EMU_indices_CY(1)};           
            xLabelForBar{end+1}='CY (s)';                       
            idv_known_vs_opt(end+1,1:length(idv_opt{EMU_indices_MT(1)}))=idv_opt{EMU_indices_MT(1)};
            xLabelForBar{end+1}='MT (s)';            
        elseif(~isnan(x_CY))
            idv_known_vs_opt(end+1,1:length(idv_opt{EMU_indices_CY(1)}))=idv_opt{EMU_indices_CY(1)};           
            xLabelForBar{end+1}='WC (s)';           
            idv_known_vs_opt(end+1,1:length(idv_opt{EMU_indices_CY(1)}))=idv_opt{EMU_indices_CY(1)};           
            xLabelForBar{end+1}='CY (s)';           
            idv_known_vs_opt(end+1,1)=0;
            xLabelForBar{end+1}='MT (s)';           
        elseif(~isnan(x_MT))
            idv_known_vs_opt(end+1,1:length(idv_opt{EMU_indices_MT(1)}))=idv_opt{EMU_indices_MT(1)};
            xLabelForBar{end+1}='WC (s)';           
            idv_known_vs_opt(end+1,1)=0;
            xLabelForBar{end+1}='CY (s)';           
            idv_known_vs_opt(end+1,1:length(idv_opt{EMU_indices_MT(1)}))=idv_opt{EMU_indices_MT(1)};
            xLabelForBar{end+1}='MT (s)';
        end
        
        
        x_cy_for_concentrations = WC_known_metabolites_idv{i}.index_CY_for_concentrations;
        x_mt_for_concentrations = WC_known_metabolites_idv{i}.index_MT_for_concentrations;
        cy_mt_ratio=CY_WC_VOLUME*exp(final_predicted_concentrations(x_cy_for_concentrations))/(CY_WC_VOLUME*exp(final_predicted_concentrations(x_cy_for_concentrations))+MT_WC_VOLUME*exp(final_predicted_concentrations(x_mt_for_concentrations)));
        
        % show sum of labeling diff. If the metabolite has only mt or cy,
        % output 0
        if((~isnan(x_CY))&&(~isnan(x_MT)))
            labeling_diff_result(i) = sum(abs(idv_known_vs_opt(3,:)-idv_known_vs_opt(4,:)));
        else
            labeling_diff_result(i) = 0;
        end
    end
    




    
    
    
        
    
