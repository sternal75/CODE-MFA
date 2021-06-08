% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% Compute error of of the fit of convoluted simulated fluxes and
% concentrations into measurements
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
function [e e_match_labeling e_match_concentrations]= ComputeError(g_WC_known_metabolites_idv, idv, g_idv_known_mat, g_WC_known_metabolites_concentration, g_EMU, concentrations, num_of_experiments)

load_constants

g_STD_FACTOR = 0.2;
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


        iterator_cy_mt_ratio=CY_WC_VOLUME*exp(concentrations(x_cy_for_concentrations))/(CY_WC_VOLUME*exp(concentrations(x_cy_for_concentrations))+MT_WC_VOLUME*exp(concentrations(x_mt_for_concentrations)));
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
    convoluted_calculated_concentration = ((CY_WC_VOLUME*exp(concentrations(x_cy)))+(MT_WC_VOLUME*exp(concentrations(x_mt))));
    % STD is based on the smaller one out of the measured concentrations or
    % the calculated convoluted concentrations
    if(g_WC_known_metabolites_concentration{i}.concentration <= convoluted_calculated_concentration)
        concentration_std = g_WC_known_metabolites_concentration{i}.concentration*g_STD_FACTOR;
    else
        concentration_std = convoluted_calculated_concentration*g_STD_FACTOR;
        if(concentration_std < 0.001)
            concentration_std = 0.001;
        end
    end
    e_match_concentrations(end+1,1) =  ((g_WC_known_metabolites_concentration{i}.concentration-convoluted_calculated_concentration)/concentration_std)^2;    
    e = e + e_match_concentrations(end);
end

