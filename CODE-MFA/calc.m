clear all;
ts=clock;
rand('seed',(ts(6)*10000));

addpath('./functions/emu') 
addpath('./functions/general') 
addpath(genpath('./cobratoolbox/src'));

% load model including carbon mapping
LoadModel

% load model including themodynamic input
LoadThermodynamicsModel

% input glutamine Mass Isotopomer Distribution (MID)
file_name = 'input_glutamine';
run processIsotopicLabel/calc;
met_list_norm_glutamine = met_list_norm;
clear met_list_norm; clear sample_list_short;
% input glucose mass isotopomer distribution
file_name = 'input_glucose';
run processIsotopicLabel/calc;
met_list_norm_glucose = met_list_norm;
clear met_list_norm; clear sample_list_short;

% Stoichiometric matix size 
[n,m] = size(model.S);

% convert atom mapping format of model to a one used by FindEMU
m = [];
for i=1:length(model.rxns)   
    model.mappings_carbon{i}.mapping_r(model.mappings_carbon{i}.mapping_p) = [1:length(model.mappings_carbon{i}.mapping_p)]';    
    model.mappings_carbon{i}.mapping_mat_p = CreateAtomMappingMat(model, model.mappings_carbon{i}, i);
    m.graph_p = model.mappings_carbon{i}.graph_r;
    m.graph_r = model.mappings_carbon{i}.graph_p;  
    m.mapping_p = model.mappings_carbon{i}.mapping_r;  
    m.mapping_r = model.mappings_carbon{i}.mapping_p;  
    model.mappings_carbon{i}.mapping_mat_r = CreateAtomMappingMat(model, m, i);
end


% Change net flux to foward backward flux - add backward flux
% Keep original net reactions model
model_net_fluxes = model;
model_net_fluxes.is_net_flux = zeros(model_net_fluxes.rxn_num,1); % Mark bidirectional fluxes with '1' and one directional fluxes with 0;
added_fluxes = 0;
for i=1:length(model.rxns)
    ind = i+added_fluxes;
    if(~isempty(findstr(model.rxns{ind},'f')))
        model_net_fluxes.is_net_flux(i) = 1;
        backward_flux_name = strrep(model.rxns{ind},'f','b');
        model.rxns = {model.rxns{1:ind} backward_flux_name model.rxns{ind+1:end}};
        model.positive_direction_lb = [model.positive_direction_lb(1:ind); model.positive_direction_lb(ind); model.positive_direction_lb(ind+1:end)];
        model.positive_direction_ub = [model.positive_direction_ub(1:ind); model.positive_direction_ub(ind); model.positive_direction_ub(ind+1:end)];
        model.S = [model.S(:,1:ind) -model.S(:,ind) model.S(:,ind+1:end)];
        
        model.mappings_carbon = {model.mappings_carbon{1:ind} model.mappings_carbon{ind} model.mappings_carbon{ind+1:end}};
        model.mappings_carbon{ind+1}.graph_r = model.mappings_carbon{ind}.graph_p;
        model.mappings_carbon{ind+1}.graph_p = model.mappings_carbon{ind}.graph_r; 
        model.mappings_carbon{ind+1}.mapping_r = model.mappings_carbon{ind}.mapping_p';
        model.mappings_carbon{ind+1}.mapping_p = model.mappings_carbon{ind}.mapping_r'; 
        model.mappings_carbon{ind+1}.mapping_mat_r = model.mappings_carbon{ind}.mapping_mat_p;
        model.mappings_carbon{ind+1}.mapping_mat_p = model.mappings_carbon{ind}.mapping_mat_r;
       
        model.equality_constraints = [model.equality_constraints(:,1:ind) zeros(size(model.equality_constraints(:,1))) model.equality_constraints(:,ind+1:end)];

        added_fluxes = added_fluxes+1;
    end
end
% Add rows to equality_constraints to handle backward fluxes equalities
equality_constraints_matrix_size = size(model.equality_constraints);
for i=1:equality_constraints_matrix_size(1,1)
    fb_reaction_index = find(model.equality_constraints(i,:)==1);
    % If there are more than one, take only the first (engough to check if
    % it is a bidirectional flux)
    fb_reaction_index = fb_reaction_index(1,1); 
    % Add another row to handle backward equalities
    if(~isempty(findstr(model.rxns{fb_reaction_index},'f')))
        model.equality_constraints = [model.equality_constraints;0 model.equality_constraints(i,1:end-1)];
    end
    
end

model.rxn_num = length(model.rxns); % Update number of reactions
model.used_reactions_status = [1:length(model.rxns)]'; 
model.exchange = zeros(model.rxn_num,1);

 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% model.mets        = name of metabolites (short names)
% model.metNames    = name of metabolites (full names)
% model.rxns        = name of reactions (short names)
% model.atom_C_num  = number of carbon atoms in each metabolite
% model.met_extra   = binar vector indicating whether metabolites are
%                     external to the model. External metabolites are not
%                     mass-balanced and their isotopic labeling is
%                     assumed to be known.
% model.used_reactions_status = ?? NOT USED??
% model.exchange              = ?? NOT USED??
% model.S                     = stoichiometric matrix
% model.met_num               = number of metabolites (rows in model.S)
% model.rxn_num               = number of reactions (columns in model.S)
% model.mapping_carbon        = atom mapping (cell array whose length is
%                               the number of reactions); see
%                               EMU_implementation_notes.ppt)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
% EMU_met_known - the EMUs for which experimental mass-isotopomer
% distribution data is exeprimentally available. We also refer to mass-isotopomer
% distribution by IDV (isotopomer distribution vector)
EMU_met_known_glutamine = cell(model.met_num,1);
EMU_met_known_glucose = cell(model.met_num,1);
for i=1:model.met_num
    EMU_met_known_glutamine{i} = zeros(model.atom_C_num(i), 0);
    EMU_met_known_glucose{i} = zeros(model.atom_C_num(i), 0);
end


% for WC metabolites with known IDV - glutamine labeling
WC_known_metabolites_idv_glutamine=cell(0);
for(i=1:length(met_list_norm_glutamine))
    met_name = met_list_norm_glutamine{i}.met_name;
    met_name_CY = strcat(met_name,'_CY');
    met_name_MT = strcat(met_name,'_MT');    
    % Cytosolic (CY) metabolites
    index_CY=strfind(model.mets,met_name_CY);
    index_CY=find(~cellfun(@isempty,index_CY));
    if(~isempty(index_CY))
        EMU_met_known_glutamine{index_CY}=ones(model.atom_C_num(index_CY), 1);
    end
    index_CY_for_concentrations=strfind(model_thermodynamics.mets,met_name_CY);    
    index_CY_for_concentrations=find(~cellfun(@isempty,index_CY_for_concentrations));            

    % Mitochondrial (MT) metabolites   
    index_MT=strfind(model.mets,met_name_MT);
    index_MT=find(~cellfun(@isempty,index_MT));
    % Added to handle Lactate which exists only in MT and not in CY. 
    % This was added just to handle it for non compartmenalized MID
    % measurements in a compartmenalized model
    if(~isempty(index_MT))
        EMU_met_known_glutamine{index_MT}=ones(model.atom_C_num(index_MT), 1);        
    end
    index_MT_for_concentrations=strfind(model_thermodynamics.mets,met_name_MT);    
    index_MT_for_concentrations=find(~cellfun(@isempty,index_MT_for_concentrations));        
    
   
    % keep all known idvs and known idv variances of WC metabolites
    [num_of_mass_isotopomers num_of_timepoints]=size(met_list_norm_glutamine{i}.data);
    WC_known_metabolites_idv_glutamine{end+1}.idv         = met_list_norm_glutamine{i}.data(:,num_of_timepoints)'; 
    WC_known_metabolites_idv_glutamine{end}.idv_variance  = met_list_norm_glutamine{i}.var(:,num_of_timepoints)';
    WC_known_metabolites_idv_glutamine{end}.contaminated_m0 = met_list_norm_glutamine{i}.contaminated_m0;
    WC_known_metabolites_idv_glutamine{end}.met_name        = met_list_norm_glutamine{i}.met_name;    
    WC_known_metabolites_idv_glutamine{end}.index_CY_for_concentrations      = index_CY_for_concentrations;    
    WC_known_metabolites_idv_glutamine{end}.index_MT_for_concentrations      = index_MT_for_concentrations;        
    if(isempty(index_CY))
        WC_known_metabolites_idv_glutamine{end}.index_CY       = nan;    
    else
        WC_known_metabolites_idv_glutamine{end}.index_CY       = index_CY;    
    end
    if(isempty(index_MT))
        WC_known_metabolites_idv_glutamine{end}.index_MT       = nan;    
    else
        WC_known_metabolites_idv_glutamine{end}.index_MT       = index_MT;    
    end    
end

% for WC metabolites with known IDV - glucose labeling
WC_known_metabolites_idv_glucose=cell(0);
for(i=1:length(met_list_norm_glucose))
    met_name = met_list_norm_glucose{i}.met_name;
    met_name_CY = strcat(met_name,'_CY');
    met_name_MT = strcat(met_name,'_MT');    
    % Cytosolic metabolites
    index_CY=strfind(model.mets,met_name_CY);
    index_CY=find(~cellfun(@isempty,index_CY));
    if(~isempty(index_CY))
        EMU_met_known_glucose{index_CY}=ones(model.atom_C_num(index_CY), 1);
    end
    index_CY_for_concentrations=strfind(model_thermodynamics.mets,met_name_CY);    
    index_CY_for_concentrations=find(~cellfun(@isempty,index_CY_for_concentrations));            

    % Mitochondrial metabolites   
    index_MT=strfind(model.mets,met_name_MT);
    index_MT=find(~cellfun(@isempty,index_MT));
    % Added to handle Lactate which exists only in MT and not in CY. 
    % This was added just to handle it for non compartmenalized MID
    % measurements in a compartmenalized model
    if(~isempty(index_MT))
        EMU_met_known_glucose{index_MT}=ones(model.atom_C_num(index_MT), 1);        
    end
    index_MT_for_concentrations=strfind(model_thermodynamics.mets,met_name_MT);    
    index_MT_for_concentrations=find(~cellfun(@isempty,index_MT_for_concentrations));        
      
    % Keep all known idvs and known idv variances of WC metabolites
    [num_of_mass_isotopomers num_of_timepoints]=size(met_list_norm_glucose{i}.data);
    WC_known_metabolites_idv_glucose{end+1}.idv         = met_list_norm_glucose{i}.data(:,num_of_timepoints)'; 
    WC_known_metabolites_idv_glucose{end}.idv_variance  = met_list_norm_glucose{i}.var(:,num_of_timepoints)';
    WC_known_metabolites_idv_glucose{end}.contaminated_m0   = met_list_norm_glucose{i}.contaminated_m0;
    WC_known_metabolites_idv_glucose{end}.met_name          = met_list_norm_glucose{i}.met_name;    
    WC_known_metabolites_idv_glucose{end}.index_CY_for_concentrations      = index_CY_for_concentrations;    
    WC_known_metabolites_idv_glucose{end}.index_MT_for_concentrations      = index_MT_for_concentrations;        
    if(isempty(index_CY))
        WC_known_metabolites_idv_glucose{end}.index_CY       = nan;    
    else
        WC_known_metabolites_idv_glucose{end}.index_CY       = index_CY;    
    end
    if(isempty(index_MT))
        WC_known_metabolites_idv_glucose{end}.index_MT       = nan;    
    else
        WC_known_metabolites_idv_glucose{end}.index_MT       = index_MT;    
    end    
end

% For WC metabolites with known concentrations
WC_known_metabolites_concentration=cell(0);
for(i=1:length(model_thermodynamics.WC.met_name))
    met_name = model_thermodynamics.WC.met_name{i};
    met_name_CY = strcat(met_name,'_CY');
    met_name_MT = strcat(met_name,'_MT');    
        
    % Cytosolic metabolites
    index_CY=strcmp(model_thermodynamics.mets,met_name_CY);
    index_CY=find(index_CY);
    % Mitochondrial metabolites    
    index_MT=strcmp(model_thermodynamics.mets,met_name_MT);
    index_MT=find(index_MT);
    
    WC_known_metabolites_concentration{end+1}.concentration         = model_thermodynamics.WC.Concentrations(i); 
    WC_known_metabolites_concentration{end}.concentrations_STD  = model_thermodynamics.WC.Concentrations_STD(i); 
    WC_known_metabolites_concentration{end}.met_name      = met_name;        
    WC_known_metabolites_concentration{end}.index_CY      = index_CY;    
    WC_known_metabolites_concentration{end}.index_MT      = index_MT;    
end

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% Find EMUs by going backwards from the starting EMUs
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
fprintf('Find EMU\n');
[EMU_list_glutamine EMU_met_glutamine EMU_reactions_glutamine]  = FindEMU (model, EMU_met_known_glutamine); 
[EMU_list_glucose   EMU_met_glucose   EMU_reactions_glucose]    = FindEMU (model, EMU_met_known_glucose); 

%EMU_list - the list of identified EMUs (20 EMUs in the example network)
%           Each row corresponds to one EMU. The value in the first column
%           is the metabolite in which the EMU is defined. The value in the
%           second column is a running index of the EMU in that metabolite.
%           In the example network there are 6 EMUs for metabolite C (i.e.
%           6 rows in EMU_list whose value in the first column in 3
%           (representing C).
%
% EMU_met - the specific carbons in each EMU. This is a cell array whose
%           length is the number of metabolites in the network. For
%           metabolite i, EMU_met{i} is a matrix whose columns represent
%           the EMUs of metabolite i and the rows represent the carbons.
%           metabolites (a value of 1 represent that the carbon is part of
%           the EMU (see EMU_implementation_notes.ppt).
%
% EMU_reactions - the reactions that produce each EMU. This is a matrix
%                 with 4 columns. Each row i repreesnts a reaction producing
%                 a certain EMU (referred to as an "EMU reaction")
%         EMU_reactions(i,4) - an index of a reaction in the network
%         EMU_reactions(i,3) - an index of an EMU produced by this reaction
%         EMU_reactions(i,1) - an index (row in EMU_list) of an EMU used as
%                              substrate 
%         EMU_reactions(i,2) - an index (row in EMU_list) of a 2nd EMU used as
%                              substrate (in case the combines two EMUs)
%                             
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Cluster the EMUs, create DAG of clusters, and prepare matrices for efficient computation of IDVs given flux rates
fprintf('Analyze EMU\n');
EMU_glutamine   = AnalyzeEMU(model, EMU_list_glutamine, EMU_met_glutamine, EMU_reactions_glutamine);
EMU_glucose     = AnalyzeEMU(model, EMU_list_glucose, EMU_met_glucose, EMU_reactions_glucose);

% handle extracellular metabolite labeling from gluatmine and glucose
extra_met_isotopomers_glutamine = cell(model.met_num,1);
for(i=1:model.met_num)
    extra_met_isotopomers_glutamine{i}  = model.met_extra_labeling_glutamine{i};
    extra_met_isotopomers_glucose{i}    = model.met_extra_labeling_glucose{i};
end

idv_glutamine   = CreateIDV(EMU_glutamine, extra_met_isotopomers_glutamine);
idv_glucose     = CreateIDV(EMU_glucose, extra_met_isotopomers_glucose);


WC_known_metabolites_idv{1} = WC_known_metabolites_idv_glutamine;
WC_known_metabolites_idv{2} = WC_known_metabolites_idv_glucose;
EMU{1} = EMU_glutamine;
EMU{2} = EMU_glucose;
idv{1}=idv_glutamine;
idv{2}=idv_glucose;

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% Save MAT files
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
save('mat_files/model.mat','model');
save('mat_files/model_net_fluxes.mat','model_net_fluxes');
save('mat_files/model_thermodynamics.mat','model_thermodynamics');
save('mat_files/EMU.mat','EMU');
save('mat_files/idv.mat','idv');
save('mat_files/WC_known_metabolites_idv.mat','WC_known_metabolites_idv');
save('mat_files/WC_known_metabolites_concentration.mat','WC_known_metabolites_concentration');

% Create the computeEmuIDV_opt function
num_of_experiments = length(WC_known_metabolites_idv);
idv_known_arr = cell(num_of_experiments,1);
idv_known_mat = cell(num_of_experiments,1);
for (i=1:num_of_experiments)
    idv_known_arr{i} = zeros(0,1);
    idv_known_mat{i}=[];
    for j=1:length(WC_known_metabolites_idv{i})
        if(isnan(WC_known_metabolites_idv{i}{j}.index_MT))
            alon=1;
        end
        idv_known_arr{i} = [idv_known_arr{i};WC_known_metabolites_idv{i}{j}.index_CY;WC_known_metabolites_idv{i}{j}.index_MT];
        idv_known_mat{i} = [idv_known_mat{i};WC_known_metabolites_idv{i}{j}.index_CY WC_known_metabolites_idv{i}{j}.index_MT];
    end    
end
for(i=1:num_of_experiments)
    fcn_name = ['ComputeEmuIDV_Opt' int2str(i)];
    ComputeEmuIDV(EMU{i}, idv{i}, idv_known_arr{i}, ones(length(model.rxns),1), fcn_name);
end

% % % % % % % % % % % % % % % % % % % % % % % % % % 
% Step I and II - see figure 2 in the paper
% % % % % % % % % % % % % % % % % % % % % % % % % % 
iterations_result = iterate_to_find_fixed_net_flux_directionalities(model, model_net_fluxes, model_thermodynamics, EMU, idv, WC_known_metabolites_idv, WC_known_metabolites_concentration)
save('mat_files/iterations_result.mat','iterations_result');

% % % % % % % % % % % % % % % % % % % % % % % % % % 
% Step III - see figure 2 in the paper
% % % % % % % % % % % % % % % % % % % % % % % % % % 
MILP_results = FindAllFluxDirectionalities(model_net_fluxes, model_thermodynamics, iterations_result.directionalities);
save('mat_files/MILP_results.mat','MILP_results');

% % % % % % % % % % % % % % % % % % % % % % % % % % 
% Step IV - see figure 2 in the paper
% % % % % % % % % % % % % % % % % % % % % % % % % % 
% Keep only dierctionality vectors with reasonable score
find_possible_directionality_vectors
% Calculate confidence intervals for metabolite concentration rations
calc_confidence_intervals_for_metabolite_ratios
% Calculate confidence intervals for metabolite concentrations
calc_confidence_intervals_for_metabolite_concentrations
% Calculate confidence intervals for Gibbs free energies
calc_confidence_intervals_for_dG
% Calculate confidence intervals for fluxes
calc_confidence_intervals_for_fluxes


