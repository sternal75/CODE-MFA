% % % % % % % % % % % % % % % % % % % % % 
% load model including carbon mapping
% % % % % % % % % % % % % % % % % % % % % 
close all;

% readCbModel - reads input model file in SBML or excell format 
% cobratoolbox
cobraToolBoxModel=readCbModel('xls_input_files/input.xlsx');
[dR,tR] = xlsread('xls_input_files/input.xlsx');
[dM,tM] = xlsread('xls_input_files/input.xlsx','Metabolite List')


% translate to the EMU code model
model.mets      = strrep(cobraToolBoxModel.mets,'[c]','')';
model.metNames  = model.mets;
model.rxns      = cobraToolBoxModel.rxns';
model.rxn_formula = tR(2:end,3);
% upload lower and upper bound values
model.positive_direction_lb=cobraToolBoxModel.lb';
model.positive_direction_lb(model.positive_direction_lb<0)=0.001;
model.positive_direction_ub=cobraToolBoxModel.ub';
if(model.positive_direction_ub<0)
    exit();
end
model.negative_direction_lb=ones(length(cobraToolBoxModel.lb),1)*0.001;
model.negative_direction_ub=abs(cobraToolBoxModel.lb');
model.negative_direction_lb(cobraToolBoxModel.lb>0)=nan;
model.negative_direction_ub(cobraToolBoxModel.lb>0)=nan;


% stoichiometric matrix S
model.S         = full(cobraToolBoxModel.S);
[model.met_num model.rxn_num] = size(model.S);
model.used_reactions_status = [1:length(model.rxns)]';  %?? NOT USED??
model.exchange = zeros(model.rxn_num,1);%[0;0;0;0;0;0;0;0];     %?? NOT USED??
model.equality_constraints = [];
model.non_equality_constraint_values = [];
model.non_equality_constraint_fluxes = [];


model.atom_C_num=[];
model.met_extra=[];
model.met_extra_labeling=cell(0);
for(i=1:length(model.mets))
    % fill number of carbons for reaction
    idx=strfind(tR(:,3),model.mets{i});
    idx = ~cellfun('isempty',idx);
    met_line_number = find(idx);
    metabolites_in_reaction = strtrim(strsplit(tR{met_line_number(1),3},{' => ',' + '}));
    met_index_in_reaction = strmatch(model.mets(i),metabolites_in_reaction);
    carbon_mapping_in_reaction = strtrim(strsplit(tR{met_line_number(1),14},{' => ',' + '}));
    number_of_carbons = length(carbon_mapping_in_reaction{met_index_in_reaction(1)});
    model.atom_C_num=[model.atom_C_num;number_of_carbons];
    % fill external metabolite = 1
    idx=strfind(tM(:,1),model.mets{i});
    idx = ~cellfun('isempty',idx);
    met_line_number = find(idx);
    if(strcmp(tM(met_line_number,5),'Extra-organism'))
        model.met_extra=[model.met_extra;1];
        % glutamine extra-organism labeling
        external_labeling = strsplit(tM{met_line_number,6},'/');
        model.met_extra_labeling_glutamine{i}.atom_idv=str2num(external_labeling{1});
        model.met_extra_labeling_glutamine{i}.enrichment=str2num(external_labeling{2});
        % glucose extra-organism labeling
        external_labeling = strsplit(tM{met_line_number,7},'/');
        model.met_extra_labeling_glucose{i}.atom_idv=str2num(external_labeling{1});
        model.met_extra_labeling_glucose{i}.enrichment=str2num(external_labeling{2});        
    else
        model.met_extra=[model.met_extra;0];
        model.met_extra_labeling_glutamine{i}   = [];
        model.met_extra_labeling_glucose{i}     = [];
    end   
end

model.mappings_carbon = cell(0);
for(i=1:length(model.rxns))
    % init all matrices for carbon mapping
    model.mappings_carbon{end+1}.graph_r.node_info=[];
    model.mappings_carbon{end}.graph_p.node_info=[];
    model.mappings_carbon{end}.graph_r.mets=[];
    model.mappings_carbon{end}.graph_p.mets=[];
    model.mappings_carbon{end}.mapping_p=[];
    temp_reactants =[];
    temp_products  =[];
    % fill carbon mapping
    reactants_and_products = strsplit(tR{i+1,3},{' => '});
    reactants = strtrim(strsplit(reactants_and_products{1},{' + '}));
    products  = strtrim(strsplit(reactants_and_products{2},{' + '}));
    reactants_and_products_carbon_map = strsplit(tR{i+1,14},{' => '});
    reactants_carbon_map = strtrim(strsplit(reactants_and_products_carbon_map{1},{' + '}));
    products_carbon_map  = strtrim(strsplit(reactants_and_products_carbon_map{2},{' + '}));    
         
   
    % go over all reactants
    for(j=1:length(reactants))
        idx = strfind(model.mets,reactants{j});
        idx = ~cellfun('isempty',idx);
        met_index = find(idx);
        model.mappings_carbon{end}.graph_r.mets = [model.mappings_carbon{end}.graph_r.mets;met_index];
        model.mappings_carbon{end}.graph_r.node_info(end+1:end+model.atom_C_num(idx),1)=j;
        
        if(j==1)
            reactant_mapping=hex2dec(regexp(num2str(reactants_carbon_map{j}),'\w','match'))';
        else
            reactant_mapping=[reactant_mapping hex2dec(regexp(num2str(reactants_carbon_map{j}),'\w','match'))'+reactant_mapping(end)];
        end        
    end
    % go over all products
    for(j=1:length(products))
        idx = strfind(model.mets,products{j});
        idx = ~cellfun('isempty',idx);
        met_index = find(idx);       
        model.mappings_carbon{end}.graph_p.mets = [model.mappings_carbon{end}.graph_p.mets;met_index];
        model.mappings_carbon{end}.graph_p.node_info(end+1:end+model.atom_C_num(idx),1)=j;
        
        if(j==1)
            product_mapping=hex2dec(regexp(num2str(products_carbon_map{j}),'\w','match'))';
        else
            product_mapping=[product_mapping hex2dec(regexp(num2str(products_carbon_map{j}),'\w','match'))'+product_mapping(end)];
        end        
    end    
    reactant_mapping
    product_mapping
    model.mappings_carbon{end}.mapping_p(product_mapping,1)=reactant_mapping';
end

% handle equality constraints for reactions with equal rates
equal_reactions = unique(dR(:,10));
equal_reactions = equal_reactions(~isnan(equal_reactions))';
% equal_reactions = rmmissing(equal_reactions);
for(i=1:length(equal_reactions))
    
    reaction_indices = find(dR(:,10)==i);
    for(j=2:length(reaction_indices))
        addDuplicateReactionsToSMatrix=zeros(1,length(model.rxns));
        addDuplicateReactionsToSMatrix(reaction_indices(1)) = 1;
        addDuplicateReactionsToSMatrix(reaction_indices(j)) = -1;
        model.equality_constraints = [model.equality_constraints;addDuplicateReactionsToSMatrix];
    end
end

%force zero flux for glutamine
model.force_zero_flux{1} = dR(:,12);
%force zero flux for glucose
model.force_zero_flux{2} = dR(:,13);


model.skip_sensitivity_analysis_reaction_indices = dR(:,11);
