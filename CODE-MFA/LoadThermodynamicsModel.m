% % % % % % % % % % % % % % % % % % % % % % 
% load model including themodynamic input
% % % % % % % % % % % % % % % % % % % % % % 
load_constants

[dR,tR] = xlsread('xls_input_files/input_for_thermodynamics.xlsx');
[dM,tM] = xlsread('xls_input_files/input_for_thermodynamics.xlsx','Metabolite List');
[dMWCCon,tMWCCon] = xlsread('xls_input_files/input_for_thermodynamics.xlsx','WC Metabolite Concentration');
[dMCoFaRatios,tMCoFaRatios] = xlsread('xls_input_files/input_for_thermodynamics.xlsx','Metabolite co-factor ratios');
[dR1,tR1] = xlsread('xls_input_files/input.xlsx');
[dM1,tM1] = xlsread('xls_input_files/input.xlsx','Metabolite List')



model_thermodynamics.rxns = tR(2:end,1)';
model_thermodynamics.rxns = tR(2:end,2)';
model_thermodynamics.full_rxns = tR(2:end,3)';
model_thermodynamics.delta_G0 = dR(1:end,1);
model_thermodynamics.delta_G_low = dR(1:end,2);
model_thermodynamics.delta_G_high = dR(1:end,3);
model_thermodynamics.thermodynamics_of_reaction_defined = dR(1:end,4);
model_thermodynamics.thermodynamics_of_reaction_defined(isnan(model_thermodynamics.thermodynamics_of_reaction_defined)) = 0;
model_thermodynamics.mets = tM(2:end,1);
model_thermodynamics.mets_lb = dM(1:end,1);
model_thermodynamics.mets_ub = dM(1:end,2);
model_thermodynamics.WC.Concentrations = dMWCCon(1:end,1);
model_thermodynamics.WC.Concentrations_STD = dMWCCon(1:end,2);
model_thermodynamics.WC.met_name = tMWCCon(2:end,1);
model_thermodynamics.reactant_indexes=cell(length(model_thermodynamics.rxns),1);
model_thermodynamics.product_indexes=cell(length(model_thermodynamics.rxns),1);

for(i=1:length(model_thermodynamics.rxns))
	% two way reaction
        % full reaction exists - thermodynamic 
        if((~isempty(model_thermodynamics.full_rxns{i})))
            reactants_and_products=strsplit(model_thermodynamics.full_rxns{i},' => ');
            reactants=strsplit(reactants_and_products{1},' + ');
            products=strsplit(reactants_and_products{2},' + ');
            for(j=1:length(reactants))
                indexC  = strcmp(model_thermodynamics.mets,reactants{j});
                index = find(indexC==1);
                if(length(index)~=1)
                    disp('An error occurred - metabolite was not declared');
                    exit();
                end            
                model_thermodynamics.reactant_indexes{i}=[model_thermodynamics.reactant_indexes{i};index];
            end
            for(j=1:length(products))
                indexC  = strcmp(model_thermodynamics.mets,products{j});
                index = find(indexC==1);
                if(length(index)~=1)
                    disp('An error occurred - metabolite was not declared');
                    exit();
                end            
                model_thermodynamics.product_indexes{i}=[model_thermodynamics.product_indexes{i};index];
            end
        end
end
compartmenalized_met_indices_CY=zeros(length(model_thermodynamics.WC.Concentrations_STD), length(model_thermodynamics.mets));
compartmenalized_met_indices_MT=zeros(length(model_thermodynamics.WC.Concentrations_STD), length(model_thermodynamics.mets));
for(i=1:length(model_thermodynamics.WC.Concentrations_STD))
    met_name = model_thermodynamics.WC.met_name{i};
    met_name_CY = strcat(met_name,'_CY');
    met_name_MT = strcat(met_name,'_MT');    
    index_MT=strcmp(model_thermodynamics.mets,met_name_MT)
    index_MT=find(index_MT==1);
    index_CY=strcmp(model_thermodynamics.mets,met_name_CY)
    index_CY=find(index_CY==1);
    compartmenalized_met_indices_CY(i,index_CY)=1;
    compartmenalized_met_indices_MT(i,index_MT)=1;
    % replace the MT/CY upper bound concentration of a metabolite if the
    % measured upper bound is better(smaller) than given input upper bound
    measured_MT_ub = log((model_thermodynamics.WC.Concentrations(i)+2*model_thermodynamics.WC.Concentrations_STD(i))/MT_WC_VOLUME);
    measured_CY_ub = log((model_thermodynamics.WC.Concentrations(i)+2*model_thermodynamics.WC.Concentrations_STD(i))/CY_WC_VOLUME);   
    if(model_thermodynamics.mets_ub(index_MT) > measured_MT_ub)
        model_thermodynamics.mets_ub(index_MT) = measured_MT_ub;
    end
    if(model_thermodynamics.mets_ub(index_CY) > measured_CY_ub)
        model_thermodynamics.mets_ub(index_CY) = measured_CY_ub;
    end    
end
model_thermodynamics.WC.compartmenalized_met_indices_CY = compartmenalized_met_indices_CY;
model_thermodynamics.WC.compartmenalized_met_indices_MT = compartmenalized_met_indices_MT;

model_thermodynamics.co_factors=cell(0);
co_factors = tMCoFaRatios(2:end,1);
for(i=1:length(co_factors))
    metabolites = strsplit(co_factors{i},'/');
    index_metabolite_1=strcmp(model_thermodynamics.mets,metabolites{1});
    index_metabolite_1=find(index_metabolite_1==1);
    index_metabolite_2=strcmp(model_thermodynamics.mets,metabolites{2});
    index_metabolite_2=find(index_metabolite_2==1);
    model_thermodynamics.co_factors{end+1}.indices=[index_metabolite_1 index_metabolite_2];
    model_thermodynamics.co_factors{end}.ratio_lb=dMCoFaRatios(i,1);
    model_thermodynamics.co_factors{end}.ratio_ub=dMCoFaRatios(i,2);
end
model_thermodynamics.skip_sensitivity_analysis_metabolite_indices = dM(:,5);