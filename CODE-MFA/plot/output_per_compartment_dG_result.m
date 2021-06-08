% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% CODE-MFA derived Gibbs free energies for cytosolic and 
% mitochondrial reactions
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %

clear all;

load('../mat_files/sensitiviy_analysis_dG.mat', 'sensitiviy_analysis_dG');
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

addpath('../functions/emu') 
addpath('../functions/general') 

RT = 2.5;
% concentration in mM
CO2_con = 1.2;

[dTransporters,tTransporters] = xlsread('../xls_input_files/reactions_per_compartment.xlsx','transporters');
[dBothCompartments,tBothCompartments] = xlsread('../xls_input_files/reactions_per_compartment.xlsx','both compartments');
[dCytosolic,tCytosolic] = xlsread('../xls_input_files/reactions_per_compartment.xlsx','cytosolic');
[dMitochondrial,tMitochondrial] = xlsread('../xls_input_files/reactions_per_compartment.xlsx','mitochondrial');

all_met_names_with_WC_measured_con = model_thermodynamics.WC.met_name;
all_met_names_with_WC_measured_con{end+1} = 'CO2';
all_met_WC_measured_con = model_thermodynamics.WC.Concentrations;
all_met_WC_measured_con(end+1) = CO2_con;
% add some metabolite concentrations that I put in the compartmentalized
% section
all_met_names_with_WC_measured_con{end+1} = '6phosphogluconate';
index_compartmentalized_con = find(ismember(model_thermodynamics.mets,'6phosphogluconate_CY'))
all_met_WC_measured_con(end+1) = (exp(model_thermodynamics.mets_lb(index_compartmentalized_con))+exp(model_thermodynamics.mets_ub(index_compartmentalized_con)))/2
all_met_names_with_WC_measured_con{end+1} = 'Glc6P';
index_compartmentalized_con = find(ismember(model_thermodynamics.mets,'Glc6P_CY'))
all_met_WC_measured_con(end+1) = (exp(model_thermodynamics.mets_lb(index_compartmentalized_con))+exp(model_thermodynamics.mets_ub(index_compartmentalized_con)))/2
all_met_names_with_WC_measured_con{end+1} = 'Ribose5P';
index_compartmentalized_con = find(ismember(model_thermodynamics.mets,'Ribose5P_CY'))
all_met_WC_measured_con(end+1) = (exp(model_thermodynamics.mets_lb(index_compartmentalized_con))+exp(model_thermodynamics.mets_ub(index_compartmentalized_con)))/2
all_met_names_with_WC_measured_con{end+1} = 'G3P';
index_compartmentalized_con = find(ismember(model_thermodynamics.mets,'G3P'))
all_met_WC_measured_con(end+1) = (exp(model_thermodynamics.mets_lb(index_compartmentalized_con))+exp(model_thermodynamics.mets_ub(index_compartmentalized_con)))/2
all_met_names_with_WC_measured_con{end+1} = '13BPG';
index_compartmentalized_con = find(ismember(model_thermodynamics.mets,'13BPG'))
all_met_WC_measured_con(end+1) = (exp(model_thermodynamics.mets_lb(index_compartmentalized_con))+exp(model_thermodynamics.mets_ub(index_compartmentalized_con)))/2
all_met_names_with_WC_measured_con{end+1} = '3PG';
index_compartmentalized_con = find(ismember(model_thermodynamics.mets,'3PG'))
all_met_WC_measured_con(end+1) = (exp(model_thermodynamics.mets_lb(index_compartmentalized_con))+exp(model_thermodynamics.mets_ub(index_compartmentalized_con)))/2

 
% Take dG of code-MFA result for CY and for MT
tBothCompartments = tBothCompartments(2:end,:);
for(i=1:size(tBothCompartments,1))
    index_rxns_cy = find(contains(model_thermodynamics.full_rxns,tBothCompartments(i,2)));
    index_rxns_mt = find(contains(model_thermodynamics.full_rxns,tBothCompartments(i,3)));
     
    % check if name in figure is the same or opposite direction to the
    % model
    if(dBothCompartments(i,1))
        code_mfa_dG_cy(i,:) = length(index_rxns_cy)*[sensitiviy_analysis_dG{index_rxns_cy(1)}.low_dG sensitiviy_analysis_dG{index_rxns_cy(1)}.high_dG];
    else
        code_mfa_dG_cy(i,:) = length(index_rxns_cy)*[-sensitiviy_analysis_dG{index_rxns_cy(1)}.high_dG -sensitiviy_analysis_dG{index_rxns_cy(1)}.low_dG];
    end
    if(dBothCompartments(i,2))
        code_mfa_dG_mt(i,:) = length(index_rxns_mt)*[sensitiviy_analysis_dG{index_rxns_mt(1)}.low_dG sensitiviy_analysis_dG{index_rxns_mt(1)}.high_dG];
    else
        code_mfa_dG_mt(i,:) = length(index_rxns_mt)*[-sensitiviy_analysis_dG{index_rxns_mt(1)}.high_dG -sensitiviy_analysis_dG{index_rxns_mt(1)}.low_dG];
    end 
     
    % remove _CY _MT _source _sink from metabolite names
    str_of_substrates_and_producs_by_cy_reaction = split(model_thermodynamics.full_rxns{index_rxns_cy(1)},' => ');
    str_of_substrates_by_cy_reaction    = str_of_substrates_and_producs_by_cy_reaction(1);
    str_of_substrates_by_cy_reaction    = split(str_of_substrates_by_cy_reaction,' + ');
    str_of_substrates_by_cy_reaction    = strrep(str_of_substrates_by_cy_reaction,'_source','');
    str_of_substrates_by_cy_reaction    = strrep(str_of_substrates_by_cy_reaction,'_sink','');
    str_of_substrates_by_cy_reaction    = strrep(str_of_substrates_by_cy_reaction,'_CY','');
    str_of_products_by_cy_reaction      = str_of_substrates_and_producs_by_cy_reaction(2);
    str_of_products_by_cy_reaction      = split(str_of_products_by_cy_reaction,' + ');
    str_of_products_by_cy_reaction      = strrep(str_of_products_by_cy_reaction,'_source','');
    str_of_products_by_cy_reaction      = strrep(str_of_products_by_cy_reaction,'_sink','');
    str_of_products_by_cy_reaction      = strrep(str_of_products_by_cy_reaction,'_CY','');
    
    % calculate Gibbs free energy based on WC level metabolites
    % concentration measurements
    index_of_products_in_wc_measured_vector = find(ismember(all_met_names_with_WC_measured_con,str_of_products_by_cy_reaction));
    if(length(index_of_products_in_wc_measured_vector)==length(str_of_products_by_cy_reaction))
        products_mul = prod((all_met_WC_measured_con(index_of_products_in_wc_measured_vector)));
    else
        products_mul = nan;
    end
    index_of_substrates_in_wc_measured_vector = find(ismember(all_met_names_with_WC_measured_con,str_of_substrates_by_cy_reaction));
    if(length(index_of_substrates_in_wc_measured_vector)==length(str_of_substrates_by_cy_reaction))
        substrates_mul = prod((all_met_WC_measured_con(index_of_substrates_in_wc_measured_vector)));
    else
        substrates_mul = nan;
    end 
     
    code_mfa_dG_per_wc_measurements(i) = model_thermodynamics.delta_G0(index_rxns_cy(1))+RT*log(products_mul/substrates_mul);
    % switch direction if the string is in the opposite direction
    if(~dBothCompartments(i,1))
        code_mfa_dG_per_wc_measurements(i) = -code_mfa_dG_per_wc_measurements(i);
    end
end


[val ind]=sort(code_mfa_dG_cy(:,1));
code_mfa_dG_cy = code_mfa_dG_cy(ind,:);
code_mfa_dG_mt = code_mfa_dG_mt(ind,:);
rxns_text_for_both_compartments = tBothCompartments(ind,1)
code_mfa_dG_per_wc_measurements = code_mfa_dG_per_wc_measurements(ind);


figure;
CONST_VAL_FOR_PLOT_MT_ABOVE_CY=0.08;
subplot(2,2,1);
h1=plot(code_mfa_dG_cy',[1-CONST_VAL_FOR_PLOT_MT_ABOVE_CY:1:length(code_mfa_dG_cy);1-CONST_VAL_FOR_PLOT_MT_ABOVE_CY:1:length(code_mfa_dG_cy)],'Color',[0 0.44 0.74], 'LineWidth',4);
hold on;
h2=plot(code_mfa_dG_mt',[1+CONST_VAL_FOR_PLOT_MT_ABOVE_CY:1:length(code_mfa_dG_mt)+CONST_VAL_FOR_PLOT_MT_ABOVE_CY;1+CONST_VAL_FOR_PLOT_MT_ABOVE_CY:1:length(code_mfa_dG_mt)+CONST_VAL_FOR_PLOT_MT_ABOVE_CY],'Color',[0.4660 0.6740 0.1880], 'LineWidth',4);
h3=plot(code_mfa_dG_per_wc_measurements',[1:1:length(code_mfa_dG_per_wc_measurements);1:1:length(rxns_text_for_both_compartments(ind))],'k*','MarkerSize',15);
h4_dummy=plot(100,100,'Color',[0.8500 0.3250 0.0980],'LineWidth',4);
line([0 0], [0 length(rxns_text_for_both_compartments)+1], 'color','red','LineStyle',':');
ylim([0 size(rxns_text_for_both_compartments,1)+1]);
% xlim([-70 10]);
YLabel={''};
for(i=1:size(rxns_text_for_both_compartments,1))
    YLabel{end+1}=rxns_text_for_both_compartments{i};
    YLabel{end} = strrep(YLabel{end},'_',' ');
end
YLabel{end+1}=''; 
set(gca, 'YTickLabel', YLabel);  
set(gca, 'yTick', [0:length(YLabel)+1], 'FontSize',12);
xlabel('Gibbs free energy [kJ/mol]', 'FontSize',11);
title('Cytosolic & mitochondrial reactions', 'FontSize',15);
set(gcf,'color','w');
ax = gca;
ax.Box='on';
ax.BoxStyle = 'full';
grid;
% legend([h1(1),h2(1),h4_dummy(1),h3(1)],'Cytosolic reactions','Mitochondrial reactions','Transporrters','Measured (WC-level)','FontSize',18,'Orientation','horizontal');
legend([h1(1),h2(1),h4_dummy(1),h3(1)],'Cytosolic reactions','Mitochondrial reactions','Transporters','Measured (WC-level)','FontSize',18);
hold off;

 
% plot dG of transporters
code_mfa_dG = [];
for(i=1:size(tTransporters,1))
    index_rxns = find(contains(model_thermodynamics.full_rxns,tTransporters(i)));
    code_mfa_dG(i,:) = length(index_rxns)*[sensitiviy_analysis_dG{index_rxns(1)}.low_dG sensitiviy_analysis_dG{index_rxns(1)}.high_dG];
      
    % remove _CY _MT _source _sink from metabolite names
    str_of_substrates_and_producs_by_reaction = split(model_thermodynamics.full_rxns{index_rxns(1)},' => ');
    str_of_substrates_by_reaction    = str_of_substrates_and_producs_by_reaction(1);
    str_of_substrates_by_reaction    = split(str_of_substrates_by_reaction,' + ');
    str_of_substrates_by_reaction    = strrep(str_of_substrates_by_reaction,'_source','');
    str_of_substrates_by_reaction    = strrep(str_of_substrates_by_reaction,'_sink','');
    str_of_substrates_by_reaction    = strrep(str_of_substrates_by_reaction,'_CY','');
    str_of_substrates_by_reaction    = strrep(str_of_substrates_by_reaction,'_MT','');
    str_of_products_by_reaction      = str_of_substrates_and_producs_by_reaction(2);
    str_of_products_by_reaction      = split(str_of_products_by_reaction,' + ');
    str_of_products_by_reaction      = strrep(str_of_products_by_reaction,'_source','');
    str_of_products_by_reaction      = strrep(str_of_products_by_reaction,'_sink','');
    str_of_products_by_reaction      = strrep(str_of_products_by_reaction,'_CY','');
    str_of_products_by_reaction      = strrep(str_of_products_by_reaction,'_MT','');
    
    % calculate Gibbs free energy based on WC level metabolites
    % concentration measurements
    index_of_products_in_wc_measured_vector = find(ismember(all_met_names_with_WC_measured_con,str_of_products_by_reaction));
    if(length(index_of_products_in_wc_measured_vector)==length(str_of_products_by_reaction))
        products_mul = prod((all_met_WC_measured_con(index_of_products_in_wc_measured_vector)));
    else
        products_mul = nan;
    end
    index_of_substrates_in_wc_measured_vector = find(ismember(all_met_names_with_WC_measured_con,str_of_substrates_by_reaction));
    if(length(index_of_substrates_in_wc_measured_vector)==length(str_of_substrates_by_reaction))
        substrates_mul = prod((all_met_WC_measured_con(index_of_substrates_in_wc_measured_vector)));
    else
        substrates_mul = nan;
    end 
     
    code_mfa_dG_per_wc_measurements(i) = model_thermodynamics.delta_G0(index_rxns(1))+RT*log(products_mul/substrates_mul);
end

[val ind]=sort(code_mfa_dG(:,1));
code_mfa_dG = code_mfa_dG(ind,:);
rxns_text = tTransporters(ind,1)
code_mfa_dG_per_wc_measurements = code_mfa_dG_per_wc_measurements(ind);


subplot(2,2,2);
plot(code_mfa_dG',[1:length(code_mfa_dG);1:length(code_mfa_dG)],'Color',[0.8500 0.3250 0.0980], 'LineWidth',4);
hold on;
plot(code_mfa_dG_per_wc_measurements',[1:1:length(code_mfa_dG_per_wc_measurements);1:1:length(code_mfa_dG_per_wc_measurements)],'k*','MarkerSize',15);
line([0 0], [0 length(code_mfa_dG_per_wc_measurements)+1], 'color','red','LineStyle',':');
ylim([0 length(code_mfa_dG_per_wc_measurements)+1]);
% xlim([-70 10]);
YLabel={''};
for(i=1:size(rxns_text,1))
    YLabel{end+1}=rxns_text{i};
    YLabel{end} = strrep(YLabel{end},'_',' ');
end
YLabel{end+1}='';
set(gca, 'YTickLabel', YLabel);  
set(gca, 'yTick', [0:length(YLabel)+1], 'FontSize',12);
xlabel('Gibbs free energy [kJ/mol]', 'FontSize',11);
title('Transporters', 'FontSize',15);
set(gcf,'color','w');
ax = gca;
ax.Box='on';
ax.BoxStyle = 'full';
grid;
hold off;



% plot dG of cytosolic reactions
tCytosolic = tCytosolic(2:end,:);
code_mfa_dG = [];
for(i=1:size(tCytosolic,1))
    index_rxns = find(contains(model_thermodynamics.full_rxns,tCytosolic(i,2)));
    code_mfa_dG(i,:) = length(index_rxns)*[sensitiviy_analysis_dG{index_rxns(1)}.low_dG sensitiviy_analysis_dG{index_rxns(1)}.high_dG];
      
    % remove _CY _MT _source _sink from metabolite names
    str_of_substrates_and_producs_by_reaction = split(model_thermodynamics.full_rxns{index_rxns(1)},' => ');
    str_of_substrates_by_reaction    = str_of_substrates_and_producs_by_reaction(1);
    str_of_substrates_by_reaction    = split(str_of_substrates_by_reaction,' + ');
    str_of_substrates_by_reaction    = strrep(str_of_substrates_by_reaction,'_source','');
    str_of_substrates_by_reaction    = strrep(str_of_substrates_by_reaction,'_sink','');
    str_of_substrates_by_reaction    = strrep(str_of_substrates_by_reaction,'_CY','');
    str_of_substrates_by_reaction    = strrep(str_of_substrates_by_reaction,'_MT','');
    str_of_products_by_reaction      = str_of_substrates_and_producs_by_reaction(2);
    str_of_products_by_reaction      = split(str_of_products_by_reaction,' + ');
    str_of_products_by_reaction      = strrep(str_of_products_by_reaction,'_source','');
    str_of_products_by_reaction      = strrep(str_of_products_by_reaction,'_sink','');
    str_of_products_by_reaction      = strrep(str_of_products_by_reaction,'_CY','');
    str_of_products_by_reaction      = strrep(str_of_products_by_reaction,'_MT','');
    
    % calculate Gibbs free energy based on WC level metabolites
    % concentration measurements
    index_of_products_in_wc_measured_vector = find(ismember(all_met_names_with_WC_measured_con,str_of_products_by_reaction));
    if(length(index_of_products_in_wc_measured_vector)==length(str_of_products_by_reaction))
        products_mul = prod((all_met_WC_measured_con(index_of_products_in_wc_measured_vector)));
    else
        products_mul = nan;
    end
    index_of_substrates_in_wc_measured_vector = find(ismember(all_met_names_with_WC_measured_con,str_of_substrates_by_reaction));
    if(length(index_of_substrates_in_wc_measured_vector)==length(str_of_substrates_by_reaction))
        substrates_mul = prod((all_met_WC_measured_con(index_of_substrates_in_wc_measured_vector)));
    else
        substrates_mul = nan;
    end  
     
    code_mfa_dG_per_wc_measurements(i) = model_thermodynamics.delta_G0(index_rxns(1))+RT*log(products_mul/substrates_mul);
end

[val ind]=sort(code_mfa_dG(:,1));
code_mfa_dG = code_mfa_dG(ind,:);
rxns_text = tCytosolic(ind,1)
code_mfa_dG_per_wc_measurements = code_mfa_dG_per_wc_measurements(ind);


subplot(2,2,3);
plot(code_mfa_dG',[1:length(code_mfa_dG);1:length(code_mfa_dG)],'Color',[0 0.44 0.74], 'LineWidth',4);
hold on;
plot(code_mfa_dG_per_wc_measurements',[1:1:length(code_mfa_dG_per_wc_measurements);1:1:length(code_mfa_dG_per_wc_measurements)],'k*','MarkerSize',15);
line([0 0], [0 length(code_mfa_dG_per_wc_measurements)+1], 'color','red','LineStyle',':');
ylim([0 length(code_mfa_dG_per_wc_measurements)+1]);
% xlim([-70 10]);
YLabel={''};
for(i=1:size(rxns_text,1))
    YLabel{end+1}=rxns_text{i};
    YLabel{end} = strrep(YLabel{end},'_',' ');
end
YLabel{end+1}='';
set(gca, 'YTickLabel', YLabel);  
set(gca, 'yTick', [0:length(YLabel)+1], 'FontSize',12);
xlabel('Gibbs free energy [kJ/mol]', 'FontSize',11);
title('Cytosolic reactions', 'FontSize',15);
set(gcf,'color','w');
ax = gca;
ax.Box='on';
ax.BoxStyle = 'full';
grid;
hold off;



% plot dG of mitochondrial reactions
tMitochondrial = tMitochondrial(2:end,:);
code_mfa_dG = [];
for(i=1:size(tMitochondrial,1))
    index_rxns = find(contains(model_thermodynamics.full_rxns,tMitochondrial(i,2)));
    code_mfa_dG(i,:) = length(index_rxns)*[sensitiviy_analysis_dG{index_rxns(1)}.low_dG sensitiviy_analysis_dG{index_rxns(1)}.high_dG];
      
    % remove _CY _MT _source _sink from metabolite names
    str_of_substrates_and_producs_by_reaction = split(model_thermodynamics.full_rxns{index_rxns(1)},' => ');
    str_of_substrates_by_reaction    = str_of_substrates_and_producs_by_reaction(1);
    str_of_substrates_by_reaction    = split(str_of_substrates_by_reaction,' + ');
    str_of_substrates_by_reaction    = strrep(str_of_substrates_by_reaction,'_source','');
    str_of_substrates_by_reaction    = strrep(str_of_substrates_by_reaction,'_sink','');
    str_of_substrates_by_reaction    = strrep(str_of_substrates_by_reaction,'_CY','');
    str_of_substrates_by_reaction    = strrep(str_of_substrates_by_reaction,'_MT','');
    str_of_products_by_reaction      = str_of_substrates_and_producs_by_reaction(2);
    str_of_products_by_reaction      = split(str_of_products_by_reaction,' + ');
    str_of_products_by_reaction      = strrep(str_of_products_by_reaction,'_source','');
    str_of_products_by_reaction      = strrep(str_of_products_by_reaction,'_sink','');
    str_of_products_by_reaction      = strrep(str_of_products_by_reaction,'_CY','');
    str_of_products_by_reaction      = strrep(str_of_products_by_reaction,'_MT','');
    
    % calculate Gibbs free energy based on WC level metabolites
    % concentration measurements
    index_of_products_in_wc_measured_vector = find(ismember(all_met_names_with_WC_measured_con,str_of_products_by_reaction));
    if(length(index_of_products_in_wc_measured_vector)==length(str_of_products_by_reaction))
        products_mul = prod((all_met_WC_measured_con(index_of_products_in_wc_measured_vector)));
    else
        products_mul = nan;
    end
    index_of_substrates_in_wc_measured_vector = find(ismember(all_met_names_with_WC_measured_con,str_of_substrates_by_reaction));
    if(length(index_of_substrates_in_wc_measured_vector)==length(str_of_substrates_by_reaction))
        substrates_mul = prod((all_met_WC_measured_con(index_of_substrates_in_wc_measured_vector)));
    else
        substrates_mul = nan;
    end  
     
    code_mfa_dG_per_wc_measurements(i) = model_thermodynamics.delta_G0(index_rxns(1))+RT*log(products_mul/substrates_mul);
end

[val ind]=sort(code_mfa_dG(:,1));
code_mfa_dG = code_mfa_dG(ind,:);
rxns_text = tMitochondrial(ind,1)
code_mfa_dG_per_wc_measurements = code_mfa_dG_per_wc_measurements(ind);


subplot(2,2,4);
plot(code_mfa_dG',[1:length(code_mfa_dG);1:length(code_mfa_dG)],'Color',[0.4660 0.6740 0.1880], 'LineWidth',4);
hold on;
plot(code_mfa_dG_per_wc_measurements',[1:1:length(code_mfa_dG_per_wc_measurements);1:1:length(code_mfa_dG_per_wc_measurements)],'k*','MarkerSize',15);
line([0 0], [0 length(code_mfa_dG_per_wc_measurements)+1], 'color','red','LineStyle',':');
ylim([0 length(code_mfa_dG_per_wc_measurements)+1]);
% xlim([-70 10]);
YLabel={''};
for(i=1:size(rxns_text,1))
    YLabel{end+1}=rxns_text{i};
    YLabel{end} = strrep(YLabel{end},'_',' ');
end
YLabel{end+1}='';
set(gca, 'YTickLabel', YLabel);  
set(gca, 'yTick', [0:length(YLabel)+1], 'FontSize',12);
xlabel('Gibbs free energy [kJ/mol]', 'FontSize',11);
title('Mitochondial reactions', 'FontSize',15);
set(gcf,'color','w');
ax = gca;
ax.Box='on';
ax.BoxStyle = 'full';
grid;
hold off;






