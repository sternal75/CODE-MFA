% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% CODE-MFA derived Gibbs free energies and net fluxes for cytosolic and 
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
load('../mat_files/net_fluxes.mat','net_fluxes');


addpath('../functions/emu') 
addpath('../functions/general') 

run ../load_constants;

RT = 2.5;
% concentration in mM
CO2_con = 1.2;

[dTransporters,tTransporters] = xlsread('../xls_input_files/all_reactions_per_compartment.xlsx','transporters');
[dBothCompartments,tBothCompartments] = xlsread('../xls_input_files/all_reactions_per_compartment.xlsx','both compartments');
[dCytosolic,tCytosolic] = xlsread('../xls_input_files/all_reactions_per_compartment.xlsx','cytosolic');
[dMitochondrial,tMitochondrial] = xlsread('../xls_input_files/all_reactions_per_compartment.xlsx','mitochondrial');


% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% Code-MFA dG
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
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
    index_rxns_cy = find(contains(model_thermodynamics.full_rxns,tBothCompartments(i,3)));
    index_rxns_mt = find(contains(model_thermodynamics.full_rxns,tBothCompartments(i,4)));
     
    % check if name in figure is the same or opposite direction to the
    % model
    if(model_thermodynamics.thermodynamics_of_reaction_defined(index_rxns_cy))
        if(dBothCompartments(i,1))
            code_mfa_dG_cy(i,:) = length(index_rxns_cy)*[sensitiviy_analysis_dG{index_rxns_cy(1)}.low_dG sensitiviy_analysis_dG{index_rxns_cy(1)}.high_dG];
        else
            code_mfa_dG_cy(i,:) = length(index_rxns_cy)*[-sensitiviy_analysis_dG{index_rxns_cy(1)}.high_dG -sensitiviy_analysis_dG{index_rxns_cy(1)}.low_dG];
        end
    else
        code_mfa_dG_cy(i,:) = [nan nan];
    end
    
    if(model_thermodynamics.thermodynamics_of_reaction_defined(index_rxns_mt))
        if(dBothCompartments(i,2))
            code_mfa_dG_mt(i,:) = length(index_rxns_mt)*[sensitiviy_analysis_dG{index_rxns_mt(1)}.low_dG sensitiviy_analysis_dG{index_rxns_mt(1)}.high_dG];
        else
            code_mfa_dG_mt(i,:) = length(index_rxns_mt)*[-sensitiviy_analysis_dG{index_rxns_mt(1)}.high_dG -sensitiviy_analysis_dG{index_rxns_mt(1)}.low_dG];
        end 
    else
        code_mfa_dG_mt(i,:) = [nan nan];
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


[val ind_sort_cy_mt_reactions]=sort(code_mfa_dG_cy(:,1));
code_mfa_dG_cy = code_mfa_dG_cy(ind_sort_cy_mt_reactions,:);
code_mfa_dG_mt = code_mfa_dG_mt(ind_sort_cy_mt_reactions,:);
rxns_text_for_both_compartments = tBothCompartments(ind_sort_cy_mt_reactions,1)
enzyme_text = tBothCompartments(ind_sort_cy_mt_reactions,2);
code_mfa_dG_per_wc_measurements = code_mfa_dG_per_wc_measurements(ind_sort_cy_mt_reactions);


figure('Position',[1 1 800 1000]);
CONST_VAL_FOR_PLOT_MT_ABOVE_CY=0.08;
subplot(4,1,1);
h1=plot(code_mfa_dG_cy',[1-CONST_VAL_FOR_PLOT_MT_ABOVE_CY:1:length(code_mfa_dG_cy);1-CONST_VAL_FOR_PLOT_MT_ABOVE_CY:1:length(code_mfa_dG_cy)],'Color',BLUE_COLOR, 'LineWidth',4);
hold on;
h2=plot(code_mfa_dG_mt',[1+CONST_VAL_FOR_PLOT_MT_ABOVE_CY:1:length(code_mfa_dG_mt)+CONST_VAL_FOR_PLOT_MT_ABOVE_CY;1+CONST_VAL_FOR_PLOT_MT_ABOVE_CY:1:length(code_mfa_dG_mt)+CONST_VAL_FOR_PLOT_MT_ABOVE_CY],'Color',GREEN_COLOR, 'LineWidth',4);
h3=plot(code_mfa_dG_per_wc_measurements',[1:1:length(code_mfa_dG_per_wc_measurements);1:1:length(rxns_text_for_both_compartments(ind_sort_cy_mt_reactions))],'k*','MarkerSize',12);
h4_dummy=plot(100,100,'Color',RED_COLOR,'LineWidth',4);
line([0 0], [0 length(rxns_text_for_both_compartments)+1], 'color','red','LineStyle',':');
ylim([0 size(rxns_text_for_both_compartments,1)+1]);
% xlim([-70 10]);
YLabel={''};
YLabel_right={''};
for(i=1:size(rxns_text_for_both_compartments,1))
    YLabel{end+1}=rxns_text_for_both_compartments{i};
    YLabel{end} = strrep(YLabel{end},'_',' ');
    YLabel_right{end+1} = enzyme_text{i};
end
YLabel{end+1}=''; 
set(gca, 'YTickLabel', YLabel);  
set(gca, 'yTick', [0:length(YLabel)+1]);
xlabel('Gibbs free energy [kJ/mol]');
title('Cytosolic & mitochondrial reactions');
set(gcf,'color','w');
ax = gca;
ax.Box='on';
ax.BoxStyle = 'full';
grid;
% legend([h1(1),h2(1),h4_dummy(1),h3(1)],'Cytosolic reactions','Mitochondrial reactions','Transporrters','Measured (WC-level)','FontSize',18,'Orientation','horizontal');
legend([h1(1),h2(1),h4_dummy(1),h3(1)],'Cytosolic','Mitochondrial','Transporters','Measured (WC)','FontSize',12, 'Orientation','horizontal', 'Location', 'north');
xlim([-40 40]);
xticks([-40:10:40]);
set(gca, 'FontSize',10);
yyaxis right;
ylim([0 length(code_mfa_dG_per_wc_measurements)+1]);
set(gca, 'YTickLabel', YLabel_right,'YColor','black');
set(gca, 'yTick', [0:length(YLabel)+1]);
hold off;

 


% plot dG of cytosolic reactions
tCytosolic = tCytosolic(2:end,:);
code_mfa_dG = [];
for(i=1:size(tCytosolic,1))
    index_rxns = find(contains(model_thermodynamics.full_rxns,tCytosolic(i,3)));
    
    % change the net flux of the reaction to be positive
    direction_of_reaction = 1;
    if(net_fluxes{index_rxns(1)}.high_flux < 0)
        direction_of_reaction = 0;
    end
    
    if(model_thermodynamics.thermodynamics_of_reaction_defined(index_rxns))
        code_mfa_dG(i,:) = length(index_rxns)*[sensitiviy_analysis_dG{index_rxns(1)}.low_dG sensitiviy_analysis_dG{index_rxns(1)}.high_dG];
    else
        code_mfa_dG(i,:) = [nan nan];
    end
      
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
    
    if(direction_of_reaction == 0)
        code_mfa_dG(i,[1 2]) = -code_mfa_dG(i,[2 1]);
        code_mfa_dG_per_wc_measurements(i) = -code_mfa_dG_per_wc_measurements(i);
        str_of_substrates_and_producs_by_reaction = split(tCytosolic(i,1),' => ');
        tCytosolic(i,1)={[str_of_substrates_and_producs_by_reaction{2} ' => ' str_of_substrates_and_producs_by_reaction{1}]};
    end
end

[val ind_sort_cy_reactions]=sort(code_mfa_dG(:,1));
code_mfa_dG = code_mfa_dG(ind_sort_cy_reactions,:);
rxns_text = tCytosolic(ind_sort_cy_reactions,1);
enzyme_text = tCytosolic(ind_sort_cy_reactions,2);
code_mfa_dG_per_wc_measurements = code_mfa_dG_per_wc_measurements(ind_sort_cy_reactions);


subplot(4,1,2);
plot(code_mfa_dG',[1:length(code_mfa_dG);1:length(code_mfa_dG)],'Color',BLUE_COLOR, 'LineWidth',4);
hold on;
plot(code_mfa_dG_per_wc_measurements',[1:1:length(code_mfa_dG_per_wc_measurements);1:1:length(code_mfa_dG_per_wc_measurements)],'k*','MarkerSize',12);
line([0 0], [0 length(code_mfa_dG_per_wc_measurements)+1], 'color','red','LineStyle',':');
ylim([0 length(code_mfa_dG_per_wc_measurements)+1]);
% xlim([-70 10]);
YLabel={''};
YLabel_right={''};
for(i=1:size(rxns_text,1))
    YLabel{end+1}=rxns_text{i};
    YLabel{end} = strrep(YLabel{end},'_',' ');
    YLabel_right{end+1} = enzyme_text{i};
end
YLabel{end+1}='';
YLabel_right{end+1}='';
set(gca, 'YTickLabel', YLabel);  
set(gca, 'yTick', [0:length(YLabel)+1]);
xlabel('Gibbs free energy [kJ/mol]');
title('Cytosolic reactions');
set(gcf,'color','w');
ax = gca;
ax.Box='on';
ax.BoxStyle = 'full';
grid;
xlim([-40 40]);
xticks([-40:10:40]);
set(gca, 'FontSize',10);
yyaxis right;
ylim([0 length(code_mfa_dG_per_wc_measurements)+1]);
set(gca, 'YTickLabel', YLabel_right,'YColor','black');
set(gca, 'yTick', [0:length(YLabel)+1]);
hold off;



% plot dG of mitochondrial reactions
tMitochondrial = tMitochondrial(2:end,:);
code_mfa_dG = [];
for(i=1:size(tMitochondrial,1))
    index_rxns = find(contains(model_thermodynamics.full_rxns,tMitochondrial(i,3)));
    
    % change the net flux of the reaction to be positive
    direction_of_reaction = 1;
    if(net_fluxes{index_rxns(1)}.high_flux < 0)
        direction_of_reaction = 0;
    end
    
    
    if(model_thermodynamics.thermodynamics_of_reaction_defined(index_rxns))
        code_mfa_dG(i,:) = length(index_rxns)*[sensitiviy_analysis_dG{index_rxns(1)}.low_dG sensitiviy_analysis_dG{index_rxns(1)}.high_dG];
    else
        code_mfa_dG(i,:) = [nan nan];
    end
    
    
      
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
    
    if(direction_of_reaction == 0)
        code_mfa_dG(i,[1 2]) = -code_mfa_dG(i,[2 1]);
        code_mfa_dG_per_wc_measurements(i) = -code_mfa_dG_per_wc_measurements(i);
        str_of_substrates_and_producs_by_reaction = split(tMitochondrial(i,1),' => ');
        tMitochondrial(i,1)={[str_of_substrates_and_producs_by_reaction{2} ' => ' str_of_substrates_and_producs_by_reaction{1}]};      
    end    
end

[val ind_sort_mt_reactions]=sort(code_mfa_dG(:,1));
code_mfa_dG = code_mfa_dG(ind_sort_mt_reactions,:);
rxns_text = tMitochondrial(ind_sort_mt_reactions,1)
enzyme_text = tMitochondrial(ind_sort_mt_reactions,2);
code_mfa_dG_per_wc_measurements = code_mfa_dG_per_wc_measurements(ind_sort_mt_reactions);


subplot(4,1,3);
plot(code_mfa_dG',[1:length(code_mfa_dG);1:length(code_mfa_dG)],'Color',GREEN_COLOR, 'LineWidth',4);
hold on;
plot(code_mfa_dG_per_wc_measurements',[1:1:length(code_mfa_dG_per_wc_measurements);1:1:length(code_mfa_dG_per_wc_measurements)],'k*','MarkerSize',12);
line([0 0], [0 length(code_mfa_dG_per_wc_measurements)+1], 'color','red','LineStyle',':');
ylim([0 length(code_mfa_dG_per_wc_measurements)+1]);
% xlim([-70 10]);
YLabel={''};
YLabel_right={''};
for(i=1:size(rxns_text,1))
    YLabel{end+1}=rxns_text{i};
    YLabel{end} = strrep(YLabel{end},'_',' ');
    YLabel_right{end+1} = enzyme_text{i};
end
YLabel{end+1}='';
set(gca, 'YTickLabel', YLabel);  
set(gca, 'yTick', [0:length(YLabel)+1]);
xlabel('Gibbs free energy [kJ/mol]');
title('Mitochondial reactions');
set(gcf,'color','w');
ax = gca;
ax.Box='on';
ax.BoxStyle = 'full';
grid;
xlim([-40 40]);
xticks([-40:10:40]);
set(gca, 'FontSize',10);
yyaxis right;
ylim([0 length(code_mfa_dG_per_wc_measurements)+1]);
set(gca, 'YTickLabel', YLabel_right,'YColor','black');
set(gca, 'yTick', [0:length(YLabel)+1]);
hold off;


% plot dG of transporters
tTransporters = tTransporters(2:end,:);
code_mfa_dG = [];
for(i=1:size(tTransporters,1))
    index_rxns = find(contains(model_thermodynamics.full_rxns,tTransporters(i,3)));
    
    % change the net flux of the reaction to be positive
    direction_of_reaction = 1;
    if(net_fluxes{index_rxns(1)}.high_flux < 0)
        direction_of_reaction = 0;
    end
    
    
    if(model_thermodynamics.thermodynamics_of_reaction_defined(index_rxns))
        code_mfa_dG(i,:) = length(index_rxns)*[sensitiviy_analysis_dG{index_rxns(1)}.low_dG sensitiviy_analysis_dG{index_rxns(1)}.high_dG];
    else
        code_mfa_dG(i,:) = [nan nan];
    end        
      
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
    
    if(direction_of_reaction == 0)
        code_mfa_dG(i,[1 2]) = -code_mfa_dG(i,[2 1]);
        code_mfa_dG_per_wc_measurements(i) = -code_mfa_dG_per_wc_measurements(i);
        str_of_substrates_and_producs_by_reaction = split(tTransporters(i,1),' => ');
        tTransporters(i,1)={[str_of_substrates_and_producs_by_reaction{2} ' => ' str_of_substrates_and_producs_by_reaction{1}]};
    end    
end

[val ind_sort_transporters]=sort(code_mfa_dG(:,1));
code_mfa_dG = code_mfa_dG(ind_sort_transporters,:);
rxns_text = tTransporters(ind_sort_transporters,1)
enzyme_text = tTransporters(ind_sort_transporters,2);
code_mfa_dG_per_wc_measurements = code_mfa_dG_per_wc_measurements(ind_sort_transporters);



subplot(4,1,4);
plot(code_mfa_dG',[1:length(code_mfa_dG);1:length(code_mfa_dG)],'Color',RED_COLOR, 'LineWidth',4);
hold on;
plot(code_mfa_dG_per_wc_measurements',[1:1:length(code_mfa_dG_per_wc_measurements);1:1:length(code_mfa_dG_per_wc_measurements)],'k*','MarkerSize',12);
line([0 0], [0 length(code_mfa_dG_per_wc_measurements)+1], 'color','red','LineStyle',':');
ylim([0 length(code_mfa_dG_per_wc_measurements)+1]);
% xlim([-70 10]);
YLabel={''};
YLabel_right={''};
for(i=1:size(rxns_text,1))
    YLabel{end+1}=rxns_text{i};
    YLabel{end} = strrep(YLabel{end},'_',' ');
    YLabel_right{end+1} = enzyme_text{i};
end
YLabel{end+1}='';
set(gca, 'YTickLabel', YLabel);  
set(gca, 'yTick', [0:length(YLabel)+1]);
xlabel('Gibbs free energy [kJ/mol]');
title('Mitochondrial transporters');
set(gcf,'color','w');
ax = gca;
ax.Box='on';
ax.BoxStyle = 'full';
grid;
xlim([-40 40]);
xticks([-40:10:40]);
set(gca, 'FontSize',10);
yyaxis right;
ylim([0 length(code_mfa_dG_per_wc_measurements)+1]);
set(gca, 'YTickLabel', YLabel_right,'YColor','black');
set(gca, 'yTick', [0:length(YLabel)+1]);
hold off;

set(gcf,'PaperSize',[50 30]);
s = sprintf('./output_images/4abc-gibbs.pdf');
saveas(gcf, s);


% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% Code-MFA fluxes
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% Take flux of code-MFA result for CY and for MT
for(i=1:size(tBothCompartments,1))
    index_rxns_cy = find(contains(model_thermodynamics.full_rxns,tBothCompartments(i,3)));
    index_rxns_mt = find(contains(model_thermodynamics.full_rxns,tBothCompartments(i,4)));
     
    % check if name in figure is the same or opposite direction to the
    % model
        if(dBothCompartments(i,1))
            %code_mfa_flux_cy(i,:) = length(index_rxns_cy)*[sensitiviy_analysis_dG{index_rxns_cy(1)}.low_dG sensitiviy_analysis_dG{index_rxns_cy(1)}.high_dG];
            code_mfa_flux_cy(i,:) = length(index_rxns_cy)*[net_fluxes{index_rxns_cy(1)}.low_flux net_fluxes{index_rxns_cy(1)}.high_flux];
        else
            %code_mfa_flux_cy(i,:) = length(index_rxns_cy)*[-sensitiviy_analysis_dG{index_rxns_cy(1)}.high_dG -sensitiviy_analysis_dG{index_rxns_cy(1)}.low_dG];
            code_mfa_flux_cy(i,:) = length(index_rxns_cy)*[-net_fluxes{index_rxns_cy(1)}.high_flux -net_fluxes{index_rxns_cy(1)}.low_flux];
        end
    
        if(dBothCompartments(i,2))
%             code_mfa_dG_mt(i,:) = length(index_rxns_mt)*[sensitiviy_analysis_dG{index_rxns_mt(1)}.low_dG sensitiviy_analysis_dG{index_rxns_mt(1)}.high_dG];
            code_mfa_flux_mt(i,:) = length(index_rxns_mt)*[net_fluxes{index_rxns_mt(1)}.low_flux net_fluxes{index_rxns_mt(1)}.high_flux];
        else
%             code_mfa_dG_mt(i,:) = length(index_rxns_mt)*[-sensitiviy_analysis_dG{index_rxns_mt(1)}.high_dG -sensitiviy_analysis_dG{index_rxns_mt(1)}.low_dG];
            code_mfa_flux_mt(i,:) = length(index_rxns_mt)*[-net_fluxes{index_rxns_mt(1)}.high_flux -net_fluxes{index_rxns_mt(1)}.low_flux];
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
end

code_mfa_flux_cy = code_mfa_flux_cy(ind_sort_cy_mt_reactions,:);
code_mfa_flux_mt = code_mfa_flux_mt(ind_sort_cy_mt_reactions,:);
rxns_text_for_both_compartments = tBothCompartments(ind_sort_cy_mt_reactions,1)
enzyme_text = tBothCompartments(ind_sort_cy_mt_reactions,2);

figure('Position',[1 1 800 1000]);
CONST_VAL_FOR_PLOT_MT_ABOVE_CY=0.08;
subplot(4,1,1);
h1=plot(code_mfa_flux_cy',[1-CONST_VAL_FOR_PLOT_MT_ABOVE_CY:1:length(code_mfa_flux_cy);1-CONST_VAL_FOR_PLOT_MT_ABOVE_CY:1:length(code_mfa_flux_cy)],'Color',BLUE_COLOR, 'LineWidth',4);
hold on;
h2=plot(code_mfa_flux_mt',[1+CONST_VAL_FOR_PLOT_MT_ABOVE_CY:1:length(code_mfa_flux_mt)+CONST_VAL_FOR_PLOT_MT_ABOVE_CY;1+CONST_VAL_FOR_PLOT_MT_ABOVE_CY:1:length(code_mfa_flux_mt)+CONST_VAL_FOR_PLOT_MT_ABOVE_CY],'Color',GREEN_COLOR, 'LineWidth',4);
h3_dummy=plot(100,100,'Color',RED_COLOR,'LineWidth',4);
line([0 0], [0 length(rxns_text_for_both_compartments)+1], 'color','red','LineStyle',':');
ylim([0 size(rxns_text_for_both_compartments,1)+1]);
% xlim([-70 10]);
YLabel={''};
YLabel_right={''};
for(i=1:size(rxns_text_for_both_compartments,1))
    YLabel{end+1}=rxns_text_for_both_compartments{i};
    YLabel{end} = strrep(YLabel{end},'_',' ');
    YLabel_right{end+1} = enzyme_text{i};
end
YLabel{end+1}=''; 
set(gca, 'YTickLabel', YLabel);  
set(gca, 'yTick', [0:length(YLabel)+1]);
xlabel('Net flux [mM/h]');
title('Cytosolic & mitochondrial reactions');
set(gcf,'color','w');
ax = gca;
ax.Box='on';
ax.BoxStyle = 'full';
grid;
% legend([h1(1),h2(1),h4_dummy(1),h3(1)],'Cytosolic reactions','Mitochondrial reactions','Transporrters','Measured (WC-level)','FontSize',18,'Orientation','horizontal');
% legend([h1(1),h2(1),h3_dummy(1),h3(1)],'Cytosolic','Mitochondrial','Transporters','Measured (WC)','FontSize',12, 'Orientation','horizontal', 'Location', 'north');
xlim([-60 60]);
xticks([-100:20:100]);
set(gca, 'FontSize',10);
yyaxis right;
ylim([0 size(rxns_text_for_both_compartments,1)+1]);
set(gca, 'YTickLabel', YLabel_right,'YColor','black');
set(gca, 'yTick', [0:length(YLabel)+1]);
hold off;

 


% plot flux of cytosolic reactions
code_mfa_flux = [];
for(i=1:size(tCytosolic,1))
    index_rxns = find(contains(model_thermodynamics.full_rxns,tCytosolic(i,3)));
    
    % change the net flux of the reaction to be positive
    direction_of_reaction = 1;
    if(net_fluxes{index_rxns(1)}.high_flux < 0)
        direction_of_reaction = 0;
    end
    
    code_mfa_flux(i,:) = length(index_rxns)*[net_fluxes{index_rxns(1)}.low_flux net_fluxes{index_rxns(1)}.high_flux];
          
    if(direction_of_reaction == 0)
        code_mfa_flux(i,[1 2]) = -code_mfa_flux(i,[2 1]);
    end    
end

code_mfa_flux = code_mfa_flux(ind_sort_cy_reactions,:);
rxns_text = tCytosolic(ind_sort_cy_reactions,1)
enzyme_text = tCytosolic(ind_sort_cy_reactions,2);


subplot(4,1,2);
plot(code_mfa_flux',[1:length(code_mfa_flux);1:length(code_mfa_flux)],'Color',BLUE_COLOR, 'LineWidth',4);
hold on;
line([0 0], [0 length(rxns_text)+1], 'color','red','LineStyle',':');
ylim([0 length(rxns_text)+1]);
% xlim([-70 10]);
YLabel={''};
YLabel_right={''};
for(i=1:size(rxns_text,1))
    YLabel{end+1}=rxns_text{i};
    YLabel{end} = strrep(YLabel{end},'_',' ');
    YLabel_right{end+1} = enzyme_text{i};
end
YLabel{end+1}='';
set(gca, 'YTickLabel', YLabel);  
set(gca, 'yTick', [0:length(YLabel)+1]);
xlabel('Net flux [mM/h]');
title('Cytosolic reactions');
set(gcf,'color','w');
ax = gca;
ax.Box='on';
ax.BoxStyle = 'full';
grid;
xlim([-60 60]);
xticks([-100:20:100]);
set(gca, 'FontSize',10);
yyaxis right;
ylim([0 length(rxns_text)+1]);
set(gca, 'YTickLabel', YLabel_right,'YColor','black');
set(gca, 'yTick', [0:length(YLabel)+1]);
hold off;



% plot flux of mitochondrial reactions
code_mfa_flux = [];
for(i=1:size(tMitochondrial,1))
    index_rxns = find(contains(model_thermodynamics.full_rxns,tMitochondrial(i,3)));
    
    % change the net flux of the reaction to be positive
    direction_of_reaction = 1;
    if(net_fluxes{index_rxns(1)}.high_flux < 0)
        direction_of_reaction = 0;
    end
    
    code_mfa_flux(i,:) = length(index_rxns)*[net_fluxes{index_rxns(1)}.low_flux net_fluxes{index_rxns(1)}.high_flux];
    
    if(direction_of_reaction == 0)
        code_mfa_flux(i,[1 2]) = -code_mfa_flux(i,[2 1]);
    end    
     
end

code_mfa_flux = code_mfa_flux(ind_sort_mt_reactions,:);
rxns_text = tMitochondrial(ind_sort_mt_reactions,1)
enzyme_text = tMitochondrial(ind_sort_mt_reactions,2);


subplot(4,1,3);
plot(code_mfa_flux',[1:length(code_mfa_flux);1:length(code_mfa_flux)],'Color',GREEN_COLOR, 'LineWidth',4);
hold on;
line([0 0], [0 length(rxns_text)+1], 'color','red','LineStyle',':');
ylim([0 length(rxns_text)+1]);
% xlim([-70 10]);
YLabel={''};
YLabel_right={''};
for(i=1:size(rxns_text,1))
    YLabel{end+1}=rxns_text{i};
    YLabel{end} = strrep(YLabel{end},'_',' ');
    YLabel_right{end+1} = enzyme_text{i};
end
YLabel{end+1}='';
set(gca, 'YTickLabel', YLabel);  
set(gca, 'yTick', [0:length(YLabel)+1]);
xlabel('Net flux [mM/h]');
title('Mitochondial reactions');
set(gcf,'color','w');
ax = gca;
ax.Box='on';
ax.BoxStyle = 'full';
grid;
xlim([-60 60]);
xticks([-100:20:100]);
set(gca, 'FontSize',10);
yyaxis right;
ylim([0 length(rxns_text)+1]);
set(gca, 'YTickLabel', YLabel_right,'YColor','black');
set(gca, 'yTick', [0:length(YLabel)+1]);
hold off;



% plot flux of transporters
code_mfa_flux = [];
for(i=1:size(tTransporters,1))
    index_rxns = find(contains(model_thermodynamics.full_rxns,tTransporters(i,3)));
    
    % change the net flux of the reaction to be positive
    direction_of_reaction = 1;
    if(net_fluxes{index_rxns(1)}.high_flux < 0)
        direction_of_reaction = 0;
    end    
    
    code_mfa_flux(i,:) = length(index_rxns)*[net_fluxes{index_rxns(1)}.low_flux net_fluxes{index_rxns(1)}.high_flux];
    
    if(direction_of_reaction == 0)
        code_mfa_flux(i,[1 2]) = -code_mfa_flux(i,[2 1]);
    end    
    
end

code_mfa_flux = code_mfa_flux(ind_sort_transporters,:);
rxns_text = tTransporters(ind_sort_transporters,1)
enzyme_text = tTransporters(ind_sort_transporters,2);



subplot(4,1,4);
plot(code_mfa_flux',[1:length(code_mfa_flux);1:length(code_mfa_flux)],'Color',RED_COLOR, 'LineWidth',4);
hold on;
line([0 0], [0 length(rxns_text)+1], 'color','red','LineStyle',':');
ylim([0 length(rxns_text)+1]);
% xlim([-70 10]);
YLabel={''};
YLabel_right={''};
for(i=1:size(rxns_text,1))
    YLabel{end+1}=rxns_text{i};
    YLabel{end} = strrep(YLabel{end},'_',' ');
    YLabel_right{end+1} = enzyme_text{i};
end
YLabel{end+1}='';
set(gca, 'YTickLabel', YLabel);  
set(gca, 'yTick', [0:length(YLabel)+1]);
xlabel('Net flux [mM/h]');
title('Mitochondrial transporters');
set(gcf,'color','w');
ax = gca;
ax.Box='on';
ax.BoxStyle = 'full';
grid;
xlim([-60 60]);
xticks([-100:20:100]);
set(gca, 'FontSize',10);
yyaxis right;
ylim([0 length(rxns_text)+1]);
set(gca, 'YTickLabel', YLabel_right,'YColor','black');
set(gca, 'yTick', [0:length(YLabel)+1]);
hold off;


set(gcf,'PaperSize',[50 30]);
s = sprintf('./output_images/4abc-fluxes.pdf');
saveas(gcf, s);


% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% plot flux of mitochondrial reactions
code_mfa_flux = [];

mitochondrial_reaction{1,1} = 'Malate_MT + NAD_MT => Pyruvate_MT + CO2_sink + NADH_MT';
mitochondrial_reaction{2,1} = 'AKG_MT + CO2_source + NADH_MT => Citrate_MT + NAD_MT';
mitochondrial_reaction{3,1} = 'Pyruvate_MT + NAD_MT + CoA_MT => CO2_sink + Acetyl_CoA_MT + NADH_MT';
mitochondrial_reaction{4,1} = 'Malate_MT + NAD_MT => OAA_MT + NADH_MT';
mitochondrial_reaction{5,1} = 'AKG_MT => CO2_sink + Fumarate_MT';
for(i=1:size(mitochondrial_reaction,1))
    index_rxns = find(contains(model_thermodynamics.full_rxns,mitochondrial_reaction{i}));
%    mitochondrial_reaction{end+1,1} = model_thermodynamics.full_rxns{index_rxns};
    
    % change the net flux of the reaction to be positive
    direction_of_reaction = 1;
    if(net_fluxes{index_rxns(1)}.high_flux < 0)
        direction_of_reaction = 0;
    end
    
    code_mfa_flux(i,:) = length(index_rxns)*[net_fluxes{index_rxns(1)}.low_flux net_fluxes{index_rxns(1)}.high_flux];
    
    if(direction_of_reaction == 0)
        code_mfa_flux(i,[1 2]) = -code_mfa_flux(i,[2 1]);
        
        str_of_substrates_and_producs_by_reaction = split(mitochondrial_reaction{i,1},' => ');
        mitochondrial_reaction{i,1}=[str_of_substrates_and_producs_by_reaction{2} ' => ' str_of_substrates_and_producs_by_reaction{1}];
        
    end    
     
end

%code_mfa_flux = code_mfa_flux(ind_sort_mt_reactions,:);
rxns_text = mitochondrial_reaction;
%enzyme_text = tMitochondrial(:,2);

figure('Position',[1 1 800 340]);
plot(code_mfa_flux',[1:length(code_mfa_flux);1:length(code_mfa_flux)],'Color',GREEN_COLOR, 'LineWidth',4);
hold on;
line([0 0], [0 length(rxns_text)+1], 'color','red','LineStyle',':');
ylim([0 length(rxns_text)+1]);
% xlim([-70 10]);
YLabel={''};
YLabel_right={''};
for(i=1:size(rxns_text,1))
    YLabel{end+1}=rxns_text{i};
    YLabel{end} = strrep(YLabel{end},'_',' ');
%    YLabel_right{end+1} = enzyme_text{i};
end
YLabel{end+1}='';
set(gca, 'YTickLabel', YLabel);  
set(gca, 'yTick', [0:length(YLabel)+1]);
xlabel('Net flux [mM/h]');
title('Mitochondial reactions producint NADH');
set(gcf,'color','w');
ax = gca;
ax.Box='on';
ax.BoxStyle = 'full';
grid;
xlim([-60 60]);
xticks([-100:20:100]);
set(gca, 'FontSize',10);
yyaxis right;
ylim([0 length(rxns_text)+1]);
set(gca, 'YTickLabel', YLabel_right,'YColor','black');
set(gca, 'yTick', [0:length(YLabel)+1]);
hold off;

