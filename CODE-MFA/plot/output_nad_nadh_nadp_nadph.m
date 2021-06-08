% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% The measured NADP+/NADPH ratio and NADP+ and NADPH concentrations (gray) 
% and inferred compartmentalized ratio and concentrations based on simple 
% thermodynamics and deconvolution analysis (blue), and via CODE-MFA (red). 
% Asterisks represent previously published ratios (Table S12). 
% The measured NAD+/NADH ratio and NAD+ and NADH concentrations (gray) 
% and inferred compartmentalized values based on simple thermodynamics and 
% deconvolution analysis (blue), and via CODE-MFA (red)
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 

clear all;
close all;

run ../load_constants;

% insert meassured values for this cell line
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Lactate_m3_from_lactate_labeling_experiment  = 0.025;
Pyruvate_m3_from_lactate_labeling_experiment = 0.01;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% plot all papers together or saparetely in legend
PLOT_ALL_PAPERS = 0;

% load model and sensitivity analysis results from Code-MFA
load('../mat_files/sensitiviy_analysis_dG.mat', 'sensitiviy_analysis_dG');
load('../mat_files/model_thermodynamics.mat','model_thermodynamics');    
load('../mat_files/sensitiviy_analysis_concentration.mat', 'sensitiviy_analysis_concentration');
load('../mat_files/sensitiviy_analysis_cofactors_ratio.mat','sensitiviy_analysis_cofactors_ratio');
load('../mat_files/WC_known_metabolites_idv.mat','WC_known_metabolites_idv');

all_met_WC_measured_con = [model_thermodynamics.WC.Concentrations model_thermodynamics.WC.Concentrations_STD];
all_met_names_with_WC_measured_con = model_thermodynamics.WC.met_name;

all_met_names_with_WC_measured_con{end+1} = '6phosphogluconate';
index_compartmentalized_con = find(ismember(model_thermodynamics.mets,'6phosphogluconate_CY'))
all_met_WC_measured_con(end+1,:) = [(exp(model_thermodynamics.mets_lb(index_compartmentalized_con))+exp(model_thermodynamics.mets_ub(index_compartmentalized_con)))/2 (exp(model_thermodynamics.mets_ub(index_compartmentalized_con))-exp(model_thermodynamics.mets_lb(index_compartmentalized_con)))/4];
all_met_names_with_WC_measured_con{end+1} = 'Glc6P';
index_compartmentalized_con = find(ismember(model_thermodynamics.mets,'Glc6P_CY'))
all_met_WC_measured_con(end+1,:) = [(exp(model_thermodynamics.mets_lb(index_compartmentalized_con))+exp(model_thermodynamics.mets_ub(index_compartmentalized_con)))/2 (exp(model_thermodynamics.mets_ub(index_compartmentalized_con))-exp(model_thermodynamics.mets_lb(index_compartmentalized_con)))/4];
all_met_names_with_WC_measured_con{end+1} = 'Ribose5P';
index_compartmentalized_con = find(ismember(model_thermodynamics.mets,'Ribose5P_CY'))
all_met_WC_measured_con(end+1,:) = [(exp(model_thermodynamics.mets_lb(index_compartmentalized_con))+exp(model_thermodynamics.mets_ub(index_compartmentalized_con)))/2 (exp(model_thermodynamics.mets_ub(index_compartmentalized_con))-exp(model_thermodynamics.mets_lb(index_compartmentalized_con)))/4];
all_met_names_with_WC_measured_con{end+1} = 'G3P';
index_compartmentalized_con = find(ismember(model_thermodynamics.mets,'G3P'))
all_met_WC_measured_con(end+1,:) = [(exp(model_thermodynamics.mets_lb(index_compartmentalized_con))+exp(model_thermodynamics.mets_ub(index_compartmentalized_con)))/2 (exp(model_thermodynamics.mets_ub(index_compartmentalized_con))-exp(model_thermodynamics.mets_lb(index_compartmentalized_con)))/4];
all_met_names_with_WC_measured_con{end+1} = '13BPG';
index_compartmentalized_con = find(ismember(model_thermodynamics.mets,'13BPG'))
all_met_WC_measured_con(end+1,:) = [(exp(model_thermodynamics.mets_lb(index_compartmentalized_con))+exp(model_thermodynamics.mets_ub(index_compartmentalized_con)))/2 (exp(model_thermodynamics.mets_ub(index_compartmentalized_con))-exp(model_thermodynamics.mets_lb(index_compartmentalized_con)))/4];
all_met_names_with_WC_measured_con{end+1} = '3PG';
index_compartmentalized_con = find(ismember(model_thermodynamics.mets,'3PG'))
all_met_WC_measured_con(end+1,:) = [(exp(model_thermodynamics.mets_lb(index_compartmentalized_con))+exp(model_thermodynamics.mets_ub(index_compartmentalized_con)))/2 (exp(model_thermodynamics.mets_ub(index_compartmentalized_con))-exp(model_thermodynamics.mets_lb(index_compartmentalized_con)))/4];



MINIMUM_CONCENTRATION = 1e-5; %in mM
RT = 2.478;
cyto_volume = 0.8;
mito_volume = 1-cyto_volume;


for(i=1:length(sensitiviy_analysis_cofactors_ratio))
    if(strcmp(sensitiviy_analysis_cofactors_ratio{i}.co_factors_name,'NAD_CY/NADH_CY'))
        code_mfa_cy_nad_nadh_ratio_low  = sensitiviy_analysis_cofactors_ratio{i}.low_cofactor_ratio;
        code_mfa_cy_nad_nadh_ratio_high = sensitiviy_analysis_cofactors_ratio{i}.high_cofactor_ratio;
    elseif(strcmp(sensitiviy_analysis_cofactors_ratio{i}.co_factors_name,'NAD_MT/NADH_MT'))
        code_mfa_mt_nad_nadh_ratio_low  = sensitiviy_analysis_cofactors_ratio{i}.low_cofactor_ratio;
        code_mfa_mt_nad_nadh_ratio_high = sensitiviy_analysis_cofactors_ratio{i}.high_cofactor_ratio;        
    elseif(strcmp(sensitiviy_analysis_cofactors_ratio{i}.co_factors_name,'NADP_CY/NADPH_CY'))
        code_mfa_cy_nadp_nadph_ratio_low  = sensitiviy_analysis_cofactors_ratio{i}.low_cofactor_ratio;
        code_mfa_cy_nadp_nadph_ratio_high = sensitiviy_analysis_cofactors_ratio{i}.high_cofactor_ratio;                
    elseif(strcmp(sensitiviy_analysis_cofactors_ratio{i}.co_factors_name,'NADP_MT/NADPH_MT'))
        code_mfa_mt_nadp_nadph_ratio_low  = sensitiviy_analysis_cofactors_ratio{i}.low_cofactor_ratio;
        code_mfa_mt_nadp_nadph_ratio_high = sensitiviy_analysis_cofactors_ratio{i}.high_cofactor_ratio;                        
    end
end
for(i=1:length(sensitiviy_analysis_concentration))
    if(strcmp(sensitiviy_analysis_concentration{i}.metabolite_name,'NAD_CY'))
        code_mfa_cy_nad_low      = exp(sensitiviy_analysis_concentration{i}.low_concentration);
        code_mfa_cy_nad_high     = exp(sensitiviy_analysis_concentration{i}.high_concentration);
    elseif(strcmp(sensitiviy_analysis_concentration{i}.metabolite_name,'NADH_CY'))
        code_mfa_cy_nadh_low     = exp(sensitiviy_analysis_concentration{i}.low_concentration);
        code_mfa_cy_nadh_high    = exp(sensitiviy_analysis_concentration{i}.high_concentration);
    elseif(strcmp(sensitiviy_analysis_concentration{i}.metabolite_name,'NAD_MT'))
        code_mfa_mt_nad_low      = exp(sensitiviy_analysis_concentration{i}.low_concentration);
        code_mfa_mt_nad_high     = exp(sensitiviy_analysis_concentration{i}.high_concentration);
    elseif(strcmp(sensitiviy_analysis_concentration{i}.metabolite_name,'NADH_MT'))
        code_mfa_mt_nadh_low     = exp(sensitiviy_analysis_concentration{i}.low_concentration);        
        code_mfa_mt_nadh_high    = exp(sensitiviy_analysis_concentration{i}.high_concentration);        
    elseif(strcmp(sensitiviy_analysis_concentration{i}.metabolite_name,'NADP_CY'))
        code_mfa_cy_nadp_low      = exp(sensitiviy_analysis_concentration{i}.low_concentration);
        code_mfa_cy_nadp_high     = exp(sensitiviy_analysis_concentration{i}.high_concentration);
    elseif(strcmp(sensitiviy_analysis_concentration{i}.metabolite_name,'NADPH_CY'))
        code_mfa_cy_nadph_low     = exp(sensitiviy_analysis_concentration{i}.low_concentration);
        code_mfa_cy_nadph_high    = exp(sensitiviy_analysis_concentration{i}.high_concentration);
    elseif(strcmp(sensitiviy_analysis_concentration{i}.metabolite_name,'NADP_MT'))
        code_mfa_mt_nadp_low      = exp(sensitiviy_analysis_concentration{i}.low_concentration);
        code_mfa_mt_nadp_high     = exp(sensitiviy_analysis_concentration{i}.high_concentration);
    elseif(strcmp(sensitiviy_analysis_concentration{i}.metabolite_name,'NADPH_MT'))
        code_mfa_mt_nadph_low     = exp(sensitiviy_analysis_concentration{i}.low_concentration);        
        code_mfa_mt_nadph_high    = exp(sensitiviy_analysis_concentration{i}.high_concentration);        
    end    
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% NAD/NADH
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Ratio from literature - Mitochondria
% Lee, W. D., Mukha, D., Aizenshtein, E., & Shlomi, T. (2019). Spatial-fluxomics provides a subcellular-compartmentalized view of reductive glutamine metabolism in cancer cells. Nature Communications. https://doi.org/10.1038/s41467-019-09352-1
mt_NAD_NADH_ratio_won = 1/0.015;
% Chen, W. W., Freinkman, E., Wang, T., Birsoy, K., & Sabatini, D. M. (2016). Absolute Quantification of Matrix Metabolites Reveals the Dynamics of Mitochondrial Metabolism. Cell. https://doi.org/10.1016/j.cell.2016.07.040
mt_NAD_NADH_ratio_sabatini = 1/0.009;
% Krebs, H. A. (1967). The redox state of nicotinamide adenine dinucleotide in the cytoplasm and mitochondria of rat liver. Advances in Enzyme Regulation. https://doi.org/10.1016/0065-2571(67)90029-5
mt_NAD_NADH_ratio_krebs = 8;
% Siess, E. A., Brocks, D. G., Lattke, H. K., & Wieland, O. H. (1977). Effect of glucagon on metabolite compartmentation in isolated rat liver cells during gluconeogenesis from lactate. Biochemical Journal. https://doi.org/10.1042/bj1660225
mt_NAD_NADH_ratio_elmar = 60;


% Ratio from literature - Cytosol
% Chen, W. W., Freinkman, E., Wang, T., Birsoy, K., & Sabatini, D. M. (2016). Absolute Quantification of Matrix Metabolites Reveals the Dynamics of Mitochondrial Metabolism. Cell. https://doi.org/10.1016/j.cell.2016.07.040
cy_NAD_NADH_ratio_sabatini = 1/0.15;
% Krebs, H. A. (1967). The redox state of nicotinamide adenine dinucleotide in the cytoplasm and mitochondria of rat liver. Advances in Enzyme Regulation. https://doi.org/10.1016/0065-2571(67)90029-5
cy_NAD_NADH_ratio_krebs = 725;
% Sun, F., Dai, C., Xie, J., & Hu, X. (2012). Biochemical issues in estimation of cytosolic free NAD/NADH ratio. PLoS ONE. https://doi.org/10.1371/journal.pone.0034525
cy_NAD_NADH_ratio_sun = 88;
% Hedeskov, C. J., Capito, K., & Thams, P. (1987). Cytosolic ratios of free [NADPH]/[NADP+] and [NADH]/[NAD+] in mouse pancreatic islets, and nutrient-induced insulin secretion. The Biochemical Journal, 241(1), 161–167. https://doi.org/10.1042/bj2410161
cy_NAD_NADH_ratio_hedeskov = 1000;
% Zhao, Y., Hu, Q., Cheng, F., Su, N., Wang, A., Zou, Y., … Yang, Y. (2015). SoNar, a Highly Responsive NAD+/NADH Sensor, Allows High-Throughput Metabolic Screening of Anti-tumor Agents. Cell Metabolism. https://doi.org/10.1016/j.cmet.2015.04.009
cy_NAD_NADH_ratio_zhao = 800;

% Metabolite concentrations [mM]
index_of_metabolite_in_wc_measured_vector = find(ismember(all_met_names_with_WC_measured_con,'NAD'));
NAD_con = all_met_WC_measured_con(index_of_metabolite_in_wc_measured_vector,:);
NAD_con_min = NAD_con(1)-NAD_con(2)*2;
NAD_con_max = NAD_con(1)+NAD_con(2)*2;

index_of_metabolite_in_wc_measured_vector = find(ismember(all_met_names_with_WC_measured_con,'NADH'));
NADH_con = all_met_WC_measured_con(index_of_metabolite_in_wc_measured_vector,:);
NADH_con_min = NADH_con(1)-NADH_con(2)*2;
NADH_con_max = NADH_con(1)+NADH_con(2)*2;

mito_NAD_con_max    = NAD_con_max/mito_volume;
mito_NADH_con_max   = NADH_con_max/mito_volume;
cyto_NAD_con_max    = NAD_con_max/cyto_volume;
cyto_NADH_con_max   = NADH_con_max/cyto_volume;

index_of_metabolite_in_wc_measured_vector = find(ismember(all_met_names_with_WC_measured_con,'Pyruvate'));
Pyruvate_con = all_met_WC_measured_con(index_of_metabolite_in_wc_measured_vector,:);
index_of_metabolite_in_wc_measured_vector = find(ismember(all_met_names_with_WC_measured_con,'Lactate'));
Lactate_con = all_met_WC_measured_con(index_of_metabolite_in_wc_measured_vector,:);
G0_Pyruvate_mt_to_cy_transporter    = [2.9 0.29];
G0_Pyruvate_to_Lactate = [-27.5 0.45];
% G_Pyruvate_to_Lactate = -2; % LDH is close to chemical equilibrium, based on Lactate labeling measurements

Pyruvate_to_Lactate_forward_div_backward_flux = Lactate_m3_from_lactate_labeling_experiment/Pyruvate_m3_from_lactate_labeling_experiment;
G_Pyruvate_to_Lactate = -RT*log(Pyruvate_to_Lactate_forward_div_backward_flux(1));


% solve a set of linear equations to find lower bound to cytosolic pyruvate
syms Pyruvate_cy_syms_con Pyruvate_mt_syms_con
% assuming net flux from cy to my
% -G0_Pyruvate_mt_to_cy_transporter+RT*log(Pyruvate_mt_syms_con/Pyruvate_cy_syms_con)<=0 
Pyruvate_cy_con_vector_min = [];
for i=1:1000
    eqn1 = Pyruvate_mt_syms_con - Pyruvate_cy_syms_con*exp((G0_Pyruvate_mt_to_cy_transporter(1)+randn(1)*G0_Pyruvate_mt_to_cy_transporter(2))/RT) == 0;
    eqn2 = Pyruvate_mt_syms_con*mito_volume + Pyruvate_cy_syms_con*cyto_volume == Pyruvate_con(1)+randn(1)*Pyruvate_con(2);
    [A,B] = equationsToMatrix([eqn1, eqn2], [Pyruvate_cy_syms_con, Pyruvate_mt_syms_con]);
    solver_res = linsolve(A,B);
    Pyruvate_cy_con_vector_min(i) = double(solver_res(1));
end
Pyruvate_cy_con_min  = prctile(Pyruvate_cy_con_vector_min,5);
Pyruvate_cy_con_max = Pyruvate_con(1)+Pyruvate_con(2)*2;

% Pyruvate_CY + NADH_CY => Lactate_CY + NAD_CY
n_Lactate_con                 = Lactate_con(1) + randn(1,1000) * Lactate_con(2);
n_Pyruvate_cy_con             = normrnd((Pyruvate_cy_con_min+Pyruvate_cy_con_max)/2,(Pyruvate_cy_con_max+Pyruvate_cy_con_min)/4,1,1000);
n_Pyruvate_cy_con(n_Pyruvate_cy_con<1e-5)=1e-5;
n_G0_Pyruvate_to_Lactate = G0_Pyruvate_to_Lactate(1) + randn(1,1000) * G0_Pyruvate_to_Lactate(2);

cyto_NAD_NADH_ratio = exp((G_Pyruvate_to_Lactate-n_G0_Pyruvate_to_Lactate)/RT - log(n_Lactate_con) + log(n_Pyruvate_cy_con));
cyto_NADH_NAD_ratio = 1./cyto_NAD_NADH_ratio;


figure('units','normalized','outerposition',[0 0 1 1])


% figure for NAD/NADH ratio
subplot(1,2,1);
grid on;
hold on;

% WC measured
wc_nad_nadh_ratio_low   = NAD_con_min/NADH_con_max;
wc_nad_nadh_ratio_high  = NAD_con_max/NADH_con_min;

% Cytosolic, based on G0 of Pyruvate transporter, L/P ratio, and LDH being
% in chemical equilibrium
cy_nadh_nad_ratio_low  = prctile(cyto_NADH_NAD_ratio,5);
cy_nadh_nad_ratio_high = prctile(cyto_NADH_NAD_ratio,95);
cy_nad_nadh_ratio_low   = 1/cy_nadh_nad_ratio_high;
cy_nad_nadh_ratio_high  = 1/cy_nadh_nad_ratio_low;
  
 

% find min and max for cytosolic NAD and NADH
% linear programming [cyto_NAD cyto_NADH mito_NAD mito_NADH]
lb = [MINIMUM_CONCENTRATION MINIMUM_CONCENTRATION MINIMUM_CONCENTRATION MINIMUM_CONCENTRATION];
ub = [cyto_NAD_con_max cyto_NADH_con_max mito_NAD_con_max mito_NADH_con_max];
% use cytosolic NAD/NADH ratio, found by Pyruvate transporter, LDH, and L/P ratio
A1=[-1 cy_nad_nadh_ratio_low 0 0; 1 -cy_nad_nadh_ratio_high 0 0];
b1=[0;0];
A2=[cyto_volume 0 mito_volume 0; -cyto_volume 0 -mito_volume 0; 0 cyto_volume 0 mito_volume; 0 -cyto_volume 0 -mito_volume];
b2=[NAD_con_max; -NAD_con_min; NADH_con_max; -NADH_con_min];
% assuming mallate aspartate shuttle is working in the canonical way

A=[A1;A2];
b=[b1;b2];


f = [1 0 0 0];  % min cyto NAD
[initial_values_new,fval,exitflag] = linprog(f, A, b, [], [], lb, ub);
cy_nad_low = fval;
f = [-1 0 0 0];  % max cyto NAD
[initial_values_new,fval,exitflag] = linprog(f, A, b, [], [], lb, ub);
cy_nad_high = -fval;
f = [0 1 0 0];  % min cyto NADH
[initial_values_new,fval,exitflag] = linprog(f, A, b, [], [], lb, ub);
cy_nadh_low = fval;
f = [0 -1 0 0];  % max cyto NADH
[initial_values_new,fval,exitflag] = linprog(f, A, b, [], [], lb, ub);
cy_nadh_high = -fval;
f = [0 0 1 0];  % min mito NAD
[initial_values_new,fval,exitflag] = linprog(f, A, b, [], [], lb, ub);
mt_nad_low = fval;
f = [0 0 -1 0];  % max mito NAD
[initial_values_new,fval,exitflag] = linprog(f, A, b, [], [], lb, ub);
mt_nad_high = -fval;
f = [0 0 0 1];  % min mito NADH
[initial_values_new,fval,exitflag] = linprog(f, A, b, [], [], lb, ub);
mt_nadh_low = fval;                        
f = [0 0 0 -1];  % max mito NADH
[initial_values_new,fval,exitflag] = linprog(f, A, b, [], [], lb, ub);
mt_nadh_high = -fval;
 
% check upper bound on mitochondrial NADP/NADPH ratio
Beq=0;
q_factor = 1.1; 
min_ratio = 1e-6;
max_ratio = 1e6;
mt_nad_nadh_ratio_low = min_ratio;
f = [0 0 0 0];  % min mito NADPH
while(mt_nad_nadh_ratio_low<max_ratio)
    Aeq=[0 0 1 -mt_nad_nadh_ratio_low];
    [initial_values_new,fval,exitflag] = linprog(f, A, b, Aeq, Beq, lb, ub);
    % no solution was found
    if(exitflag==1)
        break;
    end 
    mt_nad_nadh_ratio_low=mt_nad_nadh_ratio_low*q_factor;
end 
mt_nad_nadh_ratio_high=mt_nad_nadh_ratio_low;
while(mt_nad_nadh_ratio_high<max_ratio)        
    Aeq=[0 0 1 -mt_nad_nadh_ratio_high];
    [initial_values_new,fval,exitflag] = linprog(f, A, b, Aeq, Beq, lb, ub);
    % no solution was found
    if(exitflag~=1)
        break;
    end
    mt_nad_nadh_ratio_high=mt_nad_nadh_ratio_high*q_factor;
end


xlim([0.5 3.5]); 
  
w=0.25;
%p=1; plot(polyshape([p-w p-w p+w p+w], log10([wc_nad_nadh_ratio_low wc_nad_nadh_ratio_high wc_nad_nadh_ratio_high wc_nad_nadh_ratio_low])),'FaceColor', RED_COLOR, 'FaceAlpha',0.7);
p=1; plot(polyshape([p-w p-w p+w p+w], log10([wc_nad_nadh_ratio_low wc_nad_nadh_ratio_high wc_nad_nadh_ratio_high wc_nad_nadh_ratio_low])),'FaceColor', GRAY_COLOR, 'FaceAlpha',0.7);
p=2; plot(polyshape([p-w p-w p+w p+w], log10([cy_nad_nadh_ratio_low cy_nad_nadh_ratio_high cy_nad_nadh_ratio_high cy_nad_nadh_ratio_low])),'FaceColor', BLUE_COLOR, 'FaceAlpha',0.7);
%mean_val=(log10(code_mfa_cy_nad_nadh_ratio_low)+log10(code_mfa_cy_nad_nadh_ratio_high))/2; errorbar(p,mean_val,log10(code_mfa_cy_nad_nadh_ratio_high)-mean_val,log10(code_mfa_cy_nad_nadh_ratio_high)-mean_val,'color',BLUE_COLOR,'linewidth',3, 'LineStyle', 'none');
mean_val=(log10(code_mfa_cy_nad_nadh_ratio_low)+log10(code_mfa_cy_nad_nadh_ratio_high))/2; errorbar(p,mean_val,log10(code_mfa_cy_nad_nadh_ratio_high)-mean_val,log10(code_mfa_cy_nad_nadh_ratio_high)-mean_val,'color',RED_COLOR,'linewidth',3, 'LineStyle', 'none');
if(PLOT_ALL_PAPERS)
    plot(p,log10(cy_NAD_NADH_ratio_sabatini), 'g+', 'MarkerSize',10, 'LineWidth', 1.5);
    plot(p,log10(cy_NAD_NADH_ratio_krebs),'go', 'MarkerSize',10, 'LineWidth', 1.5);
    plot(p,log10(cy_NAD_NADH_ratio_sun),'gs', 'MarkerSize',10, 'LineWidth', 1.5);
    plot(p,log10(cy_NAD_NADH_ratio_hedeskov),'gd', 'MarkerSize',10, 'LineWidth', 1.5);
    plot(p,log10(cy_NAD_NADH_ratio_zhao),'g*', 'MarkerSize',10, 'LineWidth', 1.5);
else
    all_literature_vals = [cy_NAD_NADH_ratio_sabatini cy_NAD_NADH_ratio_krebs cy_NAD_NADH_ratio_sun cy_NAD_NADH_ratio_hedeskov cy_NAD_NADH_ratio_zhao];
    %plot(p*ones(1,length(all_literature_vals)),log10(all_literature_vals), 'MarkerSize',12, 'LineStyle', 'none', 'Color', BLUE_COLOR, 'Marker', '*');
    plot(p*ones(1,length(all_literature_vals)),log10(all_literature_vals), 'MarkerSize',12, 'LineWidth', 1.6, 'LineStyle', 'none', 'Color', GRAY_DARK_COLOR, 'Marker', '*');
end
p=3; plot(polyshape([p-w p-w p+w p+w], log10([mt_nad_nadh_ratio_low mt_nad_nadh_ratio_high mt_nad_nadh_ratio_high mt_nad_nadh_ratio_low])),'FaceColor', BLUE_COLOR, 'FaceAlpha',0.7);
%mean_val=(log10(code_mfa_mt_nad_nadh_ratio_low)+log10(code_mfa_mt_nad_nadh_ratio_high))/2; errorbar(p,mean_val,log10(code_mfa_mt_nad_nadh_ratio_high)-mean_val,log10(code_mfa_mt_nad_nadh_ratio_high)-mean_val,'color',BLUE_COLOR,'linewidth',3, 'LineStyle', 'none');
mean_val=(log10(code_mfa_mt_nad_nadh_ratio_low)+log10(code_mfa_mt_nad_nadh_ratio_high))/2; errorbar(p,mean_val,log10(code_mfa_mt_nad_nadh_ratio_high)-mean_val,log10(code_mfa_mt_nad_nadh_ratio_high)-mean_val,'color',RED_COLOR,'linewidth',3, 'LineStyle', 'none');
if(PLOT_ALL_PAPERS)
    plot(p,log10(mt_NAD_NADH_ratio_won),'bx', 'MarkerSize',10, 'LineWidth', 1.5);
    plot(p,log10(mt_NAD_NADH_ratio_sabatini),'b+', 'MarkerSize',10, 'LineWidth', 1.5);
    plot(p,log10(mt_NAD_NADH_ratio_krebs),'bo', 'MarkerSize',10, 'LineWidth', 1.5);
    plot(p,log10(mt_NAD_NADH_ratio_elmar),'b^', 'MarkerSize',10, 'LineWidth', 1.5);
    legend({'WC measured','CY TD','CY Code-MFA', 'CY Sabatini (Hela)', 'CY Krebs (Rat liver)', 'CY Sun (Hela)', 'CY Hedeskov (Mouse pancreas)', 'CY zhao (BEAS-2B)', 'MT TD (WC&CY)', 'MT Code-MFA', 'MT Won (Hela)', 'MT Sabatini (Hela)', 'MT Krebs (Rat liver)', 'MT Elmar (Rat liver)'}, 'FontSize', 11, 'Location', 'southwest');
else
    all_literature_vals = [mt_NAD_NADH_ratio_won mt_NAD_NADH_ratio_sabatini mt_NAD_NADH_ratio_krebs mt_NAD_NADH_ratio_elmar];
    plot(p*ones(1,length(all_literature_vals)),log10(all_literature_vals), 'MarkerSize',12,  'LineWidth', 1.6, 'LineStyle', 'none', 'Color', GRAY_DARK_COLOR, 'Marker', '*');
    %legend({'WC measured','CY TD','CY Code-MFA', 'CY Literature', 'MT TD (WC&CY)', 'MT Code-MFA', 'MT Literature'}, 'FontSize', 20, 'Location', 'southwest');
    f=get(gca,'Children');
    legend([f(7), f(6), f(5), f(4)],{' WC measured',' Thermodynamics',' CoDe-MFA', ' Literature'}, 'FontSize', 30, 'Location', 'southwest');
end
 
xlabel('');
ylabel('NAD+/NADH Ratio [log10]');
xticks([1 2 3]);
xticklabels({'WC', 'CY', 'MT'});
set(gcf,'color','w');
set(gca, 'FontSize', 30);
ylim([-5 5]);
yticks([-5:1:5]);
box on

 

% figure for NAD and NADH concentrations
subplot(1,2,2);
grid on;
hold on;

xlim([0.5 8.5]);

w=0.25;
%p=1; plot(polyshape([p-w p-w p+w p+w], log10([NAD_con_min NAD_con_max NAD_con_max NAD_con_min])),'FaceColor', RED_COLOR, 'FaceAlpha',0.7);
p=1; plot(polyshape([p-w p-w p+w p+w], log10([NAD_con_min NAD_con_max NAD_con_max NAD_con_min])),'FaceColor', GRAY_COLOR, 'FaceAlpha',0.7);
p=2; plot(polyshape([p-w p-w p+w p+w], log10([cy_nad_low cy_nad_high cy_nad_high cy_nad_low])),'FaceColor', BLUE_COLOR, 'FaceAlpha',0.7);
%mean_val=(log10(code_mfa_cy_nad_low)+log10(code_mfa_cy_nad_high))/2; errorbar(p,mean_val,log10(code_mfa_cy_nad_high)-mean_val,log10(code_mfa_cy_nad_high)-mean_val,'color',BLUE_COLOR,'linewidth',3, 'LineStyle', 'none');
mean_val=(log10(code_mfa_cy_nad_low)+log10(code_mfa_cy_nad_high))/2; errorbar(p,mean_val,log10(code_mfa_cy_nad_high)-mean_val,log10(code_mfa_cy_nad_high)-mean_val,'color',RED_COLOR,'linewidth',3, 'LineStyle', 'none');
p=3; plot(polyshape([p-w p-w p+w p+w], log10([mt_nad_low mt_nad_high mt_nad_high mt_nad_low])),'FaceColor', BLUE_COLOR, 'FaceAlpha',0.7);
%mean_val=(log10(code_mfa_mt_nad_low)+log10(code_mfa_mt_nad_high))/2; errorbar(p,mean_val,log10(code_mfa_mt_nad_high)-mean_val,log10(code_mfa_mt_nad_high)-mean_val,'color',BLUE_COLOR,'linewidth',3, 'LineStyle', 'none');
mean_val=(log10(code_mfa_mt_nad_low)+log10(code_mfa_mt_nad_high))/2; errorbar(p,mean_val,log10(code_mfa_mt_nad_high)-mean_val,log10(code_mfa_mt_nad_high)-mean_val,'color',RED_COLOR,'linewidth',3, 'LineStyle', 'none');

w=0.25;
%p=5; plot(polyshape([p-w p-w p+w p+w], log10([NADH_con_min NADH_con_max NADH_con_max NADH_con_min])),'FaceColor', RED_COLOR, 'FaceAlpha',0.7);
p=6; plot(polyshape([p-w p-w p+w p+w], log10([NADH_con_min NADH_con_max NADH_con_max NADH_con_min])),'FaceColor', GRAY_COLOR, 'FaceAlpha',0.7);
p=7; plot(polyshape([p-w p-w p+w p+w], log10([cy_nadh_low cy_nadh_high cy_nadh_high cy_nadh_low])),'FaceColor', BLUE_COLOR, 'FaceAlpha',0.7);
%mean_val=(log10(code_mfa_cy_nadh_low)+log10(code_mfa_cy_nadh_high))/2; errorbar(p,mean_val,log10(code_mfa_cy_nadh_high)-mean_val,log10(code_mfa_cy_nadh_high)-mean_val,'color',BLUE_COLOR,'linewidth',3, 'LineStyle', 'none');
mean_val=(log10(code_mfa_cy_nadh_low)+log10(code_mfa_cy_nadh_high))/2; errorbar(p,mean_val,log10(code_mfa_cy_nadh_high)-mean_val,log10(code_mfa_cy_nadh_high)-mean_val,'color',RED_COLOR,'linewidth',3, 'LineStyle', 'none');
p=8; plot(polyshape([p-w p-w p+w p+w], log10([mt_nadh_low mt_nadh_high mt_nadh_high mt_nadh_low])),'FaceColor', BLUE_COLOR, 'FaceAlpha',0.7);
%mean_val=(log10(code_mfa_mt_nadh_low)+log10(code_mfa_mt_nadh_high))/2; errorbar(p,mean_val,log10(code_mfa_mt_nadh_high)-mean_val,log10(code_mfa_mt_nadh_high)-mean_val,'color',BLUE_COLOR,'linewidth',3, 'LineStyle', 'none');
mean_val=(log10(code_mfa_mt_nadh_low)+log10(code_mfa_mt_nadh_high))/2; errorbar(p,mean_val,log10(code_mfa_mt_nadh_high)-mean_val,log10(code_mfa_mt_nadh_high)-mean_val,'color',RED_COLOR,'linewidth',3, 'LineStyle', 'none');


xlabel('');
ylabel('Concentrations [log10(mM)]');
xticks([1 2 3 6 7 8]);
xticklabels({'WC', 'CY', 'MT', 'WC', 'CY', 'MT'});
title('NAD+                         NADH', 'Units', 'normalized', 'Position', [0.5, -0.14, 0]); % Set Title with correct Position
f=get(gca,'Children');
legend([f(10), f(9), f(8)],{' WC measured',' Thermodynamics',' CoDe-MFA'}, 'FontSize', 30, 'Location', 'northeast');
set(gcf,'color','w');
set(gca, 'FontSize', 30);
ylim([-5 1.5]);
yticks([-5:1:1]);
box on

set(gcf,'PaperSize',[50 30]);
s = sprintf('./output_images/1gh.pdf');
saveas(gcf, s);
%suptitle('Without constraining Malate Aspartate shuttle')


% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% do the same, but now constraint the malate aspartate shuttle to work in
% the canonical view
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
figure('units','normalized','outerposition',[0 0 1 1])

% figure for NAD/NADH ratio
subplot(1,2,1);
grid on;
hold on;


% Mitochondrial lower bound, based on Malate Aspartate shuttle - assuming it is working in the
% canonical view
G0_Malate_Aspartate_shuttle_cycle   = [-20.2 2];
n_G0_Malate_Aspartate_shuttle_cycle = G0_Malate_Aspartate_shuttle_cycle(1) + randn(1,1000) * G0_Malate_Aspartate_shuttle_cycle(2);
mito_NAD_NADH_ratio = (1./cyto_NADH_NAD_ratio)./(exp(-n_G0_Malate_Aspartate_shuttle_cycle/RT));
mt_nad_nadh_ratio_based_on_malate_aspartate_shuttle_low = prctile(mito_NAD_NADH_ratio,5);


% assuming mallate aspartate shuttle is working in the canonical way
A3=[0 0 -1 ((cy_nad_nadh_ratio_low)/(exp((-G0_Malate_Aspartate_shuttle_cycle(1)+2*G0_Malate_Aspartate_shuttle_cycle(2))/RT)))];
b3=0;

A=[A1;A2;A3];
b=[b1;b2;b3];

f = [1 0 0 0];  % min cyto NAD
[initial_values_new,fval,exitflag] = linprog(f, A, b, [], [], lb, ub);
cy_nad_low = fval;
f = [-1 0 0 0];  % max cyto NAD
[initial_values_new,fval,exitflag] = linprog(f, A, b, [], [], lb, ub);
cy_nad_high = -fval;
f = [0 1 0 0];  % min cyto NADH
[initial_values_new,fval,exitflag] = linprog(f, A, b, [], [], lb, ub);
cy_nadh_low = fval;
f = [0 -1 0 0];  % max cyto NADH
[initial_values_new,fval,exitflag] = linprog(f, A, b, [], [], lb, ub);
cy_nadh_high = -fval;
f = [0 0 1 0];  % min mito NAD
[initial_values_new,fval,exitflag] = linprog(f, A, b, [], [], lb, ub);
mt_nad_low = fval;
f = [0 0 -1 0];  % max mito NAD
[initial_values_new,fval,exitflag] = linprog(f, A, b, [], [], lb, ub);
mt_nad_high = -fval;
f = [0 0 0 1];  % min mito NADH
[initial_values_new,fval,exitflag] = linprog(f, A, b, [], [], lb, ub);
mt_nadh_low = fval;                        
f = [0 0 0 -1];  % max mito NADH
[initial_values_new,fval,exitflag] = linprog(f, A, b, [], [], lb, ub);
mt_nadh_high = -fval;
 

% check upper bound on mitochondrial NAD/NADH ratio
Beq=0;
q_factor = 1.1; 
max_ratio = 1e6;
f = [0 0 0 0];  % min mito NADH
mt_nad_nadh_ratio_constraint_ma_shuttle_low = mt_nad_nadh_ratio_based_on_malate_aspartate_shuttle_low;
mt_nad_nadh_ratio_constraint_ma_shuttle_high=mt_nad_nadh_ratio_constraint_ma_shuttle_low;
while(mt_nad_nadh_ratio_constraint_ma_shuttle_high<max_ratio)        
    Aeq=[0 0 1 -mt_nad_nadh_ratio_constraint_ma_shuttle_high];
    [initial_values_new,fval,exitflag] = linprog(f, A, b, Aeq, Beq, lb, ub);
    % no solution was found
    if(exitflag~=1)
        break;
    end
    mt_nad_nadh_ratio_constraint_ma_shuttle_high=mt_nad_nadh_ratio_constraint_ma_shuttle_high*q_factor;
end


xlim([0.5 3.5]); 
  
w=0.25;
p=1; plot(polyshape([p-w p-w p+w p+w], log10([wc_nad_nadh_ratio_low wc_nad_nadh_ratio_high wc_nad_nadh_ratio_high wc_nad_nadh_ratio_low])),'FaceColor', RED_COLOR, 'FaceAlpha',0.7);
p=2; plot(polyshape([p-w p-w p+w p+w], log10([cy_nad_nadh_ratio_low cy_nad_nadh_ratio_high cy_nad_nadh_ratio_high cy_nad_nadh_ratio_low])),'FaceColor', GREEN_COLOR, 'FaceAlpha',0.7);
mean_val=(log10(code_mfa_cy_nad_nadh_ratio_low)+log10(code_mfa_cy_nad_nadh_ratio_high))/2; errorbar(p,mean_val,log10(code_mfa_cy_nad_nadh_ratio_high)-mean_val,log10(code_mfa_cy_nad_nadh_ratio_high)-mean_val,'color',GREEN_COLOR,'linewidth',3, 'LineStyle', 'none');
plot(p,log10(cy_NAD_NADH_ratio_sabatini),'g+', 'MarkerSize',10, 'LineWidth', 1.5);
plot(p,log10(cy_NAD_NADH_ratio_krebs),'go', 'MarkerSize',10, 'LineWidth', 1.5);
plot(p,log10(cy_NAD_NADH_ratio_sun),'gs', 'MarkerSize',10, 'LineWidth', 1.5);
plot(p,log10(cy_NAD_NADH_ratio_hedeskov),'gd', 'MarkerSize',10, 'LineWidth', 1.5);
plot(p,log10(cy_NAD_NADH_ratio_zhao),'g*', 'MarkerSize',10, 'LineWidth', 1.5);
p=3; plot(polyshape([p-w p-w p+w p+w], log10([mt_nad_nadh_ratio_constraint_ma_shuttle_low mt_nad_nadh_ratio_constraint_ma_shuttle_high mt_nad_nadh_ratio_constraint_ma_shuttle_high mt_nad_nadh_ratio_constraint_ma_shuttle_low])),'FaceColor', BLUE_COLOR, 'FaceAlpha',0.7);
mean_val=(log10(code_mfa_mt_nad_nadh_ratio_low)+log10(code_mfa_mt_nad_nadh_ratio_high))/2; errorbar(p,mean_val,log10(code_mfa_mt_nad_nadh_ratio_high)-mean_val,log10(code_mfa_mt_nad_nadh_ratio_high)-mean_val,'color',BLUE_COLOR,'linewidth',3, 'LineStyle', 'none');
plot(p,log10(mt_NAD_NADH_ratio_won),'bx', 'MarkerSize',10, 'LineWidth', 1.5);
plot(p,log10(mt_NAD_NADH_ratio_sabatini),'b+', 'MarkerSize',10, 'LineWidth', 1.5);
plot(p,log10(mt_NAD_NADH_ratio_krebs),'bo', 'MarkerSize',10, 'LineWidth', 1.5);
plot(p,log10(mt_NAD_NADH_ratio_elmar),'b^', 'MarkerSize',10, 'LineWidth', 1.5);
legend({'WC measured','CY TD','CY Code-MFA', 'CY Sabatini (Hela)', 'CY Krebs (Rat liver)', 'CY Sun (Hela)', 'CY Hedeskov (Mouse pancreas)', 'CY zhao (BEAS-2B)', 'MT TD (WC&CY)', 'MT Code-MFA', 'MT Won (Hela)', 'MT Sabatini (Hela)', 'MT Krebs (Rat liver)', 'MT Elmar (Rat liver)'}, 'FontSize', 11, 'Location', 'southwest');



xlabel('');
ylabel('NAD+/NADH Ratio [log10]');
xticks([1 2 3]);
xticklabels({'WC', 'CY', 'MT'});
set(gcf,'color','w');
set(gca, 'FontSize', 16);
box on



% figure for NAD and NADH concentrations
subplot(1,2,2);
grid on;
hold on; 

xlim([0.5 8.5]);

w=0.25;
p=1; plot(polyshape([p-w p-w p+w p+w], log10([NAD_con_min NAD_con_max NAD_con_max NAD_con_min])),'FaceColor', RED_COLOR, 'FaceAlpha',0.7);
p=2; plot(polyshape([p-w p-w p+w p+w], log10([cy_nad_low cy_nad_high cy_nad_high cy_nad_low])),'FaceColor', GREEN_COLOR, 'FaceAlpha',0.7);
mean_val=(log10(code_mfa_cy_nad_low)+log10(code_mfa_cy_nad_high))/2; errorbar(p,mean_val,log10(code_mfa_cy_nad_high)-mean_val,log10(code_mfa_cy_nad_high)-mean_val,'color',GREEN_COLOR,'linewidth',3, 'LineStyle', 'none');
p=3; plot(polyshape([p-w p-w p+w p+w], log10([mt_nad_low mt_nad_high mt_nad_high mt_nad_low])),'FaceColor', BLUE_COLOR, 'FaceAlpha',0.7);
mean_val=(log10(code_mfa_mt_nad_low)+log10(code_mfa_mt_nad_high))/2; errorbar(p,mean_val,log10(code_mfa_mt_nad_high)-mean_val,log10(code_mfa_mt_nad_high)-mean_val,'color',BLUE_COLOR,'linewidth',3, 'LineStyle', 'none');

w=0.25;
p=6; plot(polyshape([p-w p-w p+w p+w], log10([NADH_con_min NADH_con_max NADH_con_max NADH_con_min])),'FaceColor', RED_COLOR, 'FaceAlpha',0.7);
p=7; plot(polyshape([p-w p-w p+w p+w], log10([cy_nadh_low cy_nadh_high cy_nadh_high cy_nadh_low])),'FaceColor', GREEN_COLOR, 'FaceAlpha',0.7);
mean_val=(log10(code_mfa_cy_nadh_low)+log10(code_mfa_cy_nadh_high))/2; errorbar(p,mean_val,log10(code_mfa_cy_nadh_high)-mean_val,log10(code_mfa_cy_nadh_high)-mean_val,'color',GREEN_COLOR,'linewidth',3, 'LineStyle', 'none');
p=8; plot(polyshape([p-w p-w p+w p+w], log10([mt_nadh_low mt_nadh_high mt_nadh_high mt_nadh_low])),'FaceColor', BLUE_COLOR, 'FaceAlpha',0.7);
mean_val=(log10(code_mfa_mt_nadh_low)+log10(code_mfa_mt_nadh_high))/2; errorbar(p,mean_val,log10(code_mfa_mt_nadh_high)-mean_val,log10(code_mfa_mt_nadh_high)-mean_val,'color',BLUE_COLOR,'linewidth',3, 'LineStyle', 'none');



xlabel('');
ylabel('Concentrations [log10(mM)]');
xticks([1 2 3 6 7 8]);
xticklabels({'WC', 'CY', 'MT', 'WC', 'CY', 'MT'});
title('NAD+                                                             NADH', 'Units', 'normalized', 'Position', [0.5, -0.08, 0]); % Set Title with correct Position
set(gcf,'color','w');
set(gca, 'FontSize', 16);
box on
suptitle('Constraining Malate Aspartate shuttle')







%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% NADP/NADPH
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
CO2_con = 1.2;

% Ratio from literature - Mitochondria
% SIES, H., AKERBOOM, T. P. M., & TAGER, J. M. (1977). Mitochondrial and Cytosolic NADPH Systems and Isocitrate Dehydrogenase Indicator Metabolites during Ureogenesis from Ammonia in Isolated Rat Hepatocytes. European Journal of Biochemistry. https://doi.org/10.1111/j.1432-1033.1977.tb11253.x
mt_NADP_NADPH_ratio_sies = 1/83.5;
% Sallin, O., Reymond, L., Gondrand, C., Raith, F., Koch, B., & Johnsson, K. (2018). Semisynthetic biosensors for mapping cellular concentrations of nicotinamide adenine dinucleotides. ELife. https://doi.org/10.7554/eLife.32638
mt_NADP_NADPH_ratio_sallin = 1/200;

% Ratio from literature - Cytosol
% Hedeskov, C. J., Capito, K., & Thams, P. (1987). Cytosolic ratios of free [NADPH]/[NADP+] and [NADH]/[NAD+] in mouse pancreatic islets, and nutrient-induced insulin secretion. The Biochemical Journal, 241(1), 161–167. https://doi.org/10.1042/bj2410161
cy_NADP_NADPH_ratio_hedeskov = 1/30;
% Veech, R. L., Eggleston, L. V, & Krebs, H. a. (1969). The redox state of free nicotinamide-adenine dinucleotide phosphate in the cytoplasm of rat liver. The Biochemical Journal, 115(4), 609–619. https://doi.org/10.1042/bj1150609a
cy_NADP_NADPH_ratio_krebs = 1/100;
% Siess, E. A., Brocks, D. G., Lattke, H. K., & Wieland, O. H. (1977). Effect of glucagon on metabolite compartmentation in isolated rat liver cells during gluconeogenesis from lactate. Biochemical Journal. https://doi.org/10.1042/bj1660225
cy_NADP_NADPH_ratio_elmar = 0.038;
% Sallin, O., Reymond, L., Gondrand, C., Raith, F., Koch, B., & Johnsson, K. (2018). Semisynthetic biosensors for mapping cellular concentrations of nicotinamide adenine dinucleotides. ELife. https://doi.org/10.7554/eLife.32638
cy_NADP_NADPH_ratio_sallin = 1/75;
% Park, J. O., Rubin, S. A., Xu, Y. F., Amador-Noguez, D., Fan, J., Shlomi, T., & Rabinowitz, J. D. (2016). Metabolite concentrations, fluxes and free energies imply efficient enzyme usage. Nature Chemical Biology, 12(7), 482–489. https://doi.org/10.1038/nchembio.2077
dG_park                         = -11.15;     % iBMK cells
PG6_con_park                    = 1.65*10^-2; % iBMK cells
ribulose_5_phosphate_con_park   = 5.27*10^-3; % iBMK cells
G0_PG6_to_ribulose_5_phosphate  = -5.9;
cy_NADPH_NADP_ratio_park = (exp((dG_park-G0_PG6_to_ribulose_5_phosphate)/RT))/((ribulose_5_phosphate_con_park*CO2_con)/(PG6_con_park));
cy_NADP_NADPH_ratio_park = 1/cy_NADPH_NADP_ratio_park;


% Metabolite concentrations [mM]
index_of_metabolite_in_wc_measured_vector = find(ismember(all_met_names_with_WC_measured_con,'NADP'));
NADP_con = all_met_WC_measured_con(index_of_metabolite_in_wc_measured_vector,:);
NADP_con_min = NADP_con(1)-NADP_con(2)*2;
NADP_con_max = NADP_con(1)+NADP_con(2)*2;

index_of_metabolite_in_wc_measured_vector = find(ismember(all_met_names_with_WC_measured_con,'NADPH'));
NADPH_con = all_met_WC_measured_con(index_of_metabolite_in_wc_measured_vector,:);
NADPH_con_min = NADPH_con(1)-NADPH_con(2)*2;
NADPH_con_max = NADPH_con(1)+NADPH_con(2)*2;

mito_NADP_con_max    = NADP_con_max/mito_volume;
mito_NADPH_con_max   = NADPH_con_max/mito_volume;
cyto_NADP_con_max    = NADP_con_max/cyto_volume;
cyto_NADPH_con_max   = NADPH_con_max/cyto_volume;


% 6phosphogluconate_CY + NADP_CY => Ribose5P_CY + NADPH_CY + CO2
index_of_metabolite_in_wc_measured_vector = find(ismember(all_met_names_with_WC_measured_con,'6phosphogluconate'));
PG6_con = all_met_WC_measured_con(index_of_metabolite_in_wc_measured_vector,:);
index_of_metabolite_in_wc_measured_vector = find(ismember(all_met_names_with_WC_measured_con,'Ribose5P'));
R5P_con = all_met_WC_measured_con(index_of_metabolite_in_wc_measured_vector,:);
G0_PG6_to_R5P = [-7.8 3.2];

% find m+5 labeling for PG6
for(i=1:length(WC_known_metabolites_idv{2}))
    if(strcmp(WC_known_metabolites_idv{2}{i}.met_name,'6phosphogluconate'))
        PG6_m5 = [WC_known_metabolites_idv{2}{i}.idv(6) max(0.01,WC_known_metabolites_idv{2}{i}.idv_variance(6))];
    end
    if(strcmp(WC_known_metabolites_idv{2}{i}.met_name,'Ribose5P'))
        R5P_m5 = [WC_known_metabolites_idv{2}{i}.idv(6) max(0.01,WC_known_metabolites_idv{2}{i}.idv_variance(6))];
    end    
end

PG6_to_R5P_forward_div_backward_flux = [PG6_m5(1)/R5P_m5(1) max(PG6_m5(2),R5P_m5(2))]; % b/(f-b)=0.04/0.96

n_PG6_con                 = PG6_con(1) + randn(1,1000) * PG6_con(2);
n_R5P_con                 = R5P_con(1) + randn(1,1000) * R5P_con(2);
n_G0_PG6_to_R5P           = G0_PG6_to_R5P(1) + randn(1,1000) * G0_PG6_to_R5P(2);
n_PG6_to_R5P_forward_div_backward_flux = PG6_to_R5P_forward_div_backward_flux(1) + randn(1,1000) * PG6_to_R5P_forward_div_backward_flux(2);

cyto_NADPH_NADP_ratio = exp(log(n_PG6_to_R5P_forward_div_backward_flux) - n_G0_PG6_to_R5P/RT - log(n_R5P_con) - log(CO2_con) + log(n_PG6_con));
cyto_NADP_NADPH_ratio = 1./cyto_NADPH_NADP_ratio;

figure('units','normalized','outerposition',[0 0 1 1])

% figure for NADP/NADPH ratio
subplot(1,2,1);
grid on;
hold on;
 
% WC measured
wc_nadp_nadph_ratio_low   = NADP_con_min/NADPH_con_max;
wc_nadp_nadph_ratio_high  = NADP_con_max/NADPH_con_min;

% Cytosolic, based on 6Pg => R5P, while gibbs free energy is calculated
% based on flux force relationship
cy_nadph_nadp_ratio_low  = prctile(cyto_NADPH_NADP_ratio,5);
cy_nadph_nadp_ratio_high = prctile(cyto_NADPH_NADP_ratio,95);
cy_nadp_nadph_ratio_low   = 1/cy_nadph_nadp_ratio_high;
cy_nadp_nadph_ratio_high  = 1/cy_nadph_nadp_ratio_low;

  
 
% find min and max for cytosolic NADP and NADPH
% linear programming [cyto_NADP cyto_NADPH mito_NADP mito_NADPH]
lb = [MINIMUM_CONCENTRATION MINIMUM_CONCENTRATION MINIMUM_CONCENTRATION MINIMUM_CONCENTRATION];
ub = [cyto_NADP_con_max cyto_NADPH_con_max mito_NADP_con_max mito_NADPH_con_max];
% use cytosolic NADP/NADPH ratio
A1=[-1 cy_nadp_nadph_ratio_low 0 0; 1 -cy_nadp_nadph_ratio_high 0 0];
b1=[0;0];
% cyto_volume*cyto_NADP+mito_volume*mito_NADP=WC_NADP
A2=[cyto_volume 0 mito_volume 0; -cyto_volume 0 -mito_volume 0; 0 cyto_volume 0 mito_volume; 0 -cyto_volume 0 -mito_volume];
b2=[NADP_con_max; -NADP_con_min; NADPH_con_max; -NADPH_con_min];

A=[A1;A2];
b=[b1;b2];

f = [1 0 0 0];  % min cyto NADP
[initial_values_new,fval,exitflag] = linprog(f, A, b, [], [], lb, ub);
cy_nadp_low = fval;
f = [-1 0 0 0];  % max cyto NADP
[initial_values_new,fval,exitflag] = linprog(f, A, b, [], [], lb, ub);
cy_nadp_high = -fval;
f = [0 1 0 0];  % min cyto NADPH
[initial_values_new,fval,exitflag] = linprog(f, A, b, [], [], lb, ub);
cy_nadph_low = fval;
f = [0 -1 0 0];  % max cyto NADPH
[initial_values_new,fval,exitflag] = linprog(f, A, b, [], [], lb, ub);
cy_nadph_high = -fval;
f = [0 0 1 0];  % min mito NADP
[initial_values_new,fval,exitflag] = linprog(f, A, b, [], [], lb, ub);
mt_nadp_low = fval;
f = [0 0 -1 0];  % max mito NADP
[initial_values_new,fval,exitflag] = linprog(f, A, b, [], [], lb, ub);
mt_nadp_high = -fval;
f = [0 0 0 1];  % min mito NADPH
[initial_values_new,fval,exitflag] = linprog(f, A, b, [], [], lb, ub);
mt_nadph_low = fval;                        
f = [0 0 0 -1];  % max mito NADPH
[initial_values_new,fval,exitflag] = linprog(f, A, b, [], [], lb, ub);
mt_nadph_high = -fval;
 

% check upper bound on mitochondrial NADP/NADPH ratio
Beq=0;
q_factor = 1.1; 
min_ratio = 1e-6;
max_ratio = 1e6;
mt_nadp_nadph_ratio_low = min_ratio;
f = [0 0 0 0];  % min mito NADPH
while(mt_nadp_nadph_ratio_low<max_ratio)
    Aeq=[0 0 1 -mt_nadp_nadph_ratio_low];
    [initial_values_new,fval,exitflag] = linprog(f, A, b, Aeq, Beq, lb, ub);
    % no solution was found
    if(exitflag==1)
        break;
    end 
    mt_nadp_nadph_ratio_low=mt_nadp_nadph_ratio_low*q_factor;
end 
mt_nadp_nadph_ratio_high=mt_nadp_nadph_ratio_low;
while(mt_nadp_nadph_ratio_high<max_ratio)        
    Aeq=[0 0 1 -mt_nadp_nadph_ratio_high];
    [initial_values_new,fval,exitflag] = linprog(f, A, b, Aeq, Beq, lb, ub);
    % no solution was found
    if(exitflag~=1)
        break;
    end
    mt_nadp_nadph_ratio_high=mt_nadp_nadph_ratio_high*q_factor;
end


xlim([0.5 3.5]); 
  
w=0.25;
%p=1; plot(polyshape([p-w p-w p+w p+w], log10([wc_nadp_nadph_ratio_low wc_nadp_nadph_ratio_high wc_nadp_nadph_ratio_high wc_nadp_nadph_ratio_low])),'FaceColor', RED_COLOR, 'FaceAlpha',0.7);
p=1; plot(polyshape([p-w p-w p+w p+w], log10([wc_nadp_nadph_ratio_low wc_nadp_nadph_ratio_high wc_nadp_nadph_ratio_high wc_nadp_nadph_ratio_low])),'FaceColor', GRAY_COLOR, 'FaceAlpha',0.7);
p=2; plot(polyshape([p-w p-w p+w p+w], log10([cy_nadp_nadph_ratio_low cy_nadp_nadph_ratio_high cy_nadp_nadph_ratio_high cy_nadp_nadph_ratio_low])),'FaceColor', BLUE_COLOR, 'FaceAlpha',0.7);
% mean_val=(log10(code_mfa_cy_nadp_nadph_ratio_low)+log10(code_mfa_cy_nadp_nadph_ratio_high))/2; errorbar(p,mean_val,log10(code_mfa_cy_nadp_nadph_ratio_high)-mean_val,log10(code_mfa_cy_nadp_nadph_ratio_high)-mean_val,'color',BLUE_COLOR,'linewidth',3, 'LineStyle', 'none');
mean_val=(log10(code_mfa_cy_nadp_nadph_ratio_low)+log10(code_mfa_cy_nadp_nadph_ratio_high))/2; errorbar(p,mean_val,log10(code_mfa_cy_nadp_nadph_ratio_high)-mean_val,log10(code_mfa_cy_nadp_nadph_ratio_high)-mean_val,'color',RED_COLOR,'linewidth',3, 'LineStyle', 'none');
if(PLOT_ALL_PAPERS)
    plot(p,log10(cy_NADP_NADPH_ratio_hedeskov),'gd', 'MarkerSize',10, 'LineWidth', 1.5);
    plot(p,log10(cy_NADP_NADPH_ratio_krebs),'go', 'MarkerSize',10, 'LineWidth', 1.5);
    plot(p,log10(cy_NADP_NADPH_ratio_elmar),'g+', 'MarkerSize',10, 'LineWidth', 1.5);
    plot(p,log10(cy_NADP_NADPH_ratio_sallin),'gs', 'MarkerSize',10, 'LineWidth', 1.5);
    plot(p,log10(cy_NADP_NADPH_ratio_park),'g*', 'MarkerSize',10, 'LineWidth', 1.5);
else
    all_literature_vals = [cy_NADP_NADPH_ratio_hedeskov cy_NADP_NADPH_ratio_krebs cy_NADP_NADPH_ratio_elmar cy_NADP_NADPH_ratio_sallin cy_NADP_NADPH_ratio_park];
%     plot(p*ones(1,length(all_literature_vals)),log10([cy_NADP_NADPH_ratio_hedeskov cy_NADP_NADPH_ratio_krebs cy_NADP_NADPH_ratio_elmar cy_NADP_NADPH_ratio_sallin cy_NADP_NADPH_ratio_park]), 'MarkerSize',12, 'LineWidth', 1.5, 'LineStyle', 'none', 'Color', BLUE_COLOR, 'Marker', '*');
    plot(p*ones(1,length(all_literature_vals)),log10([cy_NADP_NADPH_ratio_hedeskov cy_NADP_NADPH_ratio_krebs cy_NADP_NADPH_ratio_elmar cy_NADP_NADPH_ratio_sallin cy_NADP_NADPH_ratio_park]), 'MarkerSize',12, 'LineWidth', 1.6, 'LineStyle', 'none', 'Color', GRAY_DARK_COLOR, 'Marker', '*');
end
p=3; plot(polyshape([p-w p-w p+w p+w], log10([mt_nadp_nadph_ratio_low mt_nadp_nadph_ratio_high mt_nadp_nadph_ratio_high mt_nadp_nadph_ratio_low])),'FaceColor', BLUE_COLOR, 'FaceAlpha',0.7);
%mean_val=(log10(code_mfa_mt_nadp_nadph_ratio_low)+log10(code_mfa_mt_nadp_nadph_ratio_high))/2; errorbar(p,mean_val,log10(code_mfa_mt_nadp_nadph_ratio_high)-mean_val,log10(code_mfa_mt_nadp_nadph_ratio_high)-mean_val,'color',BLUE_COLOR,'linewidth',3, 'LineStyle', 'none');
mean_val=(log10(code_mfa_mt_nadp_nadph_ratio_low)+log10(code_mfa_mt_nadp_nadph_ratio_high))/2; errorbar(p,mean_val,log10(code_mfa_mt_nadp_nadph_ratio_high)-mean_val,log10(code_mfa_mt_nadp_nadph_ratio_high)-mean_val,'color',RED_COLOR,'linewidth',3, 'LineStyle', 'none');
if(PLOT_ALL_PAPERS)
    plot(p,log10(mt_NADP_NADPH_ratio_sies),'b+', 'MarkerSize',10, 'LineWidth', 1.5);
    plot(p,log10(mt_NADP_NADPH_ratio_sallin),'bs', 'MarkerSize',10, 'LineWidth', 1.5);
    legend({'WC measured','CY TD (6PGD ratio)','CY Code-MFA', 'CY Hedeskov (Mouse pancreas)', 'CY Krebs (Rat liver)', 'CY Elmar (Rat liver)', 'CY Sallin (Hela)', 'CY Park (iBMK)','MT TD (WC&CY)', 'MT Code-MFA', 'MT Sies (Rat Hepatocytes)', 'MT Sallin (U20S)'}, 'FontSize', 11, 'Location', 'southwest');
else
    all_literature_vals = [mt_NADP_NADPH_ratio_sies mt_NADP_NADPH_ratio_sallin];
    %plot(p*ones(1,length(all_literature_vals)),log10(all_literature_vals), 'MarkerSize',12, 'LineWidth', 1.5, 'LineStyle', 'none', 'Color', BLUE_COLOR, 'Marker', '*');
    plot(p*ones(1,length(all_literature_vals)),log10(all_literature_vals), 'MarkerSize',12, 'LineWidth', 1.6, 'LineStyle', 'none', 'Color', GRAY_DARK_COLOR, 'Marker', '*');
%     legend({'WC measured','CY TD (6PGD ratio)','CY Code-MFA', 'CY Literature','MT TD (WC&CY)', 'MT Code-MFA', 'MT Literature'}, 'FontSize', 20, 'Location', 'southwest');
    f=get(gca,'Children');
    legend([f(7), f(6), f(5), f(4)],{' WC measured',' Thermodynamics',' CoDe-MFA', ' Literature'}, 'FontSize', 30, 'Location', 'southwest');
end




xlabel('');
ylabel('NADP+/NADPH Ratio [log10]');
xticks([1 2 3]);
xticklabels({'WC', 'CY', 'MT'});
set(gcf,'color','w');
set(gca, 'FontSize', 30);
ylim([-5 5]);
yticks([-5:1:5]);
box on



% figure for NADP and NADPH concentrations
subplot(1,2,2);
grid on;
hold on;

xlim([0.5 8.5]);

w=0.25;
%p=1; plot(polyshape([p-w p-w p+w p+w], log10([NADP_con_min NADP_con_max NADP_con_max NADP_con_min])),'FaceColor', RED_COLOR, 'FaceAlpha',0.7);
p=1; plot(polyshape([p-w p-w p+w p+w], log10([NADP_con_min NADP_con_max NADP_con_max NADP_con_min])),'FaceColor', GRAY_COLOR, 'FaceAlpha',0.7);
p=2; plot(polyshape([p-w p-w p+w p+w], log10([cy_nadp_low cy_nadp_high cy_nadp_high cy_nadp_low])),'FaceColor', BLUE_COLOR, 'FaceAlpha',0.7);
%mean_val=(log10(code_mfa_cy_nadp_low)+log10(code_mfa_cy_nadp_high))/2; errorbar(p,mean_val,log10(code_mfa_cy_nadp_high)-mean_val,log10(code_mfa_cy_nadp_high)-mean_val,'color',BLUE_COLOR,'linewidth',3, 'LineStyle', 'none');
mean_val=(log10(code_mfa_cy_nadp_low)+log10(code_mfa_cy_nadp_high))/2; errorbar(p,mean_val,log10(code_mfa_cy_nadp_high)-mean_val,log10(code_mfa_cy_nadp_high)-mean_val,'color',RED_COLOR,'linewidth',3, 'LineStyle', 'none');
p=3; plot(polyshape([p-w p-w p+w p+w], log10([mt_nadp_low mt_nadp_high mt_nadp_high mt_nadp_low])),'FaceColor', BLUE_COLOR, 'FaceAlpha',0.7);
%mean_val=(log10(code_mfa_mt_nadp_low)+log10(code_mfa_mt_nadp_high))/2; errorbar(p,mean_val,log10(code_mfa_mt_nadp_high)-mean_val,log10(code_mfa_mt_nadp_high)-mean_val,'color',BLUE_COLOR,'linewidth',3, 'LineStyle', 'none');
mean_val=(log10(code_mfa_mt_nadp_low)+log10(code_mfa_mt_nadp_high))/2; errorbar(p,mean_val,log10(code_mfa_mt_nadp_high)-mean_val,log10(code_mfa_mt_nadp_high)-mean_val,'color',RED_COLOR,'linewidth',3, 'LineStyle', 'none');

w=0.25;
%p=5; plot(polyshape([p-w p-w p+w p+w], log10([NADPH_con_min NADPH_con_max NADPH_con_max NADPH_con_min])),'FaceColor', RED_COLOR, 'FaceAlpha',0.7);
p=6; plot(polyshape([p-w p-w p+w p+w], log10([NADPH_con_min NADPH_con_max NADPH_con_max NADPH_con_min])),'FaceColor', GRAY_COLOR, 'FaceAlpha',0.7);
p=7; plot(polyshape([p-w p-w p+w p+w], log10([cy_nadph_low cy_nadph_high cy_nadph_high cy_nadph_low])),'FaceColor', BLUE_COLOR, 'FaceAlpha',0.7);
%mean_val=(log10(code_mfa_cy_nadph_low)+log10(code_mfa_cy_nadph_high))/2; errorbar(p,mean_val,log10(code_mfa_cy_nadph_high)-mean_val,log10(code_mfa_cy_nadph_high)-mean_val,'color',BLUE_COLOR,'linewidth',3, 'LineStyle', 'none');
mean_val=(log10(code_mfa_cy_nadph_low)+log10(code_mfa_cy_nadph_high))/2; errorbar(p,mean_val,log10(code_mfa_cy_nadph_high)-mean_val,log10(code_mfa_cy_nadph_high)-mean_val,'color',RED_COLOR,'linewidth',3, 'LineStyle', 'none');
p=8; plot(polyshape([p-w p-w p+w p+w], log10([mt_nadph_low mt_nadph_high mt_nadph_high mt_nadph_low])),'FaceColor', BLUE_COLOR, 'FaceAlpha',0.7);
%mean_val=(log10(code_mfa_mt_nadph_low)+log10(code_mfa_mt_nadph_high))/2; errorbar(p,mean_val,log10(code_mfa_mt_nadph_high)-mean_val,log10(code_mfa_mt_nadph_high)-mean_val,'color',BLUE_COLOR,'linewidth',3, 'LineStyle', 'none');
mean_val=(log10(code_mfa_mt_nadph_low)+log10(code_mfa_mt_nadph_high))/2; errorbar(p,mean_val,log10(code_mfa_mt_nadph_high)-mean_val,log10(code_mfa_mt_nadph_high)-mean_val,'color',RED_COLOR,'linewidth',3, 'LineStyle', 'none');



xlabel('');
ylabel('Concentrations [log10(mM)]');
xticks([1 2 3 6 7 8]);
xticklabels({'WC', 'CY', 'MT', 'WC', 'CY', 'MT'});
title('NADP+                          NADPH', 'Units', 'normalized', 'Position', [0.5, -0.14, 0]); % Set Title with correct Position
f=get(gca,'Children');
legend([f(10), f(9), f(8)],{' WC measured',' Thermodynamics',' CoDe-MFA'}, 'FontSize', 30, 'Location', 'northeast');
set(gcf,'color','w');
set(gca, 'FontSize', 30);
ylim([-5 1.5]);
yticks([-5:1:1]);
box on


set(gcf,'PaperSize',[50 30]);
s = sprintf('./output_images/1cd.pdf');
saveas(gcf, s);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% NAD/NADH NADP/NADPH in CY and MT based on thermodynamics, Code-MFA and
% various enzymes considered to be in chemical equilibrium and using WC
% level concentrations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure('units','normalized','outerposition',[0 0 1 1])

% lower and upper bounds of metabolites involved in the reactions that were
% used in the literature to determine co-factor ratios based on chemical
% equilibrium and WC concentrations
for(i=1:length(sensitiviy_analysis_concentration))
    if(strcmp(sensitiviy_analysis_concentration{i}.metabolite_name,'Pyruvate_CY'))
        code_mfa_cy_pyruvate_concentration = [exp(sensitiviy_analysis_concentration{i}.low_concentration) exp(sensitiviy_analysis_concentration{i}.high_concentration)];
    elseif(strcmp(sensitiviy_analysis_concentration{i}.metabolite_name,'Pyruvate_MT'))
        code_mfa_mt_pyruvate_concentration = [exp(sensitiviy_analysis_concentration{i}.low_concentration) exp(sensitiviy_analysis_concentration{i}.high_concentration)];        
    elseif(strcmp(sensitiviy_analysis_concentration{i}.metabolite_name,'Lactate_CY'))
        code_mfa_cy_lactate_concentration = [exp(sensitiviy_analysis_concentration{i}.low_concentration) exp(sensitiviy_analysis_concentration{i}.high_concentration)];
    elseif(strcmp(sensitiviy_analysis_concentration{i}.metabolite_name,'Lactate_MT'))
        code_mfa_mt_lactate_concentration = [exp(sensitiviy_analysis_concentration{i}.low_concentration) exp(sensitiviy_analysis_concentration{i}.high_concentration)];        
    elseif(strcmp(sensitiviy_analysis_concentration{i}.metabolite_name,'Malate_CY'))
        code_mfa_cy_malate_concentration = [exp(sensitiviy_analysis_concentration{i}.low_concentration) exp(sensitiviy_analysis_concentration{i}.high_concentration)];
    elseif(strcmp(sensitiviy_analysis_concentration{i}.metabolite_name,'Malate_MT'))
        code_mfa_mt_malate_concentration = [exp(sensitiviy_analysis_concentration{i}.low_concentration) exp(sensitiviy_analysis_concentration{i}.high_concentration)];        
    elseif(strcmp(sensitiviy_analysis_concentration{i}.metabolite_name,'Citrate_CY'))
        code_mfa_cy_citrate_concentration = [exp(sensitiviy_analysis_concentration{i}.low_concentration) exp(sensitiviy_analysis_concentration{i}.high_concentration)];
    elseif(strcmp(sensitiviy_analysis_concentration{i}.metabolite_name,'Citrate_MT'))
        code_mfa_mt_citrate_concentration = [exp(sensitiviy_analysis_concentration{i}.low_concentration) exp(sensitiviy_analysis_concentration{i}.high_concentration)];                
    elseif(strcmp(sensitiviy_analysis_concentration{i}.metabolite_name,'AKG_CY'))
        code_mfa_cy_akg_concentration = [exp(sensitiviy_analysis_concentration{i}.low_concentration) exp(sensitiviy_analysis_concentration{i}.high_concentration)];
    elseif(strcmp(sensitiviy_analysis_concentration{i}.metabolite_name,'AKG_MT'))
        code_mfa_mt_akg_concentration = [exp(sensitiviy_analysis_concentration{i}.low_concentration) exp(sensitiviy_analysis_concentration{i}.high_concentration)];        
    elseif(strcmp(sensitiviy_analysis_concentration{i}.metabolite_name,'OAA_CY'))
        code_mfa_cy_oaa_concentration = [exp(sensitiviy_analysis_concentration{i}.low_concentration) exp(sensitiviy_analysis_concentration{i}.high_concentration)];
    elseif(strcmp(sensitiviy_analysis_concentration{i}.metabolite_name,'OAA_MT'))
        code_mfa_mt_oaa_concentration = [exp(sensitiviy_analysis_concentration{i}.low_concentration) exp(sensitiviy_analysis_concentration{i}.high_concentration)];                
    end
end



% calculate CY NAD/NADH ratio assuming LDH in chemical equilibrium and Pyr & Lac
% WC measured concentrations reflect their cytosolic concentrations
cy_NAD_NADH_ratio_based_on_LDH_in_CE_and_WC_concentrations = exp((-G0_Pyruvate_to_Lactate(1))/RT - log(Lactate_con(1)) + log(Pyruvate_con(1)));
% find lower/upper Code-MFA dG values for LDH
index_rxns = find(contains(model_thermodynamics.full_rxns,'Pyruvate_CY + NADH_CY => Lactate_CY + NAD_CY'));
code_mfa_dG_Pyruvate_to_Lactate = [sensitiviy_analysis_dG{index_rxns}.low_dG sensitiviy_analysis_dG{index_rxns}.high_dG];
index_rxns  = find(contains(model_thermodynamics.full_rxns,'Pyruvate_CY + NADH_CY => Lactate_CY + NAD_CY'));

% calculate MT NAD/NADH ratio assuming malate dehydrogenase in chemical
% equilibrium and Malate & OAA WC measured concentrations reflect their mitochondrial concentrations
index_of_metabolite_in_wc_measured_vector = find(ismember(all_met_names_with_WC_measured_con,'OAA'));
OAA_con = all_met_WC_measured_con(index_of_metabolite_in_wc_measured_vector,:);
index_of_metabolite_in_wc_measured_vector = find(ismember(all_met_names_with_WC_measured_con,'Malate'));
Malate_con = all_met_WC_measured_con(index_of_metabolite_in_wc_measured_vector,:);
G0_Malate_to_OAA = [27.5 0.3];
mt_NADH_NAD_ratio_based_on_MDH_in_CE_and_WC_concentrations = exp((-G0_Malate_to_OAA(1))/RT - log(OAA_con(1)) + log(Malate_con(1)));
mt_NAD_NADH_ratio_based_on_MDH_in_CE_and_WC_concentrations = 1/mt_NADH_NAD_ratio_based_on_MDH_in_CE_and_WC_concentrations;
% find lower/upper Code-MFA dG values for MDH
index_rxns = find(contains(model_thermodynamics.full_rxns,'Malate_MT + NAD_MT => OAA_MT + NADH_MT'));
code_mfa_dG_Malate_to_OAA = [sensitiviy_analysis_dG{index_rxns}.low_dG sensitiviy_analysis_dG{index_rxns}.high_dG];


% calculate CY NADP/NADPH ratio assuming Malic enzyme in chemical equilibrium
% and Pyr & Mal WC measured concentrations reflect their cytosolic concentrations
G0_Malate_CY_to_Pyruvate_CY = -3.1;
cy_NADPH_NADP_ratio_based_on_ME_in_CE_and_WC_concentrations = exp((-G0_Malate_CY_to_Pyruvate_CY)/RT - log(Pyruvate_con(1)) - log(CO2_con(1)) + log(Malate_con(1)));
cy_NADP_NADPH_ratio_based_on_ME_in_CE_and_WC_concentrations = 1/cy_NADPH_NADP_ratio_based_on_ME_in_CE_and_WC_concentrations;
% find lower/upper Code-MFA dG values for Malic enzyme (ME1)
index_rxns = find(contains(model_thermodynamics.full_rxns,'Malate_CY + NADP_CY => Pyruvate_CY + CO2_sink + NADPH_CY'));
code_mfa_dG_Malate_CY_to_Pyruvate_CY = [sensitiviy_analysis_dG{index_rxns}.low_dG sensitiviy_analysis_dG{index_rxns}.high_dG];

% calculate CY NADP/NADPH ratio assuming isocitrate dehydrogenase in chemical equilibrium
% and aKG & Cit WC measured concentrations reflect their cytosolic concentrations
index_of_metabolite_in_wc_measured_vector = find(ismember(all_met_names_with_WC_measured_con,'Citrate'));
Citrate_con = all_met_WC_measured_con(index_of_metabolite_in_wc_measured_vector,:);
index_of_metabolite_in_wc_measured_vector = find(ismember(all_met_names_with_WC_measured_con,'AKG'));
aKG_con = all_met_WC_measured_con(index_of_metabolite_in_wc_measured_vector,:);
G0_Citrate_CY_to_aKG_CY = -3.1;
cy_NADPH_NADP_ratio_based_on_IDH1_in_CE_and_WC_concentrations = exp((-G0_Citrate_CY_to_aKG_CY)/RT - log(aKG_con(1)) - log(CO2_con(1)) + log(Citrate_con(1)));
cy_NADP_NADPH_ratio_based_on_IDH1_in_CE_and_WC_concentrations = 1/cy_NADPH_NADP_ratio_based_on_IDH1_in_CE_and_WC_concentrations;
% find lower/upper Code-MFA dG values for isocitrate dehydrogenase (IDH1)
index_rxns = find(contains(model_thermodynamics.full_rxns,'Citrate_CY + NADP_CY => AKG_CY + NADPH_CY + CO2_sink'));
code_mfa_dG_Citrate_CY_to_aKG_CY = [sensitiviy_analysis_dG{index_rxns}.low_dG sensitiviy_analysis_dG{index_rxns}.high_dG];

% calculate MT NADP/NADPH ratio assuming isocitrate dehydrogenase in chemical equilibrium
% and aKG & Cit WC measured concentrations reflect their cytosolic concentrations
G0_Citrate_MT_to_aKG_CY = -3.1;
mt_NADPH_NADP_ratio_based_on_IDH2_in_CE_and_WC_concentrations = exp((-G0_Citrate_MT_to_aKG_CY)/RT - log(aKG_con(1)) - log(CO2_con(1)) + log(Citrate_con(1)));
mt_NADP_NADPH_ratio_based_on_IDH2_in_CE_and_WC_concentrations = 1/mt_NADPH_NADP_ratio_based_on_IDH2_in_CE_and_WC_concentrations;
% find lower/upper Code-MFA dG values for isocitrate dehydrogenase (IDH2)
index_rxns = find(contains(model_thermodynamics.full_rxns,'AKG_MT + CO2_source + NADPH_MT => Citrate_MT + NADP_MT'));
code_mfa_dG_Citrate_MT_to_aKG_MT = [-sensitiviy_analysis_dG{index_rxns}.high_dG -sensitiviy_analysis_dG{index_rxns}.low_dG];

min_ylim = min([code_mfa_dG_Pyruvate_to_Lactate(1) code_mfa_dG_Malate_CY_to_Pyruvate_CY(1) code_mfa_dG_Citrate_CY_to_aKG_CY(1) code_mfa_dG_Citrate_MT_to_aKG_MT(1)])-1;
max_ylim = max([code_mfa_dG_Pyruvate_to_Lactate(2) code_mfa_dG_Malate_CY_to_Pyruvate_CY(2) code_mfa_dG_Citrate_CY_to_aKG_CY(2) code_mfa_dG_Citrate_MT_to_aKG_MT(2)])+1;

    
% figure for CY NAD/NADH ratio
% % % % % % % % % % % % % % % % % % % 
figure('units','normalized','outerposition',[0 0 1 1])
subplot(2,4,[1 2]);
grid on;
hold on;
 
xlim([0.5 3.5]);   
w=0.25;
p=1; plot(polyshape([p-w p-w p+w p+w], log10([cy_nad_nadh_ratio_low cy_nad_nadh_ratio_high cy_nad_nadh_ratio_high cy_nad_nadh_ratio_low])),'FaceColor', GREEN_COLOR, 'FaceAlpha',0.7);
mean_val=(log10(code_mfa_cy_nad_nadh_ratio_low)+log10(code_mfa_cy_nad_nadh_ratio_high))/2; errorbar(p,mean_val,log10(code_mfa_cy_nad_nadh_ratio_high)-mean_val,log10(code_mfa_cy_nad_nadh_ratio_high)-mean_val,'-', 'color',GREEN_COLOR,'linewidth',3, 'LineStyle', 'none');
plot(p,log10(cy_NAD_NADH_ratio_krebs),'o', 'MarkerSize',10, 'LineWidth', 1.5, 'Color', GREEN_COLOR);
plot(p,log10(cy_NAD_NADH_ratio_sun),'s', 'MarkerSize',10, 'LineWidth', 1.5, 'Color', GREEN_COLOR);
plot(p,log10(cy_NAD_NADH_ratio_hedeskov),'d', 'MarkerSize',10, 'LineWidth', 1.5, 'Color', GREEN_COLOR);
plot(p,log10(cy_NAD_NADH_ratio_based_on_LDH_in_CE_and_WC_concentrations),'*', 'color', GREEN_DARK_COLOR, 'MarkerSize',13, 'LineWidth', 2.5);
legend({'TD','Code-MFA', 'Krebs (Rat liver)', 'Sun (Hela)', 'hedeskov (Mouse pancreas)', 'TD (CE & WC)'}, 'FontSize', 15);
xlabel('');
ylabel('Cytosolic NAD+/NADH ratio [log10]');
xticks([1]);
xticklabels({'LDH'});
set(gcf,'color','w');
set(gca, 'FontSize', 14);
box on
hold off;

subplot(2,4,3);
grid on;
hold on;
xlim([0.5 1.5]);   
w=0.125;
p=1; mean_val=((code_mfa_dG_Pyruvate_to_Lactate(2))+(code_mfa_dG_Pyruvate_to_Lactate(1)))/2; errorbar(p,mean_val,code_mfa_dG_Pyruvate_to_Lactate(2)-mean_val,code_mfa_dG_Pyruvate_to_Lactate(2)-mean_val,'color',GREEN_COLOR,'linewidth',3, 'LineStyle', 'none');
xlabel('');
ylabel('Gibbs free energy [kJ/mol]');
xticks([1]);
xticklabels({'LDH (Pyr=>Lac)'});
set(gcf,'color','w');
set(gca, 'FontSize', 14);
box on
ylim([min_ylim max_ylim]); 
line([0.5 1.5], [0 0], 'color',RED_COLOR,'LineStyle',':','linewidth',2);
text(0.55,1.8,'Chemical','color',RED_COLOR, 'FontSize', 15);
text(0.55,-1.5,'equilibrium','color',RED_COLOR, 'FontSize', 15);
hold off;

subplot(2,4,4);
grid on;
hold on;
xlim([0 3]);     
gap=0.05;
p=1; mean_val=(log10(code_mfa_cy_pyruvate_concentration(2))+log10(code_mfa_cy_pyruvate_concentration(1)))/2; errorbar(p-gap,mean_val,log10(code_mfa_cy_pyruvate_concentration(2))-mean_val,log10(code_mfa_cy_pyruvate_concentration(2))-mean_val,'color',GREEN_DARK_COLOR,'linewidth',2, 'LineStyle', 'none');
text(p-gap-0.35,log10(code_mfa_cy_pyruvate_concentration(2)),'CY','color',GREEN_DARK_COLOR, 'FontSize', 12);
p=1; mean_val=(log10(code_mfa_mt_pyruvate_concentration(2))+log10(code_mfa_mt_pyruvate_concentration(1)))/2; errorbar(p+gap,mean_val,log10(code_mfa_mt_pyruvate_concentration(2))-mean_val,log10(code_mfa_mt_pyruvate_concentration(2))-mean_val,'color',GREEN_COLOR,'linewidth',2, 'LineStyle', 'none');
text(p+gap+0.1,log10(code_mfa_mt_pyruvate_concentration(2)),'MT','color',GREEN_COLOR, 'FontSize', 12);
plot(p,log10(Pyruvate_con(1)),'*','color',[0 0.2 0]);
text(p+0.1,log10(Pyruvate_con(1)),'WC','color',[0 0.2 0], 'FontSize', 12);
p=2; mean_val=(log10(code_mfa_cy_lactate_concentration(2))+log10(code_mfa_cy_lactate_concentration(1)))/2; errorbar(p-gap,mean_val,log10(code_mfa_cy_lactate_concentration(2))-mean_val,log10(code_mfa_cy_lactate_concentration(2))-mean_val,'color',GREEN_DARK_COLOR,'linewidth',2, 'LineStyle', 'none');
text(p-gap-0.35,log10(code_mfa_cy_lactate_concentration(2)),'CY','color',GREEN_DARK_COLOR, 'FontSize', 12);
plot(p,log10(Lactate_con(1)),'*','color',[0 0.2 0]);
text(p+0.1,log10(Lactate_con(1)),'WC','color',[0 0.2 0], 'FontSize', 12);
xlabel('');
ylabel('Concentration [mM(log10)]');
xticks([1 2]);
xticklabels({'Pyr', 'Lac'});
set(gcf,'color','w');
set(gca, 'FontSize', 14);
ylim([-5 1]);
box on
hold off;


% figure for MT NAD/NADH ratio
% % % % % % % % % % % % % % % % % % % 
subplot(2,4,[5 6]);
grid on;
hold on;

xlim([0.5 3.5]);   
w=0.125;
p=1; plot(polyshape([p-w p-w p+w p+w], log10([mt_nad_nadh_ratio_low mt_nad_nadh_ratio_high mt_nad_nadh_ratio_high mt_nad_nadh_ratio_low])),'FaceColor', GREEN_COLOR, 'FaceAlpha',0.7);
mean_val=(log10(code_mfa_mt_nad_nadh_ratio_low)+log10(code_mfa_mt_nad_nadh_ratio_high))/2; errorbar(p,mean_val,log10(code_mfa_mt_nad_nadh_ratio_high)-mean_val,log10(code_mfa_mt_nad_nadh_ratio_high)-mean_val,'color',GREEN_COLOR,'linewidth',3, 'LineStyle', 'none');
plot(p,log10(mt_NAD_NADH_ratio_won), 'x', 'MarkerSize',10, 'LineWidth', 1.5, 'Color', GREEN_COLOR);
plot(p,log10(mt_NAD_NADH_ratio_sabatini), '+', 'MarkerSize',10, 'LineWidth', 1.5, 'Color', GREEN_COLOR);
plot(p,log10(mt_NAD_NADH_ratio_krebs), 'o', 'MarkerSize',10, 'LineWidth', 1.5, 'Color', GREEN_COLOR);
plot(p,log10(mt_NAD_NADH_ratio_elmar), '^', 'MarkerSize',10, 'LineWidth', 1.5, 'Color', GREEN_COLOR);
plot(p,log10(mt_NAD_NADH_ratio_based_on_MDH_in_CE_and_WC_concentrations),'*', 'color', GREEN_DARK_COLOR, 'MarkerSize',13, 'LineWidth', 2.5);

legend({'TD','Code-MFA', 'Won (Hela)', 'Sabatini (Hela)', 'Krebs (Rat liver)', 'Elmar (Rat liver)', 'TD (CE & WC)'}, 'FontSize', 15);

xlabel('');
ylabel('Mitochondrial NAD+/NADH ratio [log10]');
xticks([1]);
xticklabels({'MDH'});
set(gcf,'color','w');
set(gca, 'FontSize', 14);
box on
hold off;

subplot(2,4,7);
grid on;
hold on;
xlim([0.5 1.5]);   
w=0.125;
p=1; mean_val=((code_mfa_dG_Malate_to_OAA(2))+(code_mfa_dG_Malate_to_OAA(1)))/2; errorbar(p,mean_val,code_mfa_dG_Malate_to_OAA(2)-mean_val,code_mfa_dG_Malate_to_OAA(2)-mean_val,'color',GREEN_COLOR,'linewidth',3, 'LineStyle', 'none');
xlabel('');
ylabel('Gibbs free energy [kJ/mol]');
xticks([1]);
xticklabels({'MDH (Mal=>OAA)'});
set(gcf,'color','w');
set(gca, 'FontSize', 14);
box on
ylim([min_ylim max_ylim]); 
line([0.5 1.5], [0 0], 'color',RED_COLOR,'LineStyle',':','linewidth',2);
text(0.55,1.8,'Chemical','color',RED_COLOR, 'FontSize', 15);
text(0.55,-1.5,'equilibrium','color',RED_COLOR, 'FontSize', 15);
hold off;

subplot(2,4,8);
grid on;
hold on;
xlim([0 3]);     
gap=0.05;
p=1; mean_val=(log10(code_mfa_cy_malate_concentration(2))+log10(code_mfa_cy_malate_concentration(1)))/2; errorbar(p-gap,mean_val,log10(code_mfa_cy_malate_concentration(2))-mean_val,log10(code_mfa_cy_malate_concentration(2))-mean_val,'color',GREEN_DARK_COLOR,'linewidth',2, 'LineStyle', 'none');
text(p-gap-0.35,log10(code_mfa_cy_malate_concentration(2)),'CY','color',GREEN_DARK_COLOR, 'FontSize', 12);
p=1; mean_val=(log10(code_mfa_mt_malate_concentration(2))+log10(code_mfa_mt_malate_concentration(1)))/2; errorbar(p+gap,mean_val,log10(code_mfa_mt_malate_concentration(2))-mean_val,log10(code_mfa_mt_malate_concentration(2))-mean_val,'color',GREEN_COLOR,'linewidth',2, 'LineStyle', 'none');
text(p+gap+0.1,log10(code_mfa_mt_malate_concentration(2)),'MT','color',GREEN_COLOR, 'FontSize', 12);
plot(p,log10(Malate_con(1)),'*','color',[0 0.2 0]);
text(p+0.1,log10(Malate_con(1)),'WC','color',[0 0.2 0], 'FontSize', 12);
p=2; mean_val=(log10(code_mfa_cy_oaa_concentration(2))+log10(code_mfa_cy_oaa_concentration(1)))/2; errorbar(p-gap,mean_val,log10(code_mfa_cy_oaa_concentration(2))-mean_val,log10(code_mfa_cy_oaa_concentration(2))-mean_val,'color',GREEN_DARK_COLOR,'linewidth',2, 'LineStyle', 'none');
text(p-gap-0.35,log10(code_mfa_cy_oaa_concentration(2)),'CY','color',GREEN_DARK_COLOR, 'FontSize', 12);
p=2; mean_val=(log10(code_mfa_mt_oaa_concentration(2))+log10(code_mfa_mt_oaa_concentration(1)))/2; errorbar(p+gap,mean_val,log10(code_mfa_mt_oaa_concentration(2))-mean_val,log10(code_mfa_mt_oaa_concentration(2))-mean_val,'color',GREEN_COLOR,'linewidth',2, 'LineStyle', 'none');
text(p+gap+0.1,log10(code_mfa_mt_oaa_concentration(2)),'MT','color',GREEN_COLOR, 'FontSize', 12);
plot(p,log10(OAA_con(1)),'*','color',[0 0.2 0]);
text(p+0.1,log10(OAA_con(1)),'WC','color',[0 0.2 0], 'FontSize', 12);
xlabel('');
ylabel('Concentration [mM(log10)]');
xticks([1 2]);
xticklabels({'Mal', 'OAA'});
set(gcf,'color','w');
set(gca, 'FontSize', 14);
ylim([-5 1]);
box on
hold off;




% figure for CY NADP/NADPH ratio
% % % % % % % % % % % % % % % % % % % 
cy_NADP_NADPH_ratio_krebs_ME1   = 0.0118;
cy_NADP_NADPH_ratio_krebs_IDH1  = 0.0101;


figure('units','normalized','outerposition',[0 0 1 1])
subplot(2,4,[1 2]);
grid on;
hold on;

xlim([0.5 3.5]);   
w=0.25;
p=1; plot(polyshape([p-w p-w p+w p+w], log10([cy_nadp_nadph_ratio_low cy_nadp_nadph_ratio_high cy_nadp_nadph_ratio_high cy_nadp_nadph_ratio_low])),'FaceColor', GREEN_COLOR);
mean_val=(log10(code_mfa_cy_nadp_nadph_ratio_low)+log10(code_mfa_cy_nadp_nadph_ratio_high))/2; errorbar(p,mean_val,log10(code_mfa_cy_nadp_nadph_ratio_high)-mean_val,log10(code_mfa_cy_nadp_nadph_ratio_high)-mean_val,'color',GREEN_COLOR,'linewidth',3, 'LineStyle', 'none');
plot(p,log10(cy_NADP_NADPH_ratio_hedeskov), 'o', 'MarkerSize',10, 'LineWidth', 1.5, 'Color', GREEN_COLOR);
plot(p,log10(cy_NADP_NADPH_ratio_krebs_ME1), 's', 'MarkerSize',10, 'LineWidth', 1.5, 'Color', GREEN_COLOR);
plot(p,log10(cy_NADP_NADPH_ratio_elmar), 'd', 'MarkerSize',10, 'LineWidth', 1.5, 'Color', GREEN_COLOR);
plot(p,log10(cy_NADP_NADPH_ratio_based_on_ME_in_CE_and_WC_concentrations),'*', 'color', GREEN_DARK_COLOR, 'MarkerSize',13, 'LineWidth', 2.5);
p=1.75; plot(polyshape([p-w p-w p+w p+w], log10([cy_nadp_nadph_ratio_low cy_nadp_nadph_ratio_high cy_nadp_nadph_ratio_high cy_nadp_nadph_ratio_low])),'FaceColor', BLUE_COLOR, 'FaceAlpha',0.7);
mean_val=(log10(code_mfa_cy_nadp_nadph_ratio_low)+log10(code_mfa_cy_nadp_nadph_ratio_high))/2; errorbar(p,mean_val,log10(code_mfa_cy_nadp_nadph_ratio_high)-mean_val,log10(code_mfa_cy_nadp_nadph_ratio_high)-mean_val,'color',BLUE_COLOR,'linewidth',3, 'LineStyle', 'none');
plot(p,log10(cy_NADP_NADPH_ratio_krebs_IDH1),'s', 'MarkerSize',10, 'LineWidth', 1.5, 'Color', BLUE_COLOR);
plot(p,log10(cy_NADP_NADPH_ratio_based_on_IDH1_in_CE_and_WC_concentrations),'*', 'color', BLUE_DARK_COLOR, 'MarkerSize',13, 'LineWidth', 2.5);

legend({'TD','Code-MFA', 'hedeskov (Mouse pancreas)', 'Krebs (Rat liver)', 'Elmar (Mouse pancreas)', 'TD (CE & WC)', 'TD','Code-MFA', 'Krebs (Rat liver)', 'TD (CE & WC)'}, 'FontSize', 15);
xlabel('');
ylabel('Cytosolic NADP+/NADPH ratio [log10]');
xticks([1 1.75]);
xticklabels({'ME1', 'IDH1'});
set(gcf,'color','w');
set(gca, 'FontSize', 14);
box on
hold off;

subplot(2,4,3);
grid on;
hold on;
xlim([0.5 2.5]);   
w=0.25;
p=1; mean_val=((code_mfa_dG_Malate_CY_to_Pyruvate_CY(2))+(code_mfa_dG_Malate_CY_to_Pyruvate_CY(1)))/2; errorbar(p,mean_val,code_mfa_dG_Malate_CY_to_Pyruvate_CY(2)-mean_val,code_mfa_dG_Malate_CY_to_Pyruvate_CY(2)-mean_val,'color',GREEN_COLOR,'linewidth',3, 'LineStyle', 'none');
p=2; mean_val=((code_mfa_dG_Citrate_CY_to_aKG_CY(2))+(code_mfa_dG_Citrate_CY_to_aKG_CY(1)))/2; errorbar(p,mean_val,code_mfa_dG_Citrate_CY_to_aKG_CY(2)-mean_val,code_mfa_dG_Citrate_CY_to_aKG_CY(2)-mean_val,'color',BLUE_COLOR,'linewidth',3, 'LineStyle', 'none');
xlabel('');
ylabel('Gibbs free energy [kJ/mol]');
xticks([1 2]);
xticklabels({'ME1 (Mal=>Pyr)', 'IDH1 (Cit=>aKG)'});
set(gcf,'color','w');
set(gca, 'FontSize', 14);
box on
ylim([min_ylim max_ylim]); 
line([0.5 2.5], [0 0], 'color',RED_COLOR,'LineStyle',':','linewidth',2);
text(0.6,1.8,'Chemical','color',RED_COLOR, 'FontSize', 15);
text(0.6,-1.5,'equilibrium','color',RED_COLOR, 'FontSize', 15);
hold off;

subplot(2,4,4);
grid on;
hold on;
xlim([0 6]);     
gap=0.125;
p=1; mean_val=(log10(code_mfa_cy_malate_concentration(2))+log10(code_mfa_cy_malate_concentration(1)))/2; errorbar(p-gap,mean_val,log10(code_mfa_cy_malate_concentration(2))-mean_val,log10(code_mfa_cy_malate_concentration(2))-mean_val,'color',GREEN_DARK_COLOR,'linewidth',2, 'LineStyle', 'none');
text(p-gap-0.65,log10(code_mfa_cy_malate_concentration(2)),'CY','color',GREEN_DARK_COLOR, 'FontSize', 12);
p=1; mean_val=(log10(code_mfa_mt_malate_concentration(2))+log10(code_mfa_mt_malate_concentration(1)))/2; errorbar(p+gap,mean_val,log10(code_mfa_mt_malate_concentration(2))-mean_val,log10(code_mfa_mt_malate_concentration(2))-mean_val,'color',GREEN_COLOR,'linewidth',2, 'LineStyle', 'none');
text(p+gap+0.18,log10(code_mfa_mt_malate_concentration(2)),'MT','color',GREEN_COLOR, 'FontSize', 12);
plot(p,log10(Malate_con(1)),'*','color',[0 0.2 0]);
text(p+0.18,log10(Malate_con(1)),'WC','color',[0 0.2 0], 'FontSize', 12);
p=2; mean_val=(log10(code_mfa_cy_pyruvate_concentration(2))+log10(code_mfa_cy_pyruvate_concentration(1)))/2; errorbar(p-gap,mean_val,log10(code_mfa_cy_pyruvate_concentration(2))-mean_val,log10(code_mfa_cy_pyruvate_concentration(2))-mean_val,'color',GREEN_DARK_COLOR,'linewidth',2, 'LineStyle', 'none');
text(p-gap+0.13,0.1+log10(code_mfa_cy_pyruvate_concentration(2)),'CY','color',GREEN_DARK_COLOR, 'FontSize', 12);
p=2; mean_val=(log10(code_mfa_mt_pyruvate_concentration(2))+log10(code_mfa_mt_pyruvate_concentration(1)))/2; errorbar(p+gap,mean_val,log10(code_mfa_mt_pyruvate_concentration(2))-mean_val,log10(code_mfa_mt_pyruvate_concentration(2))-mean_val,'color',GREEN_COLOR,'linewidth',2, 'LineStyle', 'none');
text(p+gap+0.18,log10(code_mfa_mt_pyruvate_concentration(2)),'MT','color',GREEN_COLOR, 'FontSize', 12);
plot(p,log10(Pyruvate_con(1)),'*','color',[0 0.2 0]);
text(p+0.18,log10(Pyruvate_con(1)),'WC','color',[0 0.2 0], 'FontSize', 12);
p=4; mean_val=(log10(code_mfa_cy_citrate_concentration(2))+log10(code_mfa_cy_citrate_concentration(1)))/2; errorbar(p-gap,mean_val,log10(code_mfa_cy_citrate_concentration(2))-mean_val,log10(code_mfa_cy_citrate_concentration(2))-mean_val,'color',BLUE_DARK_COLOR,'linewidth',2, 'LineStyle', 'none');
text(p-gap-0.65,log10(code_mfa_cy_citrate_concentration(2)),'CY','color',BLUE_DARK_COLOR, 'FontSize', 12);
p=4; mean_val=(log10(code_mfa_mt_citrate_concentration(2))+log10(code_mfa_mt_citrate_concentration(1)))/2; errorbar(p+gap,mean_val,log10(code_mfa_mt_citrate_concentration(2))-mean_val,log10(code_mfa_mt_citrate_concentration(2))-mean_val,'color',BLUE_COLOR,'linewidth',2, 'LineStyle', 'none');
text(p+gap+0.18,log10(code_mfa_mt_citrate_concentration(2)),'MT','color',BLUE_COLOR, 'FontSize', 12);
plot(p,log10(Citrate_con(1)),'*','color',[0 0 0.2]);
text(p+0.18,log10(Citrate_con(1)),'WC','color',[0 0 0.2], 'FontSize', 12);
p=5; mean_val=(log10(code_mfa_cy_akg_concentration(2))+log10(code_mfa_cy_akg_concentration(1)))/2; errorbar(p-gap,mean_val,log10(code_mfa_cy_akg_concentration(2))-mean_val,log10(code_mfa_cy_akg_concentration(2))-mean_val,'color',BLUE_DARK_COLOR,'linewidth',2, 'LineStyle', 'none');
text(p-gap-0.65,log10(code_mfa_cy_akg_concentration(2)),'CY','color',BLUE_DARK_COLOR, 'FontSize', 12);
p=5; mean_val=(log10(code_mfa_mt_akg_concentration(2))+log10(code_mfa_mt_akg_concentration(1)))/2; errorbar(p+gap,mean_val,log10(code_mfa_mt_akg_concentration(2))-mean_val,log10(code_mfa_mt_akg_concentration(2))-mean_val,'color',BLUE_COLOR,'linewidth',2, 'LineStyle', 'none');
text(p+gap+0.18,log10(code_mfa_mt_akg_concentration(2)),'MT','color',BLUE_COLOR, 'FontSize', 12);
plot(p,log10(aKG_con(1)),'*','color',[0 0 0.2]);
text(p+0.18,log10(aKG_con(1)),'WC','color',[0 0 0.2], 'FontSize', 12);
xlabel('');
ylabel('Concentration [mM(log10)]');
xticks([1 2 4 5]);
xticklabels({'Mal', 'Pyr', 'Cit', 'AKG'});
set(gcf,'color','w');
set(gca, 'FontSize', 14);
ylim([-5 1]);
box on
hold off;



% figure for MT NADP/NADPH ratio
% % % % % % % % % % % % % % % % % % % 
subplot(2,4,[5 6]);
grid on;
hold on;

xlim([0.5 3.5]);   
w=0.125;
p=1; plot(polyshape([p-w p-w p+w p+w], log10([mt_nadp_nadph_ratio_low mt_nadp_nadph_ratio_high mt_nadp_nadph_ratio_high mt_nadp_nadph_ratio_low])),'FaceColor', GREEN_COLOR, 'FaceAlpha',0.7);
mean_val=(log10(code_mfa_mt_nadp_nadph_ratio_low)+log10(code_mfa_mt_nadp_nadph_ratio_high))/2; errorbar(p,mean_val,log10(code_mfa_mt_nadp_nadph_ratio_high)-mean_val,log10(code_mfa_mt_nadp_nadph_ratio_high)-mean_val,'color',GREEN_COLOR,'linewidth',3, 'LineStyle', 'none');
plot(p,log10(mt_NADP_NADPH_ratio_sies), 'o', 'MarkerSize',10, 'LineWidth', 1.5, 'Color', GREEN_COLOR);
plot(p,log10(mt_NADP_NADPH_ratio_based_on_IDH2_in_CE_and_WC_concentrations),'*', 'color', GREEN_DARK_COLOR, 'MarkerSize',13, 'LineWidth', 2.5);
                                            
legend({'TD','Code-MFA', 'Sies (Rat Hepatocytes)', 'TD (CE & WC)'}, 'FontSize', 15);
xlabel('');
ylabel('Mitochondrial NADP+/NADPH ratio [log10]');
xticks([1]);
xticklabels({'IDH2'});
set(gcf,'color','w');
set(gca, 'FontSize', 14);
box on
hold off;

subplot(2,4,7);
grid on;
hold on;
xlim([0.5 1.5]);   
w=0.125;
p=1; mean_val=((code_mfa_dG_Citrate_MT_to_aKG_MT(2))+(code_mfa_dG_Citrate_MT_to_aKG_MT(1)))/2; errorbar(p,mean_val,code_mfa_dG_Citrate_MT_to_aKG_MT(2)-mean_val,code_mfa_dG_Citrate_MT_to_aKG_MT(2)-mean_val,'color',GREEN_COLOR,'linewidth',3, 'LineStyle', 'none');
xlabel('');
ylabel('Gibbs free energy [kJ/mol]');
xticks([1]);
xticklabels({'IDH2 (Cit=>aKG)'});
set(gcf,'color','w');
set(gca, 'FontSize', 14);
box on
ylim([min_ylim max_ylim]); 
line([0.5 1.5], [0 0], 'color',RED_COLOR,'LineStyle',':','linewidth',2);
text(0.55,1.8,'Chemical','color',RED_COLOR, 'FontSize', 15);
text(0.55,-1.5,'equilibrium','color',RED_COLOR, 'FontSize', 15);
hold off;

subplot(2,4,8);
grid on;
hold on;
xlim([0 3]);     
gap=0.05;
p=1; mean_val=(log10(code_mfa_cy_citrate_concentration(2))+log10(code_mfa_cy_citrate_concentration(1)))/2; errorbar(p-gap,mean_val,log10(code_mfa_cy_citrate_concentration(2))-mean_val,log10(code_mfa_cy_citrate_concentration(2))-mean_val,'color',GREEN_DARK_COLOR,'linewidth',2, 'LineStyle', 'none');
text(p-gap-0.35,log10(code_mfa_cy_citrate_concentration(2)),'CY','color',GREEN_DARK_COLOR, 'FontSize', 12);
p=1; mean_val=(log10(code_mfa_mt_citrate_concentration(2))+log10(code_mfa_mt_citrate_concentration(1)))/2; errorbar(p+gap,mean_val,log10(code_mfa_mt_citrate_concentration(2))-mean_val,log10(code_mfa_mt_citrate_concentration(2))-mean_val,'color',GREEN_COLOR,'linewidth',2, 'LineStyle', 'none');
text(p+gap+0.1,log10(code_mfa_mt_citrate_concentration(2)),'MT','color',GREEN_COLOR, 'FontSize', 12);
plot(p,log10(Citrate_con(1)),'*','color',[0 0.2 0]);
text(p+0.1,log10(Citrate_con(1)),'WC','color',[0 0.2 0], 'FontSize', 12);
p=2; mean_val=(log10(code_mfa_cy_akg_concentration(2))+log10(code_mfa_cy_akg_concentration(1)))/2; errorbar(p-gap,mean_val,log10(code_mfa_cy_akg_concentration(2))-mean_val,log10(code_mfa_cy_akg_concentration(2))-mean_val,'color',GREEN_DARK_COLOR,'linewidth',2, 'LineStyle', 'none');
text(p-gap-0.35,log10(code_mfa_cy_akg_concentration(2)),'CY','color',GREEN_DARK_COLOR, 'FontSize', 12);
p=2; mean_val=(log10(code_mfa_mt_akg_concentration(2))+log10(code_mfa_mt_akg_concentration(1)))/2; errorbar(p+gap,mean_val,log10(code_mfa_mt_akg_concentration(2))-mean_val,log10(code_mfa_mt_akg_concentration(2))-mean_val,'color',GREEN_COLOR,'linewidth',2, 'LineStyle', 'none');
text(p+gap+0.1,log10(code_mfa_mt_akg_concentration(2)),'MT','color',GREEN_COLOR, 'FontSize', 12);
plot(p,log10(aKG_con(1)),'*','color',[0 0.2 0]);
text(p+0.1,log10(aKG_con(1)),'WC','color',[0 0.2 0], 'FontSize', 12);
xlabel('');
ylabel('Concentration [mM(log10)]');
xticks([1 2]);
xticklabels({'Cit', 'AKG'});
set(gcf,'color','w');
set(gca, 'FontSize', 14);
ylim([-5 1]);
box on
hold off;


