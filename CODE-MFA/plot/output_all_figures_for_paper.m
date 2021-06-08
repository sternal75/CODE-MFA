% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% Figure 1 - Inferring mitochondrial and cytosolic co-factor concentrations 
% and ratios in HeLa cells. 
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 

% (b) The steady state fractional M+5 labeling of 6-phosphogluconate (6PG) 
% and ribulose-5-phonehate (R5P) when feeding [U-13C]-glucose
output_R5P_and_6PG_MID_from_glc_labeling                   % 1b.pdf
% (c-d) The measured NADP+/NADPH ratio and NADP+ and NADPH concentrations (gray) 
% and inferred compartmentalized ratio and concentrations based on simple 
% thermodynamics and deconvolution analysis (blue), and via CODE-MFA (red). 
% Asterisks represent previously published ratios (Table S12). 
% (g-h) The measured NAD+/NADH ratio and NAD+ and NADH concentrations (gray) 
% and inferred compartmentalized values based on simple thermodynamics and 
% deconvolution analysis (blue), and via CODE-MFA (red)
output_nad_nadh_nadp_nadph                                 % 1cdgh.pdf
% (f) The fractional M+3 labeling of lactate and pyruvate when feeding [U-13C]-lactate
output_lactate_and_pyruvate_MID_fron_lactate_labeling      % 1f.pdf

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% Figure 3 - CODE-MFA performance inferring compartmentalized fluxes, 
% Gibbs energies, and concentration in HeLa cells
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 

% (a) Number of reactions whose direction of net flux is uniquely inferred across CODE-MFA iterations. 
output_known_directionalities_per_iteration                % 3a.pdf
% (b) The fit been simulated total cellular metabolite concentration 
% (i.e. convolution of simulated mitochondrial and cytosolic labeling; 
% x-axis) and measurements (y-axis) 
output_measured_vs_computed_concentrations                 % 3b.pdf
% (c) The fit been simulated total cellular metabolite isotopic labeling (x-axis) 
% and experimental measurements (y-axis)  
output_measured_vs_computed_labeling_forms                 % 3c.pdf
% (d) Total percentage of reactions in the model whose direction of net flux 
% is inferred by CODE-MFA versus with MFA and CODE-MFA without thermodynamic considerations
% (e) Cumulative distribution of reaction net flux confidence interval sizes 
% inferred by CODE-MFA versus MFA and CODE-MFA without thermodynamic considerations 
output_compare_fluxes_to_mfa                               % 3de.pdf
% (f) Correlation between gene expression levels and most probable net flux 
% inferred by CODE-MFA through the corresponding enzyme
% (g) Correlation between gene expression levels and most probable net flux 
% inferred by MFA through the corresponding enzyme
output_gene_expression_vs_fluxes                           % 3fg.pdf
% (h) Cumulative distribution of reaction metabolite concentration confidence 
% interval sizes inferred by CODE-MFA versus with strictly thermodynamic analysis 
output_comulative_concentrations                           % 3h.pdf
% (i) Cumulative distribution of reaction Gibbs energy confidence interval 
% sizes inferred by CODE-MFA versus with strictly thermodynamic analysis
output_comulative_dG                                       % 3i.pdf

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% Figure 4 - CODE-MFA derived compartmentalized fluxes, Gibbs energies, 
% and concentration confidence intervals in HeLa cells
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 

% CODE-MFA derived Gibbs free energies and net fluxes for cytosolic and 
% mitochondrial isozyme (a), for strictly cytosolic reactions (b), 
% and for strictly mitochondrial reactions (c). 
% Blue and green bars represent cytosolic and mitochondrial fluxes and 
% Gibbs energies, respectively; asterisks represent Gibbs energies computed 
% directly based on measured cellular metabolite concentrations
output_per_compartment_dG_and_fluxes_result                % 4abc
% (d) CODE-MFA derived cytosolic and mitochondrial metabolite concentrations 
% (blue and green bars) and measured cellular concentrations (asterisk)
output_per_concentration_result                            % 4d.pdf

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% Figure 6 - Correlation of CODE-MFA derived fluxes for HeLa, 
% HCT116, A549 and LN229 cell lines
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 

% (a) Pearson correlation of CODE-MFA derived fluxes for HeLa and A549 
% (b) Pearson correlation of CODE-MFA derived fluxes for HeLa and HCT116
% (c) Pearson correlation of CODE-MFA derived fluxes for HeLa and LN229 
% (d) Pearson correlation of CODE-MFA derived fluxes for A549 and HCT116 
% (e) Pearson correlation of CODE-MFA derived fluxes for A549 and LN229 
% (f) Pearson correlation of CODE-MFA derived fluxes for HCT116 and LN229 
% (g) CODE-MFA derived net fluxes for fluxes with large difference 
% per figure f (HCT116 vs. LN229)
% output_correlation_all_cell_lines                          % 6abcdefg.pdf

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% Figure S1 - Gene expression levels obtained from the 
% Broad Institute Cancer Cell Line Encyclopedia (CCLE), with 0.5 RPKM value 
% as the cutoff to determine non-expressed genes
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
output_gene_expression_levels                              % s1.pdf

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% Figure S2 - Correlation between measured enzyme abundance in HeLa44 (x-axis) 
% and most probable net fluxes (y-axis) derived by CODE-MFA (a) and MFA (b)
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
output_enzyme_concentrations_vs_fluxes                     % s2ab.pdf



