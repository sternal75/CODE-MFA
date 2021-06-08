% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% Correlation between gene expression levels and most probable net flux 
% inferred by MFA through the corresponding enzyme
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 

clear all;
close all;

run ../load_constants;

load('../mat_files/directionalities.mat', 'directionalities');
load('../mat_files/model_thermodynamics.mat','model_thermodynamics');    
load('../mat_files/net_fluxes.mat', 'net_fluxes');


[dFluxes,tFluxes] = xlsread('../xls_input_files/mfa_and_code_mfa_without_thermodynamics_results.xlsx');
low_flux_MFA    = dFluxes(:,1)';
high_flux_MFA   = dFluxes(:,2)';


% load enzyme names and full reactions
[dEnzymes,tEnzymes] = xlsread('../xls_input_files/enzyme_gene_expressions.xlsx','enzymes');
tEnzymes = tEnzymes(2:end,:);

% get best score net fluxes
index_best_score = find(directionalities.errors==min(directionalities.errors));    
best_score_predicted_net_fluxes = directionalities.predicted_net_fluxes(:,index_best_score);

% Fluxes of isoenzymes computed by Code-MFA
enzymes         = cell(0);
best_fit_fluxes = [];
mean_fluxes     = [];
mean_fluxes_MFA = [];
enzyme_gene_expression = [];
for(i=1:size(tEnzymes,1))
    if((~isempty(tEnzymes{i,2}))&&(dEnzymes(i,2)~=1))
        full_reaction           = tEnzymes{i,1};
        
        % if the gene encodes more than one reactions in my model, sum all
        % the fluxes of all reactions encoded by this gene
        existint_index = find(ismember(enzymes,tEnzymes{i,2}));
        if(~isempty(existint_index))
            % find reaction index
            reaction_index = find(ismember(model_thermodynamics.full_rxns,full_reaction));
            best_fit_fluxes(existint_index)    = best_fit_fluxes(existint_index)+abs(best_score_predicted_net_fluxes(reaction_index(1)));
            mean_fluxes(existint_index)        = mean_fluxes(existint_index)+abs((net_fluxes{reaction_index(1)}.low_flux+net_fluxes{reaction_index(1)}.high_flux)/2);
            mean_fluxes_MFA(existint_index)    = mean_fluxes_MFA(existint_index)+abs((low_flux_MFA(reaction_index(1))+high_flux_MFA(reaction_index(1)))/2);
            if(mean_fluxes_MFA(existint_index)==0)
                mean_fluxes_MFA(existint_index) = 0.01;
            end
        % if it is the first time we get this gene
        else
            enzymes{end+1}          = tEnzymes{i,2};
            % enzyme gene expression
            enzyme_gene_expression(end+1) = dEnzymes(i,1);
            % find reaction index
            reaction_index = find(ismember(model_thermodynamics.full_rxns,full_reaction));
            best_fit_fluxes(end+1)    = abs(best_score_predicted_net_fluxes(reaction_index(1)));
            mean_fluxes(end+1)        = abs((net_fluxes{reaction_index(1)}.low_flux+net_fluxes{reaction_index(1)}.high_flux)/2);
            mean_fluxes_MFA(end+1)    = abs((low_flux_MFA(reaction_index(1))+high_flux_MFA(reaction_index(1)))/2);            
            if(mean_fluxes_MFA(end)==0)
                mean_fluxes_MFA(end) = 0.01;
            end            
        end
    end
end

gene_name='GOT1';
index           = find(ismember(enzymes,gene_name));
GOT1_best_fit_flux_abs   = abs(best_fit_fluxes(index));
GOT1_mean_flux_abs   = abs(mean_fluxes(index));
GOT1_gene_expression = enzyme_gene_expression(index);

gene_name='GOT2';
index           = find(ismember(enzymes,gene_name));
GOT2_best_fit_flux_abs   = abs(best_fit_fluxes(index));
GOT2_mean_flux_abs   = abs(mean_fluxes(index));
GOT2_gene_expression = enzyme_gene_expression(index);

gene_name='MDH1';
index           = find(ismember(enzymes,gene_name));
MDH1_best_fit_flux_abs   = abs(best_fit_fluxes(index));
MDH1_mean_flux_abs   = abs(mean_fluxes(index));
MDH1_gene_expression = enzyme_gene_expression(index);

gene_name='MDH2';
index           = find(ismember(enzymes,gene_name));
MDH2_best_fit_flux_abs   = abs(best_fit_fluxes(index));
MDH2_mean_flux_abs   = abs(mean_fluxes(index));
MDH2_gene_expression = enzyme_gene_expression(index);

gene_name='ME1';
index           = find(ismember(enzymes,gene_name));
ME1_best_fit_flux_abs   = abs(best_fit_fluxes(index));
ME1_mean_flux_abs   = abs(mean_fluxes(index));
ME1_gene_expression = enzyme_gene_expression(index);

gene_name='ME3';
index           = find(ismember(enzymes,gene_name));
ME3_best_fit_flux_abs   = abs(best_fit_fluxes(index));
ME3_mean_flux_abs   = abs(mean_fluxes(index));
ME3_gene_expression = enzyme_gene_expression(index);

gene_name='IDH1';
index           = find(ismember(enzymes,gene_name));
IDH1_best_fit_flux_abs   = abs(best_fit_fluxes(index));
IDH1_mean_flux_abs   = abs(mean_fluxes(index));
IDH1_gene_expression = enzyme_gene_expression(index);

gene_name='IDH2';
index           = find(ismember(enzymes,gene_name));
IDH2_best_fit_flux_abs   = abs(best_fit_fluxes(index));
IDH2_mean_flux_abs   = abs(mean_fluxes(index));
IDH2_gene_expression = enzyme_gene_expression(index);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot isoenzymes fold change of gene expression and flux
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot best fit fluxes
figure('units','normalized','outerposition',[0 0 1 1])
% fold change of MT/CY fluxes and gene expressions
GOT_MT_CY_fold_change_best_fit_flux_abs = GOT2_best_fit_flux_abs/GOT1_best_fit_flux_abs;
GOT_MT_CY_fold_change_gene_expression   = GOT2_gene_expression/GOT1_gene_expression;
MDH_MT_CY_fold_change_best_fit_flux_abs = MDH2_best_fit_flux_abs/MDH1_best_fit_flux_abs;
MDH_MT_CY_fold_change_gene_expression   = MDH2_gene_expression/MDH1_gene_expression;
ME_MT_CY_fold_change_best_fit_flux_abs  = ME3_best_fit_flux_abs/ME1_best_fit_flux_abs;
ME_MT_CY_fold_change_gene_expression    = ME3_gene_expression/ME1_gene_expression;
IDH_MT_CY_fold_change_best_fit_flux_abs = IDH2_best_fit_flux_abs/IDH1_best_fit_flux_abs;
IDH_MT_CY_fold_change_gene_expression   = IDH2_gene_expression/IDH1_gene_expression;

y = log10([GOT_MT_CY_fold_change_best_fit_flux_abs GOT_MT_CY_fold_change_gene_expression; MDH_MT_CY_fold_change_best_fit_flux_abs MDH_MT_CY_fold_change_gene_expression; ME_MT_CY_fold_change_best_fit_flux_abs ME_MT_CY_fold_change_gene_expression; IDH_MT_CY_fold_change_best_fit_flux_abs IDH_MT_CY_fold_change_gene_expression]);
bar(y);
xlim([0.5 4.5]); 
legend({'Flux','Gene expression'});
title('Flux is taken from optimization best score');
xlabel('');
ylabel({'Mitochondrial vs. cytosolic', 'fold change [log10]'});
xticks([1 2 3 4]);
xticklabels({'GOT', 'MDH', 'ME', 'IDH'});
set(gcf,'color','w');
set(gca, 'FontSize', 32);
box on


% plot mean fluxes take nfrom confidence intervals (higher bound+lower
% bound)/2
figure('units','normalized','outerposition',[0 0 1 1])
% fold change of MT/CY fluxes and gene expressions
GOT_MT_CY_fold_change_mean_flux_abs     = GOT2_mean_flux_abs/GOT1_mean_flux_abs;
GOT_MT_CY_fold_change_gene_expression   = GOT2_gene_expression/GOT1_gene_expression;
MDH_MT_CY_fold_change_mean_flux_abs     = MDH2_mean_flux_abs/MDH1_mean_flux_abs;
MDH_MT_CY_fold_change_gene_expression   = MDH2_gene_expression/MDH1_gene_expression;
ME_MT_CY_fold_change_mean_flux_abs      = ME3_mean_flux_abs/ME1_mean_flux_abs;
ME_MT_CY_fold_change_gene_expression    = ME3_gene_expression/ME1_gene_expression;
IDH_MT_CY_fold_change_mean_flux_abs     = IDH2_mean_flux_abs/IDH1_mean_flux_abs;
IDH_MT_CY_fold_change_gene_expression   = IDH2_gene_expression/IDH1_gene_expression;

y = log10([GOT_MT_CY_fold_change_mean_flux_abs GOT_MT_CY_fold_change_gene_expression; MDH_MT_CY_fold_change_mean_flux_abs MDH_MT_CY_fold_change_gene_expression; ME_MT_CY_fold_change_mean_flux_abs ME_MT_CY_fold_change_gene_expression; IDH_MT_CY_fold_change_mean_flux_abs IDH_MT_CY_fold_change_gene_expression]);
bar(y);
xlim([0.5 4.5]); 
legend({'Flux','Gene expression'});
title('Flux is taken from mean confidence intervals');
xlabel('');
ylabel({'Mitochondrial vs. cytosolic', 'fold change [log10]'});
xticks([1 2 3 4]);
xticklabels({'GOT', 'MDH', 'ME', 'IDH'});
set(gcf,'color','w');
set(gca, 'FontSize', 32);
box on



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% scatter plot of all genes
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% net flux of best fit
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
figure('units','normalized','outerposition',[0 0 1 1])
scatter(log10(enzyme_gene_expression),log10(best_fit_fluxes),100,'filled');
text(log10(enzyme_gene_expression)-0.03, log10(best_fit_fluxes)-0.1, enzymes);
[R1,P1]=corrcoef(log10(enzyme_gene_expression),log10(best_fit_fluxes));
[RHO1,PVAL1] = corr(enzyme_gene_expression',best_fit_fluxes','Type','Spearman');
title(sprintf('Best fit flux vs gene expresssion: RHO=%f, PVAL=%e',RHO1,PVAL1));
%title(sprintf('Best fit flux vs gene expresssion: R=%f, P=%e',R1(1,2),P1(1,2)));
xlabel('Gene expression [log10(RPKM)]');
ylabel('Flux [log10(mM/h)]');
set(gcf,'color','w');
set(gca, 'FontSize', 30);
grid on;
box on

% net flux of mean flux
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
figure('Position',[1 1 1100 950]);
% scatter(log10(enzyme_gene_expression),log10(mean_fluxes),230,'filled');
scatter(log10(enzyme_gene_expression),log10(mean_fluxes),230,'*','LineWidth',2);
% scatter(log10(enzyme_gene_expression),log10(mean_fluxes),230,'x','LineWidth',5);
text_position = [log10(enzyme_gene_expression)-0.03;log10(mean_fluxes)-0.09];
text_position(:,16) = [log10(enzyme_gene_expression(16))-0.03 log10(mean_fluxes(16))+0.1];
text_position(:,1) = [log10(enzyme_gene_expression(1))+0.03 log10(mean_fluxes(1))];
text_position(:,3) = [log10(enzyme_gene_expression(3))-0.13 log10(mean_fluxes(3))-0.09];
text_position(:,4) = [log10(enzyme_gene_expression(4))-0.05 log10(mean_fluxes(4))-0.09];
text_position(:,6) = [log10(enzyme_gene_expression(6))+0.03 log10(mean_fluxes(6))];
text_position(:,8) = [log10(enzyme_gene_expression(8))+0.03 log10(mean_fluxes(8))];
text_position(:,10) = [log10(enzyme_gene_expression(10))-0.07 log10(mean_fluxes(10))-0.09];
text_position(:,11) = [log10(enzyme_gene_expression(11))+0.01 log10(mean_fluxes(11))-0.07];
text_position(:,14) = [log10(enzyme_gene_expression(14))+0.03 log10(mean_fluxes(14))];
text_position(:,15) = [log10(enzyme_gene_expression(15))+0.03 log10(mean_fluxes(15))+0.04];
text_position(:,19) = [log10(enzyme_gene_expression(19))+0.03 log10(mean_fluxes(19))];

text(text_position(1,:), text_position(2,:), enzymes, 'FontSize', 17);
[R2,P2]=corrcoef(log10(enzyme_gene_expression),log10(mean_fluxes));
[RHO2,PVAL2] = corr(enzyme_gene_expression',mean_fluxes','Type','Spearman');
title(sprintf('Spearman Correlation: RHO=%.2f, PVAL=%.2e',RHO2,PVAL2));
% title(sprintf('Mean flux vs gene expresssion: R=%f, P=%e',R2(1,2),P2(1,2)));
xlabel('Gene expression [log10(RPKM)]');
ylabel('Mean flux [log10(mM/h)]');
set(gcf,'color','w');
set(gca, 'FontSize', 28);
grid on;
xlim([min(log10(enzyme_gene_expression))-0.1 max(log10(enzyme_gene_expression))+0.2]);
h=lsline;
h.LineWidth = 4;
h.Color=GREEN_COLOR;
box on

set(gcf,'PaperSize',[50 30]);
s = sprintf('./output_images/3f.pdf');
saveas(gcf, s);


% Traditional MFA - net flux of mean flux
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 

figure('Position',[1 1 1100 950]);
% scatter(log10(enzyme_gene_expression),log10(mean_fluxes),230,'filled');
scatter(log10(enzyme_gene_expression),log10(mean_fluxes_MFA),230,'*','LineWidth',2);
text_position = [log10(enzyme_gene_expression)-0.03;log10(mean_fluxes_MFA)-0.18];
text_position(:,1)=text_position(:,1)+[0;0.35];
text_position(:,2)=text_position(:,2)+[-0.15;0.35];
text_position(:,3)=text_position(:,3)+[-0.05;0];
text_position(:,5)=text_position(:,5)+[0.05;0.1];
text_position(:,7)=text_position(:,7)+[0.05;0.1];
text_position(:,8)=text_position(:,8)+[0.06;0.2];
text_position(:,10)=text_position(:,10)+[-0.08;0.35];
text_position(:,11)=text_position(:,11)+[0;0.35];
text_position(:,12)=text_position(:,12)+[0;0.35];
text_position(:,13)=text_position(:,13)+[0.05;0.1];
text_position(:,15)=text_position(:,15)+[-0.05;0.35];
text_position(:,16)=text_position(:,16)+[-0.05;0];
text_position(:,17)=text_position(:,17)+[-0.07;0];

text(text_position(1,:), text_position(2,:), enzymes, 'FontSize', 17);
[R2,P2]=corrcoef(log10(enzyme_gene_expression),log10(mean_fluxes_MFA));
[RHO2,PVAL2] = corr(enzyme_gene_expression',mean_fluxes_MFA','Type','Spearman');
title(sprintf('Spearman Correlation: RHO=%.2f, PVAL=%.2e',RHO2,PVAL2));
% title(sprintf('Mean flux vs gene expresssion: R=%f, P=%e',R2(1,2),P2(1,2)));
xlabel('Gene expression [log10(RPKM)]');
ylabel('Mean flux [log10(mM/h)]');
set(gcf,'color','w');
set(gca, 'FontSize', 28);
grid on;
xlim([min(log10(enzyme_gene_expression))-0.1 max(log10(enzyme_gene_expression))+0.2]);
ylim([min(log10(mean_fluxes_MFA))-0.4 max(log10(mean_fluxes_MFA))+0.2]);
h=lsline;
h.LineWidth = 4;
h.Color=GREEN_COLOR;
box on

set(gcf,'PaperSize',[50 30]);
s = sprintf('./output_images/3g.pdf');
saveas(gcf, s);


% Both CoDe-MFA and traditional MFA - net flux of mean flux
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
figure('Position',[1 1 1100 950]);
% scatter(log10(enzyme_gene_expression),log10(mean_fluxes),230,'filled');
hold on;
scatter(log10(enzyme_gene_expression),log10(mean_fluxes),230,'*','LineWidth',2);
% scatter(log10(enzyme_gene_expression),log10(mean_fluxes),230,'x','LineWidth',5);
% [R2,P2]=corrcoef(log10(enzyme_gene_expression),log10(mean_fluxes));
% [RHO2,PVAL2] = corr(enzyme_gene_expression',mean_fluxes','Type','Spearman');
% title(sprintf('Spearman Correlation: RHO=%.2f, PVAL=%.2e',RHO2,PVAL2));
% title(sprintf('Mean flux vs gene expresssion: R=%f, P=%e',R2(1,2),P2(1,2)));
xlabel('Gene expression [log10(RPKM)]');
ylabel('Mean flux [log10(mM/h)]');
set(gcf,'color','w');
set(gca, 'FontSize', 28);
grid on;
xlim([min(log10(enzyme_gene_expression))-0.1 max(log10(enzyme_gene_expression))+0.2]);
scatter(log10(enzyme_gene_expression),log10(mean_fluxes_MFA),230,'*','LineWidth',2);
h=lsline;
set(h(2),'color',BLUE_COLOR)
set(h(1),'LineWidth',4)
set(h(1),'color',RED_COLOR)
set(h(2),'LineWidth',4)
legend('CoDe-MFA','MFA')
box on


