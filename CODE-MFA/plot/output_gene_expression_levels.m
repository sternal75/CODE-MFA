% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% Gene expression levels obtained from the Broad Institute 
% Cancer Cell Line Encyclopedia (CCLE), with 0.5 RPKM value 
% as the cutoff to determine non-expressed genes
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 

close all;
clear all;
run ../load_constants;

% get all genes in our model, with their gene expression levels (RPKM)
all_genes = readtable('../xls_input_files/gene_expression_all_cell_lines.xlsx', 'ReadVariableNames', true);
figure('Position',[100 100 1500 500]);
x=1:size(all_genes,1);

[sorted_rpkm sorted_index] = sort(all_genes.GeneExpressionHCT_RPKM_);
% Set minimum RPKM to 0.1, as we want to plot the RPKM values in log10
% scale
all_genes.GeneExpressionHela_RPKM_(all_genes.GeneExpressionHela_RPKM_<0.1)=0.1;
all_genes.GeneExpressionHCT_RPKM_(all_genes.GeneExpressionHCT_RPKM_<0.1)=0.1;
all_genes.GeneExpressionA549_RPKM_(all_genes.GeneExpressionA549_RPKM_<0.1)=0.1;
all_genes.GeneExpressionLN229_RPKM_(all_genes.GeneExpressionLN229_RPKM_<0.1)=0.1;

bar(x,log10([all_genes.GeneExpressionHela_RPKM_(sorted_index) all_genes.GeneExpressionHCT_RPKM_(sorted_index) all_genes.GeneExpressionA549_RPKM_(sorted_index) all_genes.GeneExpressionLN229_RPKM_(sorted_index)]));
xticks([1:1:length(x)]);
xticklabels(all_genes.Enzyme(sorted_index));
xtickangle(90);
xlim([0 length(x)+1]);
set(gcf,'color','w');
set(gca, 'FontSize',14);
xlabel('Genes');
ylabel('Gene expression level [log10(RPKM)]');
grid on;
line([0 length(x)+1], log10([0.5 0.5]), 'color','red','LineStyle','--','LineWidth',2.5);
text(4,-0.45,'RPKM cutoff level','color',RED_COLOR, 'FontSize', 16);
legend('Hela cells','HCT116 cells  .','A549 cells','LN229 cells', 'Location', 'northwest');
ax=gca;
ax.YTickLabel{1}=['<' ax.YTickLabel{1}]

set(gcf,'PaperSize',[50 30]);
s = sprintf('./output_images/s1.pdf');
saveas(gcf, s);
