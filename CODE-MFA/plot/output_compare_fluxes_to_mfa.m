% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% Cumulative distribution of reaction net flux confidence interval 
% sizes inferred by CODE-MFA versus MFA and CODE-MFA without 
% thermodynamic considerations 
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
clear all;
close all;


load('../mat_files/net_fluxes.mat', 'net_fluxes');
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
run ../load_constants;

index_best_score = find(directionalities.errors==min(directionalities.errors));    
best_score_predicted_net_fluxes = directionalities.predicted_net_fluxes(:,index_best_score);


for(i=1:length(net_fluxes))
    low_flux_our_method(i) = net_fluxes{i}.low_flux;
    high_flux_our_method(i) = net_fluxes{i}.high_flux;
end 




% bidirectional fluxes
bidirectional_fluxes=(strfind(model_net_fluxes.rxns,'f'))
bidirectional_fluxes=find(~cellfun(@isempty,bidirectional_fluxes));

[dFluxes,tFluxes] = xlsread('../xls_input_files/mfa_and_code_mfa_without_thermodynamics_results.xlsx');


% compare to regular MFA
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
low_flux_MFA    = dFluxes(:,1)';
high_flux_MFA   = dFluxes(:,2)';


% comulative distribution - bidirectional fluxes
MAX_confidence_interval_for_axes_1000 = 1000;
MFA_confidence_intervals_size = high_flux_MFA-low_flux_MFA;
MFA_confidence_intervals_size(MFA_confidence_intervals_size>MAX_confidence_interval_for_axes_1000)=MAX_confidence_interval_for_axes_1000;
MFA_confidence_intervals_size = MFA_confidence_intervals_size(bidirectional_fluxes);
our_method_confidence_intervals_size = high_flux_our_method-low_flux_our_method;
our_method_confidence_intervals_size(our_method_confidence_intervals_size>MAX_confidence_interval_for_axes_1000)=MAX_confidence_interval_for_axes_1000;
our_method_confidence_intervals_size = our_method_confidence_intervals_size(bidirectional_fluxes);
x_our_method = unique(sort(log10(our_method_confidence_intervals_size)));
x_mfa        = unique(sort(log10(MFA_confidence_intervals_size)));
x_our_method = unique(sort(log10(our_method_confidence_intervals_size)));
x_mfa        = unique(sort(log10(MFA_confidence_intervals_size)));
our_method_num_of_sensitivity_analysis_range_smaller_than_x=[];
MFA_num_of_sensitivity_analysis_range_smaller_than_x=[];
for(i=1:length(x_our_method))
    our_method_num_of_sensitivity_analysis_range_smaller_than_x(i)  =  sum(log10(our_method_confidence_intervals_size) <= x_our_method(i));
end
for(i=1:length(x_mfa))
    MFA_num_of_sensitivity_analysis_range_smaller_than_x(i)         =  sum(log10(MFA_confidence_intervals_size) <= x_mfa(i));
end
figure('Position',[1 1 1100 950]);
hold on;
% plot(x_our_method, our_method_num_of_sensitivity_analysis_range_smaller_than_x/length(our_method_confidence_intervals_size),'-x','MarkerSize',14,'LineWidth',2);
% plot(x_mfa, MFA_num_of_sensitivity_analysis_range_smaller_than_x/length(MFA_confidence_intervals_size),'-x','MarkerSize',14,'LineWidth',2);
plot(x_our_method, our_method_num_of_sensitivity_analysis_range_smaller_than_x/length(our_method_confidence_intervals_size),'-*','MarkerSize',14,'LineWidth',2);
plot(x_mfa, MFA_num_of_sensitivity_analysis_range_smaller_than_x/length(MFA_confidence_intervals_size),'-*','MarkerSize',14,'LineWidth',2);
% plot(x_our_method, our_method_num_of_sensitivity_analysis_range_smaller_than_x/length(our_method_confidence_intervals_size),'.-','MarkerSize',26,'LineWidth',2);
% plot(x_mfa, MFA_num_of_sensitivity_analysis_range_smaller_than_x/length(MFA_confidence_intervals_size),'.-','MarkerSize',26,'LineWidth',2);
hold off;
xlabel('Net flux confidence interval size [log10(mM/h)]');
ylabel('Fraction of reactions');
% title('Comulative distribution of net flux confidence intervals');
legend('CoDe-MFA','MFA' );
grid on
set(gcf,'color','w');
set(gca, 'FontSize', 28);
ax=gca;
ax.Box='on';
ylim([0 1.02]);
xlim([0 str2num(ax.XTickLabel{end})]);
ax.XTickLabel{end}=['>' ax.XTickLabel{end}]




% print bars with number of fluxes with known directionality
figure('Position',[1 1 550 600]);
unknown_directionalities = sign(low_flux_MFA.*high_flux_MFA);
unknown_directionalities(unknown_directionalities==1)=0;
unknown_directionalities(unknown_directionalities==-1)=1;
num_of_unknown_directionaloities_mfa = sum(unknown_directionalities);
unknown_directionalities = sign(low_flux_our_method.*high_flux_our_method);
unknown_directionalities(unknown_directionalities==1)=0;
unknown_directionalities(unknown_directionalities==-1)=1;
num_of_unknown_directionaloities_our_method = sum(unknown_directionalities);
bar(1,100-num_of_unknown_directionaloities_our_method*100/sum(model_net_fluxes.is_net_flux));
hold on
bar(2,100-num_of_unknown_directionaloities_mfa*100/sum(model_net_fluxes.is_net_flux));
xticks([1 2]);
xticklabels({'CoDe-MFA','MFA'});
set(gcf,'color','w');
set(gca, 'FontSize',20);
ylabel('% of reactions with known directionality');
ylim([0 100]);
a=[cellstr(num2str(get(gca,'ytick')'))]; 
pct = char(ones(size(a,1),1)*'%'); 
new_yticks = [char(a),pct];
set(gca,'yticklabel',new_yticks)
grid on
hold off




% compare to CoDe-MFA without thermodynamics
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
low_flux_code_mfa_without_thermodynamics    = dFluxes(:,4)';
high_flux_code_mfa_without_thermodynamics   = dFluxes(:,5)';



% comulative distribution - bidirectional fluxes
MAX_confidence_interval_for_axes_1000 = 1000;
code_mfa_without_thermodynamics_confidence_intervals_size = high_flux_code_mfa_without_thermodynamics-low_flux_code_mfa_without_thermodynamics;
code_mfa_without_thermodynamics_confidence_intervals_size(code_mfa_without_thermodynamics_confidence_intervals_size>MAX_confidence_interval_for_axes_1000)=MAX_confidence_interval_for_axes_1000;
code_mfa_without_thermodynamics_confidence_intervals_size = code_mfa_without_thermodynamics_confidence_intervals_size(bidirectional_fluxes);
our_method_confidence_intervals_size = high_flux_our_method-low_flux_our_method;
our_method_confidence_intervals_size(our_method_confidence_intervals_size>MAX_confidence_interval_for_axes_1000)=MAX_confidence_interval_for_axes_1000;
our_method_confidence_intervals_size = our_method_confidence_intervals_size(bidirectional_fluxes);
x_our_method = unique(sort(log10(our_method_confidence_intervals_size)));
x_code_mfa_without_thermodynamics        = unique(sort(log10(code_mfa_without_thermodynamics_confidence_intervals_size)));
our_method_num_of_sensitivity_analysis_range_smaller_than_x=[];
code_mfa_without_thermodynamics_num_of_sensitivity_analysis_range_smaller_than_x=[];
for(i=1:length(x_our_method))
    our_method_num_of_sensitivity_analysis_range_smaller_than_x(i)  =  sum(log10(our_method_confidence_intervals_size) <= x_our_method(i));
end
for(i=1:length(x_code_mfa_without_thermodynamics))
    code_mfa_without_thermodynamics_num_of_sensitivity_analysis_range_smaller_than_x(i)         =  sum(log10(code_mfa_without_thermodynamics_confidence_intervals_size) <= x_code_mfa_without_thermodynamics(i));
end

figure('Position',[1 1 1100 950]);
hold on;
plot(x_our_method, our_method_num_of_sensitivity_analysis_range_smaller_than_x/length(our_method_confidence_intervals_size),'-*','MarkerSize',14,'LineWidth',2);
plot(x_code_mfa_without_thermodynamics, code_mfa_without_thermodynamics_num_of_sensitivity_analysis_range_smaller_than_x/length(code_mfa_without_thermodynamics_confidence_intervals_size),'-*','MarkerSize',14,'LineWidth',2);
hold off;

xlabel('Net flux confidence interval size [log10(mM/h)]');
ylabel('Fraction of reactions');
% title('Comulative distribution of net flux confidence intervals');
legend('CoDe-MFA','CoDe-MFA without thermodynamics' );
grid on
set(gcf,'color','w');
set(gca, 'FontSize',18);
ax=gca;
ax.Box='on';
ylim([0 1.02]);
xlim([0 str2num(ax.XTickLabel{end})]);
% ax.XTickLabel{end}=['>' ax.XTickLabel{end}]



% print bars with number of fluxes with known directionality
figure('Position',[1 1 550 600]);
unknown_directionalities = sign(low_flux_code_mfa_without_thermodynamics.*high_flux_code_mfa_without_thermodynamics);
unknown_directionalities(unknown_directionalities==1)=0;
unknown_directionalities(unknown_directionalities==-1)=1;
num_of_unknown_directionaloities_code_mfa_without_thermodynamics = sum(unknown_directionalities);
unknown_directionalities = sign(low_flux_our_method.*high_flux_our_method);
unknown_directionalities(unknown_directionalities==1)=0;
unknown_directionalities(unknown_directionalities==-1)=1;
num_of_unknown_directionaloities_our_method = sum(unknown_directionalities);
bar(1,100-num_of_unknown_directionaloities_our_method*100/sum(model_net_fluxes.is_net_flux));
hold on
bar(2,100-num_of_unknown_directionaloities_code_mfa_without_thermodynamics*100/sum(model_net_fluxes.is_net_flux));
xticks([1 2]);
xticklabels({'CoDe-MFA','CoDe-MFA without thermodynamics'});
set(gcf,'color','w');
set(gca, 'FontSize',20);
ylabel('% of reactions with known directionality');
ylim([0 100]);
a=[cellstr(num2str(get(gca,'ytick')'))]; 
pct = char(ones(size(a,1),1)*'%'); 
new_yticks = [char(a),pct];
set(gca,'yticklabel',new_yticks)
grid on
hold off



% compare to regular MFA and to CoDe-MFA without thermodynamics
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 

% comulative distribution - bidirectional fluxes
figure('Position',[1 1 1100 950]);
hold on;
plot(x_our_method, our_method_num_of_sensitivity_analysis_range_smaller_than_x/length(our_method_confidence_intervals_size),'-*','MarkerSize',14,'LineWidth',2);
plot(x_code_mfa_without_thermodynamics, code_mfa_without_thermodynamics_num_of_sensitivity_analysis_range_smaller_than_x/length(code_mfa_without_thermodynamics_confidence_intervals_size),'-*','MarkerSize',14,'LineWidth',2);
plot(x_mfa, MFA_num_of_sensitivity_analysis_range_smaller_than_x/length(MFA_confidence_intervals_size),'-*','MarkerSize',14,'LineWidth',2);
hold off;

xlabel('Net flux confidence interval size [log10(mM/h)]');
ylabel('Fraction of reactions');
% title('Comulative distribution of net flux confidence intervals');
legend('CoDe-MFA','CoDe-MFA w/o thermo','MFA', 'Location','northwest');
grid on
set(gcf,'color','w');
set(gca, 'FontSize', 28);
ax=gca;
ax.Box='on';
ylim([0 1.02]);
xlim([0 str2num(ax.XTickLabel{end})]);
ax.XTickLabel{end}=['>' ax.XTickLabel{end}]

set(gcf,'PaperSize',[50 30]);
s = sprintf('./output_images/3d.pdf');
saveas(gcf, s);


% print bars with number of fluxes with known directionality
% figure('Position',[1 1 550 600]);
figure('Position',[1 1 1100 950]);
unknown_directionalities = sign(low_flux_code_mfa_without_thermodynamics.*high_flux_code_mfa_without_thermodynamics);
unknown_directionalities(unknown_directionalities==1)=0;
unknown_directionalities(unknown_directionalities==-1)=1;
num_of_unknown_directionaloities_code_mfa_without_thermodynamics = sum(unknown_directionalities);
unknown_directionalities = sign(low_flux_our_method.*high_flux_our_method);
unknown_directionalities(unknown_directionalities==1)=0;
unknown_directionalities(unknown_directionalities==-1)=1;
num_of_unknown_directionaloities_our_method = sum(unknown_directionalities);
bar(1,100-num_of_unknown_directionaloities_our_method*100/sum(model_net_fluxes.is_net_flux));
hold on
bar(2,100-num_of_unknown_directionaloities_code_mfa_without_thermodynamics*100/sum(model_net_fluxes.is_net_flux));
bar(3,100-num_of_unknown_directionaloities_mfa*100/sum(model_net_fluxes.is_net_flux));
xticks([1 2 3]);
row1 = {'CoDe-MFA' 'CoDe-MFA' 'MFA'};
row2 = {'' 'w/o thermo' ''};
labelArray = [row1; row2]; 
tickLabels = strtrim(sprintf('%s\\newline%s\n', labelArray{:}));
xticklabels(tickLabels);
% xticklabels({'CoDe-MFA','CoDe-MFA w/o thermo', 'MFA'});
xtickangle(45)
set(gcf,'color','w');
set(gca, 'FontSize', 28);
ylabel('% of reactions with known directionality');
ylim([0 100]);
xlim([0.5 3.5]);
a=[cellstr(num2str(get(gca,'ytick')'))]; 
pct = char(ones(size(a,1),1)*'%'); 
new_yticks = [char(a),pct];
set(gca,'yticklabel',new_yticks)
grid on
xlim([0.2 3.8]);
hold off

set(gcf,'PaperSize',[50 30]);
s = sprintf('./output_images/3e.pdf');
saveas(gcf, s);



% calculate Wilcoxon rank sum test
p = ranksum(MFA_confidence_intervals_size,our_method_confidence_intervals_size)
fprintf('Wilcoxon rank sum test for MFA compared to CoDe-MFA. P-value=%e\n', p);
p = ranksum(code_mfa_without_thermodynamics_confidence_intervals_size,our_method_confidence_intervals_size)
fprintf('Wilcoxon rank sum test for CoDe-MFA w/o thermodynamics compared to CoDe-MFA. P-value=%e\n', p);


xlswrite('temp.xlsx',{'reaction'},'temp','A1');  
xlswrite('temp.xlsx',{'CoDe-MFA low flux'},'temp','b1');  
xlswrite('temp.xlsx',{'CoDe-MFA high flux'},'temp','c1');  
xlswrite('temp.xlsx',{'CODE MFA - best fit'},'temp','d1');  
xlswrite('temp.xlsx',model_net_fluxes.rxns','temp','A2');  
xlswrite('temp.xlsx',low_flux_our_method','temp','b2');  
xlswrite('temp.xlsx',high_flux_our_method','temp','c2');  
xlswrite('temp.xlsx',best_score_predicted_net_fluxes,'temp','d2');  


