% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% CODE-MFA derived cytosolic and mitochondrial metabolite concentrations 
% (blue and green bars) and measured cellular concentrations (asterisk)
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 

function output_per_concentration_result()
close all;
load('../../CodeMFA_Hela_FINAL/mat_files/sensitiviy_analysis_concentration.mat', 'sensitiviy_analysis_concentration');
load('../../CodeMFA_Hela_FINAL/mat_files/model_thermodynamics.mat','model_thermodynamics');    
sensitiviy_analysis_concentration_Hela = sensitiviy_analysis_concentration;
model_thermodynamics_Hela = model_thermodynamics;


load('../mat_files/sensitiviy_analysis_concentration.mat', 'sensitiviy_analysis_concentration');
load('../mat_files/model_thermodynamics.mat','model_thermodynamics');    
load('../mat_files/directionalities.mat','directionalities');    


addpath('../functions/emu') 
addpath('../functions/general') 
addpath('../') 
run ../load_constants;

%     find the index of the best score among all directionalities
best_score = min(directionalities.errors);
index_best_score = find(directionalities.errors==min(directionalities.errors));
best_score_predicted_concentrations     = directionalities.predicted_concentrations(:,index_best_score);    

for(i=1:length(model_thermodynamics.mets))
    % low/high concentration (log10)
    low_concentration(i)   = log10(exp(sensitiviy_analysis_concentration{i}.low_concentration));
    high_concentration(i)  = log10(exp(sensitiviy_analysis_concentration{i}.high_concentration));   
end 

x=[];
err=[]; 
hf1=figure('Position',[1 1 700 1000]);
hold on;
YLabel={''};
concentrations_confidence_intervals = [];
FVA_confidence_interval_sizes=[];
MFA_with_thermodynamics_confidence_interval_sizes=[];
% find error and mean for all concentrations
for(i=1:length(sensitiviy_analysis_concentration))   
    if(isempty(sensitiviy_analysis_concentration{i}))
        x(end+1)=nan;
        err(end+1)=nan;        
    else
        x(end+1)=(low_concentration(i)+high_concentration(i))/2;        
        err(end+1)=high_concentration(i)-(low_concentration(i)+high_concentration(i))/2;        
    end    
end
% sort by max concentration of CY metabolites
% metabolite_measured_high_concentration=[-0.625471606 -2.026872146 -2.301029996 -1.22184875 -0.854430702 -1.552841969 -1.522878745 -1.420216403 -1.013228266 -0.946921557 -0.795880017 -0.649751982 -0.389339837 -0.289036881 -0.283996656 -0.251811973 -0.15739076 0.243038049 0.598790507 0.794488047 1.164352856 nan 1.192009593];
all_indices_CY=strfind(model_thermodynamics.mets,'_CY');
all_indices_CY=find(~cellfun(@isempty,all_indices_CY));
[sorted_values sorted_indices] = sort(high_concentration(all_indices_CY));
% add MT near each CY metabolite
CONST_VAL_FOR_PLOT_MT_ABOVE_CY=0.1;
for(i=1:length(sorted_indices))
    % do not output confidence intervals of this meatbolite
    if(model_thermodynamics.skip_sensitivity_analysis_metabolite_indices(all_indices_CY(sorted_indices(i)))==1)
        continue;
    end        
    metabolite_CY = model_thermodynamics.mets{all_indices_CY(sorted_indices(i))};
    metabolite_CY_index = all_indices_CY(sorted_indices(i));
    metabolite_MT = strrep(metabolite_CY,'_CY','_MT');
    metabolite_MT_index = find(strcmp(model_thermodynamics.mets,metabolite_MT));
    
    metabolite_name = strrep(metabolite_CY,'_CY','');
    
    % do not show co-factors
    if((strcmp(metabolite_name,'NAD'))      ||...
       (strcmp(metabolite_name,'NADH'))     ||...
       (strcmp(metabolite_name,'NADP'))     ||...
       (strcmp(metabolite_name,'NADPH'))    ||...
       (strcmp(metabolite_name,'ATP'))      ||...
       (strcmp(metabolite_name,'ADP'))      ||...
       (strcmp(metabolite_name,'GTP'))      ||...
       (strcmp(metabolite_name,'GDP')))
       continue;
    end
    metabolite_WC_index=strcmp(model_thermodynamics.WC.met_name,metabolite_name);
    metabolite_WC_index=find(metabolite_WC_index==1);
    metabolite_measured_high_concentration_log10 = log10(model_thermodynamics.WC.Concentrations(metabolite_WC_index));

    
    Y_index=length(YLabel);
    plot([x(metabolite_CY_index)-err(metabolite_CY_index);x(metabolite_CY_index)+err(metabolite_CY_index)],[Y_index-CONST_VAL_FOR_PLOT_MT_ABOVE_CY Y_index-CONST_VAL_FOR_PLOT_MT_ABOVE_CY],'Color',BLUE_COLOR,'LineWidth',5)
    
    
    % plot mitochondrial only if should not be skipped
    plot([x(metabolite_MT_index)-err(metabolite_MT_index);x(metabolite_MT_index)+err(metabolite_MT_index)],[Y_index+CONST_VAL_FOR_PLOT_MT_ABOVE_CY Y_index+CONST_VAL_FOR_PLOT_MT_ABOVE_CY],'Color',GREEN_COLOR,'LineWidth',5)
    
    if(~isnan(metabolite_measured_high_concentration_log10))
        plot(metabolite_measured_high_concentration_log10,Y_index,'*','Color',[0 0 0],'MarkerSize',12);
    end
 
    YLabel{end+1}=sprintf('%s',model_thermodynamics.mets{metabolite_CY_index});
    YLabel{end} = strrep(YLabel{end},'_CY','');    
    YLabel{end} = strrep(YLabel{end},'_',' ');    
end

line([0 0], [0 length(YLabel)+1], 'color','red','LineStyle',':','LineWidth',1.4);
set(gca, 'yTick', [0:length(YLabel)+1]);
YLabel{end+1}='';
xlabel('Metabolite concentration confidence intervals [log10(mM)]', 'FontSize',26);
set(gca, 'YTickLabel', YLabel);    
% hleg=legend('our method','FVA+thermodynamics');
% chleg = get(hleg,'children');
% set(chleg(2),'color','g');set(chleg(3),'color','g');
% set(chleg(5),'color','b');set(chleg(6),'color','b');
% title('Metabolite measured concentrations & our method','FontSize',32);
% ylim([0 length(YLabel)+1])
% annotation(hf1, 'textbox', [0,0,0.3,0.07], 'String', ({num_of_consistent_reactions_str, num_of_smaller_confidence_intervals_in_our_method_str}), 'VerticalAlignment', 'top', 'FontSize',14);
grid;
ylim([0 length(YLabel)-1]);
xlim([-5 2]);
hold off;
legend('Cytosolic','Mitochondrial','Measured (WC)', 'FontSize',12, 'Orientation','horizontal', 'Location', 'north');
set(gca, 'FontSize',14);
set(gcf,'color','w');
ax = gca;
ax.Box='on';
ax.BoxStyle = 'full';




% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% Plot concentrations, sorted by Hela cell results, to show metabolites in
% the same order
for(i=1:length(model_thermodynamics_Hela.mets))
    % low/high concentration (log10)
    low_concentration_Hela(i)   = log10(exp(sensitiviy_analysis_concentration_Hela{i}.low_concentration));
    high_concentration_Hela(i)  = log10(exp(sensitiviy_analysis_concentration_Hela{i}.high_concentration));   
end 


x=[];
err=[]; 
hf1=figure('Position',[1 1 700 1000]);
hold on;
YLabel={''};
concentrations_confidence_intervals = [];
FVA_confidence_interval_sizes=[];
MFA_with_thermodynamics_confidence_interval_sizes=[];
% find error and mean for all concentrations
for(i=1:length(sensitiviy_analysis_concentration))   
    if(isempty(sensitiviy_analysis_concentration{i}))
        x(end+1)=nan;
        err(end+1)=nan;        
    else
        x(end+1)=(low_concentration(i)+high_concentration(i))/2;        
        err(end+1)=high_concentration(i)-(low_concentration(i)+high_concentration(i))/2;        
    end    
end
% sort by max concentration of CY metabolites
% metabolite_measured_high_concentration=[-0.625471606 -2.026872146 -2.301029996 -1.22184875 -0.854430702 -1.552841969 -1.522878745 -1.420216403 -1.013228266 -0.946921557 -0.795880017 -0.649751982 -0.389339837 -0.289036881 -0.283996656 -0.251811973 -0.15739076 0.243038049 0.598790507 0.794488047 1.164352856 nan 1.192009593];
all_indices_Hela_CY=strfind(model_thermodynamics_Hela.mets,'_CY');
all_indices_Hela_CY=find(~cellfun(@isempty,all_indices_Hela_CY));
[sorted_values sorted_indices] = sort(high_concentration_Hela(all_indices_Hela_CY));

% add MT near each CY metabolite
CONST_VAL_FOR_PLOT_MT_ABOVE_CY=0.1;
for(i=1:length(sorted_indices))
    metabolite_CY = model_thermodynamics_Hela.mets{all_indices_Hela_CY(sorted_indices(i))};
    metabolite_CY_index = strcmp(model_thermodynamics.mets,metabolite_CY);
    metabolite_CY_index=find(metabolite_CY_index)
    metabolite_MT = strrep(metabolite_CY,'_CY','_MT');
    metabolite_MT_index = strcmp(model_thermodynamics.mets,metabolite_MT);
    metabolite_MT_index=find(metabolite_MT_index);
    
    % do not output confidence intervals of this meatbolite
    if(model_thermodynamics.skip_sensitivity_analysis_metabolite_indices(metabolite_CY_index)==1)
        continue;
    end        
    
    metabolite_name = strrep(metabolite_CY,'_CY','');
    
    % do not show co-factors
    if((strcmp(metabolite_name,'NAD'))      ||...
       (strcmp(metabolite_name,'NADH'))     ||...
       (strcmp(metabolite_name,'NADP'))     ||...
       (strcmp(metabolite_name,'NADPH'))    ||...
       (strcmp(metabolite_name,'ATP'))      ||...
       (strcmp(metabolite_name,'ADP'))      ||...
       (strcmp(metabolite_name,'GTP'))      ||...
       (strcmp(metabolite_name,'GDP')))
       continue;
    end
    metabolite_WC_index=strcmp(model_thermodynamics.WC.met_name,metabolite_name);
    metabolite_WC_index=find(metabolite_WC_index==1);
    metabolite_measured_high_concentration_log10 = log10(model_thermodynamics.WC.Concentrations(metabolite_WC_index));

    
    Y_index=length(YLabel);
    plot([x(metabolite_CY_index)-err(metabolite_CY_index);x(metabolite_CY_index)+err(metabolite_CY_index)],[Y_index-CONST_VAL_FOR_PLOT_MT_ABOVE_CY Y_index-CONST_VAL_FOR_PLOT_MT_ABOVE_CY],'Color',BLUE_COLOR,'LineWidth',5)
    
    
    % plot mitochondrial only if should not be skipped
    plot([x(metabolite_MT_index)-err(metabolite_MT_index);x(metabolite_MT_index)+err(metabolite_MT_index)],[Y_index+CONST_VAL_FOR_PLOT_MT_ABOVE_CY Y_index+CONST_VAL_FOR_PLOT_MT_ABOVE_CY],'Color',GREEN_COLOR,'LineWidth',5)
    
    if(~isnan(metabolite_measured_high_concentration_log10))
        plot(metabolite_measured_high_concentration_log10,Y_index,'*','Color',[0 0 0],'MarkerSize',12);
    end
 
    YLabel{end+1}=sprintf('%s',model_thermodynamics.mets{metabolite_CY_index});
    YLabel{end} = strrep(YLabel{end},'_CY','');    
    YLabel{end} = strrep(YLabel{end},'_',' ');    
end

line([0 0], [0 length(YLabel)+1], 'color','red','LineStyle',':','LineWidth',1.4);
set(gca, 'yTick', [0:length(YLabel)+1]);
YLabel{end+1}='';
xlabel('Metabolite concentration confidence intervals [mM(log10)]', 'FontSize',26);
set(gca, 'YTickLabel', YLabel);    
grid;
ylim([0 length(YLabel)-1]);
xlim([-5 2]);
hold off;
legend('Cytosolic','Mitochondrial','Measured (WC)', 'FontSize',12, 'Orientation','horizontal', 'Location', 'north');
set(gca, 'FontSize',14);
set(gcf,'color','w');
ax = gca;
ax.Box='on';
ax.BoxStyle = 'full';


set(gcf,'PaperSize',[50 30]);
s = sprintf('./output_images/4d.pdf');
saveas(gcf, s);



% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% Plot measured concentrations confidence intervals
hf1=figure('Position',[1 1 700 1000]);
measured_lower_bound=[];
measured_upper_bound=[];
YLabel=cell(0);
metabolite_exists_in_one_compartment = {'6phosphogluconate_CY','Glc6P_CY','Ribose5P_CY','G3P','13BPG','3PG'};
for(i=1:length(metabolite_exists_in_one_compartment))
    index_compartmentalized_con = find(ismember(model_thermodynamics.mets,metabolite_exists_in_one_compartment{i}));
    measured_lower_bound(end+1,1) = exp(model_thermodynamics.mets_lb(index_compartmentalized_con));
    measured_upper_bound(end+1,1) = exp(model_thermodynamics.mets_ub(index_compartmentalized_con));
    YLabel{end+1} = sprintf('%s',metabolite_exists_in_one_compartment{i});
    YLabel{end} = strrep(YLabel{end},'_CY','');        
end
measured_lower_bound = [measured_lower_bound;model_thermodynamics.WC.Concentrations-2*model_thermodynamics.WC.Concentrations_STD];
measured_lower_bound(measured_lower_bound<1e-5)=1e-5;
measured_upper_bound = [measured_upper_bound;model_thermodynamics.WC.Concentrations+2*model_thermodynamics.WC.Concentrations_STD];
number_of_measured_mets = length(measured_lower_bound);

for(i=1:length(model_thermodynamics.WC.met_name))    
    YLabel{end+1} = sprintf('%s',model_thermodynamics.WC.met_name{i});
    YLabel{end} = strrep(YLabel{end},'_',' ');    
end
[sorted_values sorted_indices] = sort(measured_upper_bound);
measured_lower_bound_sorted = measured_lower_bound(sorted_indices);
measured_upper_bound_sorted = measured_upper_bound(sorted_indices);
YLabel_sorted = YLabel(sorted_indices);
plot(log10([measured_lower_bound_sorted';measured_upper_bound_sorted']),[1:number_of_measured_mets;1:number_of_measured_mets],'Color',BLUE_COLOR,'LineWidth',5);
ylim([0 number_of_measured_mets+1]);
yticks([1:number_of_measured_mets]);
set(gca, 'YTickLabel', YLabel_sorted);    
line([0 0], [0 number_of_measured_mets+1], 'color','red','LineStyle',':','LineWidth',1.4);
xlabel('Metabolite measured concentration confidence intervals [mM(log10)]', 'FontSize',26);
grid;
set(gca, 'FontSize',14);
set(gcf,'color','w');
ax = gca;
ax.Box='on';
ax.BoxStyle = 'full';

% Print all measured metaboilte concentrations
fprintf('******* Measured concentrations lower and upper bounds *******\n');
for(i=1:length(YLabel_sorted))
    fprintf('Measured concentration for %s (mM): [%f %f]\n',YLabel_sorted{i}, measured_lower_bound_sorted(i), measured_upper_bound_sorted(i))
end
fprintf('******************************************************************************\n');


output_for_excel = cell(0);
CODE_MFA_INFERRED_concentration_lb=[];
CODE_MFA_INFERRED_concentration_ub=[];
for(i=1:length(sensitiviy_analysis_concentration))
    output_for_excel{end+1} = sensitiviy_analysis_concentration{i}.metabolite_name;
    CODE_MFA_INFERRED_concentration_lb(end+1) = sensitiviy_analysis_concentration{i}.low_concentration;
    CODE_MFA_INFERRED_concentration_ub(end+1) = sensitiviy_analysis_concentration{i}.high_concentration;    
end

% print gln labeling to excel
xlswrite('temp.xlsx',{'met name'},'temp','A1');  
xlswrite('temp.xlsx',{'CODE-MFA inferred lb'},'temp','b1');  
xlswrite('temp.xlsx',{'CODE-MFA inferred ub'},'temp','c1');  
xlswrite('temp.xlsx',{'CODE-MFA inferred best fit'},'temp','d1');  
xlswrite('temp.xlsx',output_for_excel','temp','A2');  
xlswrite('temp.xlsx',CODE_MFA_INFERRED_concentration_lb','temp','b2');  
xlswrite('temp.xlsx',CODE_MFA_INFERRED_concentration_ub','temp','c2');  
xlswrite('temp.xlsx',best_score_predicted_concentrations,'temp','d2');  







end     


