% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% Metabolite cytosolic/mitochondrial ratio - inferred by CODE-MFA
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
function output_per_concentration_result()
load('../mat_files/model_thermodynamics.mat','model_thermodynamics');    
load('../mat_files/sensitiviy_analysis_cofactors_ratio.mat','sensitiviy_analysis_cofactors_ratio');

sensitiviy_analysis_cofactors_ratio_Hela = sensitiviy_analysis_cofactors_ratio;
model_thermodynamics_Hela = model_thermodynamics;


load('../mat_files/sensitiviy_analysis_cofactors_ratio.mat', 'sensitiviy_analysis_cofactors_ratio');
load('../mat_files/model_thermodynamics.mat','model_thermodynamics');    


addpath('../functions/emu') 
addpath('../functions/general') 
addpath('../') 
run ../load_constants;

low_cy_mt_ratio=[];
high_cy_mt_ratio=[];
metabolite_name=cell(0);
for(i=1:length(sensitiviy_analysis_cofactors_ratio))
    % cofactors
    if((contains(sensitiviy_analysis_cofactors_ratio{i}.co_factors_name,'NADP'))||...
      (contains(sensitiviy_analysis_cofactors_ratio{i}.co_factors_name,'NADPH'))||...
      (contains(sensitiviy_analysis_cofactors_ratio{i}.co_factors_name,'NAD'))||...
      (contains(sensitiviy_analysis_cofactors_ratio{i}.co_factors_name,'NADH'))||...
      (contains(sensitiviy_analysis_cofactors_ratio{i}.co_factors_name,'GDP'))||...
      (contains(sensitiviy_analysis_cofactors_ratio{i}.co_factors_name,'GTP'))||...
      (contains(sensitiviy_analysis_cofactors_ratio{i}.co_factors_name,'ADP'))||...
      (contains(sensitiviy_analysis_cofactors_ratio{i}.co_factors_name,'ATP')))
          continue;
    % metabolites that are not coffactors
    else
        metabolite_name{end+1}  = sensitiviy_analysis_cofactors_ratio{i}.co_factors_name;
        low_cy_mt_ratio(end+1)  = sensitiviy_analysis_cofactors_ratio{i}.low_cofactor_ratio;
        high_cy_mt_ratio(end+1) = sensitiviy_analysis_cofactors_ratio{i}.high_cofactor_ratio;
    end    
end


x=[];
hf1=figure('Position',[1 1 800 1000]);
hold on;
YLabel={''};

% sort from high to low cy/mt ratio
[sorted_values sorted_indices] = sort(high_cy_mt_ratio);
% add MT near each CY metabolite
for(i=1:length(sorted_indices))
    metabolite_cy_mt_ratio_name = metabolite_name(sorted_indices(i));
        
    Y_index=length(YLabel);
    plot([log10(low_cy_mt_ratio(sorted_indices(i)));log10(high_cy_mt_ratio(sorted_indices(i)))],[Y_index Y_index],'Color',BLUE_COLOR,'LineWidth',5)
     
    YLabel{end+1}=sprintf('%s',metabolite_cy_mt_ratio_name{1});
    YLabel{end} = strrep(YLabel{end},'_',' ');    
end

line([0 0], [0 length(YLabel)+1], 'color','red','LineStyle',':','LineWidth',1.4);
set(gca, 'yTick', [0:length(YLabel)+1]);
YLabel{end+1}='';
xlabel('Metabolite concentration ratio confidence intervals [log10]', 'FontSize',26);
set(gca, 'YTickLabel', YLabel);    
grid;
ylim([0 length(YLabel)-1]);
xticks([log10(min(low_cy_mt_ratio)):1:log10(max(high_cy_mt_ratio))]);
hold off;
set(gca, 'FontSize',14);
set(gcf,'color','w');
ax = gca;
ax.Box='on';
ax.BoxStyle = 'full';




% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% Plot concentrations, sorted by Hela cell results, to show metabolites in
% the same order
low_cy_mt_ratio_Hela=[];
high_cy_mt_ratio_Hela=[];
metabolite_name_Hela=cell(0);
for(i=1:length(sensitiviy_analysis_cofactors_ratio_Hela))
    % cofactors
    if((contains(sensitiviy_analysis_cofactors_ratio_Hela{i}.co_factors_name,'NADP'))||...
      (contains(sensitiviy_analysis_cofactors_ratio_Hela{i}.co_factors_name,'NADPH'))||...
      (contains(sensitiviy_analysis_cofactors_ratio_Hela{i}.co_factors_name,'NAD'))||...
      (contains(sensitiviy_analysis_cofactors_ratio_Hela{i}.co_factors_name,'NADH'))||...
      (contains(sensitiviy_analysis_cofactors_ratio_Hela{i}.co_factors_name,'GDP'))||...
      (contains(sensitiviy_analysis_cofactors_ratio_Hela{i}.co_factors_name,'GTP'))||...
      (contains(sensitiviy_analysis_cofactors_ratio_Hela{i}.co_factors_name,'ADP'))||...
      (contains(sensitiviy_analysis_cofactors_ratio_Hela{i}.co_factors_name,'ATP')))
          continue;
    % metabolites that are not coffactors
    else
        metabolite_name_Hela{end+1}  = sensitiviy_analysis_cofactors_ratio_Hela{i}.co_factors_name;
        low_cy_mt_ratio_Hela(end+1)  = sensitiviy_analysis_cofactors_ratio_Hela{i}.low_cofactor_ratio;
        high_cy_mt_ratio_Hela(end+1) = sensitiviy_analysis_cofactors_ratio_Hela{i}.high_cofactor_ratio;
    end    
end



hf1=figure('Position',[1 1 800 1000]);
hold on;
YLabel={''};


% sort from high to low cy/mt ratio per Hela cells
[sorted_values sorted_indices] = sort(high_cy_mt_ratio_Hela);
% add MT near each CY metabolite
for(i=1:length(sorted_indices))
    metabolite_cy_mt_ratio_name_Hela = metabolite_name_Hela(sorted_indices(i));
    metabolite_cy_mt_ratio_index = strfind(metabolite_name,metabolite_cy_mt_ratio_name_Hela);
    metabolite_cy_mt_ratio_index=find(~cellfun(@isempty,metabolite_cy_mt_ratio_index));
        
    Y_index=length(YLabel);
    plot([log10(low_cy_mt_ratio(metabolite_cy_mt_ratio_index));log10(high_cy_mt_ratio(metabolite_cy_mt_ratio_index))],[Y_index Y_index],'Color',BLUE_COLOR,'LineWidth',5)
     
    YLabel{end+1}=sprintf('%s',metabolite_cy_mt_ratio_name_Hela{1});
    YLabel{end} = strrep(YLabel{end},'_',' ');    
end

line([0 0], [0 length(YLabel)+1], 'color','red','LineStyle',':','LineWidth',1.4);
set(gca, 'yTick', [0:length(YLabel)+1]);
YLabel{end+1}='';
xlabel('Metabolite concentration ratio confidence intervals [log10]', 'FontSize',26);
set(gca, 'YTickLabel', YLabel);    
grid;
ylim([0 length(YLabel)-1]);
% xlim(log10([low_cy_mt_ratio(sorted_indices(1)) high_cy_mt_ratio(sorted_indices(end))]));
xticks([log10(min(low_cy_mt_ratio)):1:log10(max(high_cy_mt_ratio))]);
hold off;
set(gca, 'FontSize',14);
set(gcf,'color','w');
ax = gca;
ax.Box='on';
ax.BoxStyle = 'full';



end     


