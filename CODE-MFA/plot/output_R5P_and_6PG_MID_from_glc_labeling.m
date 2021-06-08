% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% The steady state fractional M+5 labeling of 6-phosphogluconate (6PG) 
% and ribulose-5-phonehate (R5P) when feeding [U-13C]-glucose
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 

clear all;
close all;

file_name = 'input_glucose';
addpath('../processIsotopicLabel') ;

[met_list_norm_neg met_name_arr_neg sample_list_short] = ProcessMavenIsotopicLabel(sprintf('../processIsotopicLabel/data/%s.xls', file_name), 3);
[met_list_norm_pos met_name_arr_pos sample_list_short] = ProcessMavenIsotopicLabel(sprintf('../processIsotopicLabel/data/%s.xls', file_name), 3);    


if 1==0
    met_list_norm = met_list_norm_neg;
else
    % Merge metabolite lists
    neg_only = setdiff(met_name_arr_neg, met_name_arr_pos);
    pos_only = setdiff(met_name_arr_pos, met_name_arr_neg);
    neg_and_pos = intersect(met_name_arr_pos, met_name_arr_neg);

    for i=1:length(neg_and_pos)
        neg_index = strmatch(neg_and_pos{i}, met_name_arr_neg, 'exact');    
        pos_index = strmatch(neg_and_pos{i}, met_name_arr_pos, 'exact');

        if (length(neg_index)~=1) || (length(pos_index)~=1)
            fprintf('Error..\n');
        end
        if met_list_norm_neg{neg_index}.median_intensity > met_list_norm_pos{pos_index}.median_intensity
            neg_only{end+1} = neg_and_pos{i};
        else
            pos_only{end+1} = neg_and_pos{i};
        end
    end

    met_list_norm = cell(0);
    for i=1:length(neg_only)
        neg_index = strmatch(neg_only{i}, met_name_arr_neg, 'exact');    
        met_list_norm{end+1} = met_list_norm_neg{neg_index};
    end

    for i=1:length(pos_only)
        pos_index = strmatch(pos_only{i}, met_name_arr_pos, 'exact');    
        met_list_norm{end+1} = met_list_norm_pos{pos_index};
    end
end

% Fix mas isotopomer distribution vector, by removing natural abundance
% and impurity effects
for i=1:length(met_list_norm)
    v=met_list_norm{i}.data;
    prob_labeled_carbon_is_C12 = 0.01;   % the inpurity of C13
    prob_natural_C13 = 0.011;  %C13 natual abundance
    
    u = AdjustMat(v', prob_labeled_carbon_is_C12, prob_natural_C13);    
    met_list_norm{i}.data  = u';
end





% Find major mass isotopomers
for i=1:length(met_list_norm)
    v = max(met_list_norm{i}.data, [], 2);
    t = find(v >= 0.01);
    met_list_norm{i}.mass_isotopomer_list = t;
end


legend_arr = cell(0);
for x=1:50
    legend_arr{x} = sprintf('m+%d', x-1);
end




num_of_labeling_forms=7;
last_bar = zeros(2,num_of_labeling_forms);
last_std = zeros(2,num_of_labeling_forms);
values_for_legend = [];
for i=1:length(met_list_norm)
    mat_bar = met_list_norm{i}.data(met_list_norm{i}.mass_isotopomer_list, :)';
    mat_var = met_list_norm{i}.var(met_list_norm{i}.mass_isotopomer_list, :)';
    mat_std = sqrt(mat_var);
    
    met_name = strrep(met_list_norm{i}.met_name,'_',' ');
    if(strcmp(met_name,'6phosphogluconate'))
        last_bar(1,met_list_norm{i}.mass_isotopomer_list)=mat_bar(end,:);
        last_bar(1,met_list_norm{i}.mass_isotopomer_list)=last_bar(1,met_list_norm{i}.mass_isotopomer_list)/sum(last_bar(1,met_list_norm{i}.mass_isotopomer_list));
        last_std(1,met_list_norm{i}.mass_isotopomer_list)=mat_std(end,:);
        values_for_legend = [values_for_legend;met_list_norm{i}.mass_isotopomer_list];
    end    
    if(strcmp(met_name,'Ribose5P'))
        last_bar(2,met_list_norm{i}.mass_isotopomer_list)=mat_bar(end,:);
        last_bar(2,met_list_norm{i}.mass_isotopomer_list)=last_bar(2,met_list_norm{i}.mass_isotopomer_list)/sum(last_bar(2,met_list_norm{i}.mass_isotopomer_list));
        last_std(2,met_list_norm{i}.mass_isotopomer_list)=mat_std(end,:);
        values_for_legend = [values_for_legend;met_list_norm{i}.mass_isotopomer_list];
    end
end



bar_1 = last_bar;
bar_2 = last_bar;
bar_1( :, all(~bar_1,1) ) = [];

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
figure;
values_for_legend = unique(sort(values_for_legend));
bar_label(bar_1, bar_1);
% title(sprintf('%s',met_name), 'FontSize', 18);
legend(legend_arr(values_for_legend), 'Location', 'NorthEastOutside');
xlim([0.2 2.7]);
xticks([1 2]);
xticklabels({'6PG', 'R5P'});
ylabel('Fraction of labeling','fontsize',20);
set(gca, 'Ylim', [0 1]);
set(gca, 'FontSize', 17);

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
figure('Position',[1 1 750 400]);
bar(bar_2');
hold on

x = [1:num_of_labeling_forms];
w=0.15;
mean_value = bar_2(1,:);
low_value = 1.95*last_std(1,:);
low_value((mean_value-low_value) < 0) = mean_value((mean_value-low_value) < 0);
high_value = 1.95*last_std(1,:);
high_value((mean_value+high_value) > 1) = 1-mean_value((mean_value+high_value) > 1);
er = errorbar(x(mean_value~=0)-w, mean_value(mean_value~=0),low_value(mean_value~=0),high_value(mean_value~=0), 'linewidth',2);
er.Color = [0 0 0];                            
er.LineStyle = 'none'; 
 
mean_value = bar_2(2,:);
low_value = 1.95*last_std(2,:);
low_value((mean_value-low_value) < 0) = mean_value((mean_value-low_value) < 0);
high_value = 1.95*last_std(2,:);
high_value((mean_value+high_value) > 1) = 1-mean_value((mean_value+high_value) > 1);
er = errorbar(x(mean_value~=0)+w, mean_value(mean_value~=0),low_value(mean_value~=0),high_value(mean_value~=0), 'linewidth',2);
er.Color = [0 0 0];                            
er.LineStyle = 'none'; 
hold off
yticks(0:0.2:1);
xlim([0.5 7.5]);
xticks(1:num_of_labeling_forms);
xticklabels({'m+0', 'm+1', 'm+2', 'm+3', 'm+4', 'm+5', 'm+6'});
ylabel({'Labeling fraction','(from [U-13C]-Glucose)'},'fontsize',20);
xtickangle(30);
legend({' 6PG',' R5P'}, 'Location', 'NorthWest')
set(gca, 'FontSize', 22);



% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
figure('Position',[1 1 650 400]);
bar_3 = bar_2(:,6);
last_std_3 = last_std(:,6);
bar(bar_3');
hold on

x = [1:num_of_labeling_forms];
mean_value = bar_3(1,:);
low_value = 1.95*last_std_3(1);
low_value((mean_value-low_value) < 0) = mean_value((mean_value-low_value) < 0);
high_value = 1.95*last_std_3(1);
high_value((mean_value+high_value) > 1) = 1-mean_value((mean_value+high_value) > 1);
er = errorbar(1, mean_value(mean_value~=0),low_value(mean_value~=0),high_value(mean_value~=0), 'linewidth',2);
er.Color = [0 0 0];                            
er.LineStyle = 'none'; 
 
mean_value = bar_3(2,:);
low_value = 1.95*last_std_3(2);
low_value((mean_value-low_value) < 0) = mean_value((mean_value-low_value) < 0);
high_value = 1.95*last_std_3(2);
high_value((mean_value+high_value) > 1) = 1-mean_value((mean_value+high_value) > 1);
er = errorbar(2, mean_value(mean_value~=0),low_value(mean_value~=0),high_value(mean_value~=0), 'linewidth',2);
er.Color = [0 0 0];                            
er.LineStyle = 'none'; 
hold off
yticks(0:0.2:1);
xlim([0.2 2.8]);
xticks(1:2);
xticklabels({'6PG', 'R5P'});
ylabel({'m+5 labeling fraction','(from [U-13C]-Glucose)'},'fontsize',20);
set(gca, 'FontSize', 22);

s = sprintf('./output_images/1b.pdf');
saveas(gcf, s);
