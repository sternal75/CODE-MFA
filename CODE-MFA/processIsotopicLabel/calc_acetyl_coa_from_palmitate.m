close all;
clear all;

STD_OF_MEASUREMNTS = 0.05;
CHI_SQUARE_DISTRIBUTION_ONE_DEGREE_OF_FREEDON = 3.84;

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% Glutamine labeling
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% isnert labeling of palmitic acid from gln [m+0 m+2 m+4 m+6 m+8 m+10 m+12 m+14 m+16]
palmitate_gln=[0    0.2097    0.2351    0.1952    0.1486    0.1128    0.0641    0.0345 0];

% remove m+0
palmitate_gln(1)=0;

% normalize all to the sum of the labeling pattern, so it will sum-up to 1
% (after removing m+0)
palmitate_gln=palmitate_gln/sum(palmitate_gln);

pow_distance=[];
p=0:0.001:1;
simulated_distribution=[];
for(i=1:8)
    a=binopdf(i,8,p);
    simulated_distribution=[simulated_distribution;a];
    pow_distance = [pow_distance;(a-palmitate_gln(i+1)).^2];
end
[min_val min_pos]=min(sum(pow_distance));

% to find confidence intervals, find the score (sum of square residuals of
% the std), and look for the values comply with chi square distribution
% with one degree of freedom
sum_of_square_residuals = sum(pow_distance/(STD_OF_MEASUREMNTS^2));
all_positions_within_confidence_intervals = sum_of_square_residuals<(min(sum_of_square_residuals)+3.84);
position_of_lower_bound = min(find(all_positions_within_confidence_intervals==1));
position_of_upper_bound = max(find(all_positions_within_confidence_intervals==1));
figure;
hold on;
plot(palmitate_gln(2:end),'*','MarkerSize',22, 'LineWidth', 3);
plot(simulated_distribution(:,min_pos),'x','MarkerSize',22, 'LineWidth', 3);
legend('measured','simulated')
box on;
set(gcf,'color','w');
set(gca, 'FontSize', 30);
xlim([1 8]);
title(sprintf('acetyl coa m+2 from glutamine labeling=%f\n',p(min_pos)));
grid on;
xticks([1:8])
xticklabels({'m+2', 'm+4', 'm+6', 'm+8', 'm+10', 'm+12', 'm+14', 'm+16'});
hold off;
fprintf('acetyl coa m+2 from glutamine labeling=%f, lower bound=%f, upper bound=%f\n', p(min_pos), p(position_of_lower_bound), p(position_of_upper_bound));


% % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% Glucose labeling
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % 

% isnert labeling of palmitic acid from glc [m+0 m+2 m+4 m+6 m+8 m+10 m+12 m+14 m+16]
palmitate_glc=[0    0.0871    0.0842    0.1003    0.1413    0.1960    0.2013    0.1403    0.0494];

% remove m+0
palmitate_glc(1)=0;

% normalize all to the sum of the labeling pattern, so it will sum-up to 1
% (after removing m+0)
palmitate_glc=palmitate_glc/sum(palmitate_glc);

pow_distance=[];
p=0:0.001:1;
simulated_distribution=[];
for(i=1:8)
    a=binopdf(i,8,p);
    simulated_distribution=[simulated_distribution;a];
    pow_distance = [pow_distance;(a-palmitate_glc(i+1)).^2];
end
[min_val min_pos]=min(sum(pow_distance));

sum_of_square_residuals = sum(pow_distance/(STD_OF_MEASUREMNTS^2));
all_positions_within_confidence_intervals = sum_of_square_residuals<(min(sum_of_square_residuals)+3.84);
position_of_lower_bound = min(find(all_positions_within_confidence_intervals==1));
position_of_upper_bound = max(find(all_positions_within_confidence_intervals==1));

figure;
hold on;
plot(palmitate_glc(2:end),'*','MarkerSize',22, 'LineWidth', 3);
plot(simulated_distribution(:,min_pos),'x','MarkerSize',22, 'LineWidth', 3);
legend('measured','simulated')
box on;
set(gcf,'color','w');
set(gca, 'FontSize', 30);
xlim([1 8]);
title(sprintf('acetyl coa m+2 from glucose labeling=%f\n',p(min_pos)));
grid on;
xticks([1:8])
xticklabels({'m+2', 'm+4', 'm+6', 'm+8', 'm+10', 'm+12', 'm+14', 'm+16'});
hold off;

fprintf('acetyl coa m+2 from glucose labeling=%f, lower bound=%f, upper bound=%f\n',p(min_pos), p(position_of_lower_bound), p(position_of_upper_bound));
