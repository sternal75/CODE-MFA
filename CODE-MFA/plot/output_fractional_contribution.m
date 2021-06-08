% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% This script verify the fractional contribution of Gln and Glc, and the
% sum of the two, in order to verify that there are no other carbon sources
% for the labled metabolites in the model
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
clear all;

addpath('../functions/emu') 
addpath('../functions/general') 
addpath(genpath('./cobratoolbox/src'))

% input glutamine MIDs
file_name = 'input_glutamine';
run ../processIsotopicLabel/calc;
met_list_norm_glutamine = met_list_norm;
clear met_list_norm;
% input glucose MIDs
file_name = 'input_glucose';
run ../processIsotopicLabel/calc;
met_list_norm_glucose = met_list_norm;
clear met_list_norm;
close all;

metabolites_with_MID_Glc = cell(0); fractional_contribution_Glc=[];
metabolites_with_MID_Gln = cell(0); fractional_contribution_Gln=[];
metabolites_with_MID_all     = cell(0); labeling_pattern_Gln=[];
%Glc
for(i=1:length(met_list_norm_glucose))
    metabolites_with_MID_all{end+1}         = met_list_norm_glucose{i}.met_name;
    metabolites_with_MID_Glc{end+1}         = met_list_norm_glucose{i}.met_name;
    fractional_contribution_Glc(end+1)    = sum(met_list_norm_glucose{i}.data'.*[0:length(met_list_norm_glucose{i}.data)-1])/(length(met_list_norm_glucose{i}.data)-1);
end

%Gln
for(i=1:length(met_list_norm_glutamine))
    metabolites_with_MID_all{end+1}         = met_list_norm_glutamine{i}.met_name;
    metabolites_with_MID_Gln{end+1}         = met_list_norm_glutamine{i}.met_name;
    fractional_contribution_Gln(end+1)    = sum(met_list_norm_glutamine{i}.data'.*[0:length(met_list_norm_glutamine{i}.data)-1])/(length(met_list_norm_glutamine{i}.data)-1);
end

metabolites_with_MID_all=unique(metabolites_with_MID_all);

for(i=1:length(metabolites_with_MID_all))
    index_glc = find(ismember(metabolites_with_MID_Glc,metabolites_with_MID_all{i}));
    index_gln = find(ismember(metabolites_with_MID_Gln,metabolites_with_MID_all{i}));
    fprintf('FC %s: Total=%.2f\n',metabolites_with_MID_all{i}, fractional_contribution_Glc(index_glc)+fractional_contribution_Gln(index_gln));
end


