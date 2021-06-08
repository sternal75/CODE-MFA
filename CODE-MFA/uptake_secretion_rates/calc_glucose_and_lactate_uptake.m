clear all;
file_name='data_for_uptake_secretion_for_glc_and_lactate';
addpath('../processIsotopicLabel') 
[met_list_norm met_name_arr sample_list_short] = ProcessMavenIsotopicLabel(sprintf('%s.xlsx', file_name), 3);
[d,t] = xlsread(sprintf('%s.xlsx', file_name),'PCV')
[d1,t1] = xlsread(sprintf('%s.xlsx', file_name),'metabolite_preliminary_info');
[d2,t2] = xlsread(sprintf('%s.xlsx', file_name),'experimental_constants')
exist_in_bottle_media=d1(:,1);
fully_labled_in_bottle_media=d1(:,2);

% experimental constants
% experimental constants
MEDIA_FINAL_SAMPLE_SIZE_OF_13C_IN_mL=d2(1,1);
CONVERT_13C_SAMPLE_SIZE_TO_13C_IN_PLATE_FACTOR=d2(2,1);


% get experiment time in hours
EXPERIMENT_TIME_IN_HOURS = d(1,8);
% calculate PCV
DOUBLING_TIME = d(1,7);
cell_number_per_well = d(1:3,1).*d(1:3,3);
PCV_of_well_in_uL    = d(1:3,1).*d(1:3,2).*d(1:3,3)*1e-9;
average_PCV = ((-PCV_of_well_in_uL*DOUBLING_TIME)/(EXPERIMENT_TIME_IN_HOURS*log(2)))*(2^(-EXPERIMENT_TIME_IN_HOURS/DOUBLING_TIME)-1)
AVERAGE_PCV = mean(average_PCV);  %in ul, 48 hours
STD_PCV     = std(average_PCV);

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



NUM_OF_RANDOM_SAMPLES = 1000;
uptake_secretion_rate_array=[];
R1_no_cells_array=[];
R1_with_cells_array=[];
R2_no_cells_array=[];
R2_with_cells_array=[];
for(simulation_index=1:NUM_OF_RANDOM_SAMPLES)
    met_list_norm_add_normal_distribution = met_list_norm;
    for i=1:length(met_list_norm)
        normDist = normrnd(repmat(zeros(size(met_list_norm{i}.var)),1,1),repmat(met_list_norm{i}.var.^0.5,1,1));
        met_list_norm_add_normal_distribution{i}.data(1,:) = met_list_norm{i}.data(1,:)+normDist(1,:);
        met_list_norm_add_normal_distribution{i}.data(met_list_norm_add_normal_distribution{i}.data<0)=0;
    end
    PCV = normrnd(AVERAGE_PCV,STD_PCV);
    [R1_no_cells_array(:,end+1) R1_with_cells_array(:,end+1) R2_no_cells_array(:,end+1) R2_with_cells_array(:,end+1) uptake_secretion_rate_array(:,end+1)] = uptake_secretion_rates(met_list_norm_add_normal_distribution, sample_list_short, EXPERIMENT_TIME_IN_HOURS, DOUBLING_TIME, MEDIA_FINAL_SAMPLE_SIZE_OF_13C_IN_mL, CONVERT_13C_SAMPLE_SIZE_TO_13C_IN_PLATE_FACTOR, cell_number_per_well, PCV_of_well_in_uL, PCV, exist_in_bottle_media, fully_labled_in_bottle_media);
end

mean_uptake_secretion = mean(uptake_secretion_rate_array')';
std_uptake_secretion = std(uptake_secretion_rate_array')';
mean_R1_no_cells = mean(R1_no_cells_array')';std_R1_no_cells = std(R1_no_cells_array')';
mean_R2_no_cells = mean(R2_no_cells_array')';std_R2_no_cells = std(R2_no_cells_array')';
mean_R1_with_cells = mean(R1_with_cells_array')';std_R1_with_cells = std(R1_with_cells_array')';
mean_R2_with_cells = mean(R2_with_cells_array')';std_R2_with_cells = std(R2_with_cells_array')';

for i=1:length(met_list_norm)
    fprintf('%s uptake: mean=%s STD=%s\n', met_list_norm{i}.met_name, num2str(mean_uptake_secretion(i)), num2str(std_uptake_secretion(i)));       
end

print_figures(met_list_norm, sample_list_short);

