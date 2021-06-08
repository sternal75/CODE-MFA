function [R1_no_cells_array R1_with_cells_array R2_no_cells_array R2_with_cells_array uptake_secretion_rate_array] = uptake_secretion_rates(met_list_norm, sample_list_short, EXPERIMENT_TIME_IN_HOURS, DOUBLING_TIME, MEDIA_FINAL_SAMPLE_SIZE_OF_13C_IN_mL, CONVERT_13C_SAMPLE_SIZE_TO_13C_IN_PLATE_FACTOR, cell_number_per_well, PCV_of_well_in_uL, PCV, exist_in_bottle_media, fully_labled_in_bottle_media)

    % get labeled and unlabaled data
    R=cell(0);
    for i=1:length(met_list_norm)
        unlabeled   = met_list_norm{i}.data(1,:);
%         labeled     = sum(met_list_norm{i}.data(2:end,:));
        labeled     = 1-unlabeled;
        % convert to matrix with first row contains sample media and second row
        % contains media with cells
        unlabeled=[unlabeled(1:length(unlabeled)/2);unlabeled(length(unlabeled)/2+1:length(unlabeled))];
        labeled=[labeled(1:length(labeled)/2);labeled(length(labeled)/2+1:length(labeled))];
        R{end+1}.R1=labeled(:,1)./unlabeled(:,1);
        R{end}.R2=labeled(:,2:end)./unlabeled(:,2:end);
        R2_converted = R{end}.R2';
        R2_converted(R2_converted<1)=1./R2_converted(R2_converted<1);
        [x R{end}.best_dilution_position]=min(R2_converted);
        x=x';R{end}.best_dilution_position=R{end}.best_dilution_position';
        best_dilution_position_no_cells     = R{end}.best_dilution_position(1);
        best_dilution_position_with_cells   = R{end}.best_dilution_position(2);
        best_dilution_no_cells      = R{end}.R2(1,best_dilution_position_no_cells);
        best_dilution_with_cells    = R{end}.R2(2,best_dilution_position_with_cells);
        if((best_dilution_no_cells>100)&&(exist_in_bottle_media(i)))
            error(sprintf('ERROR - outside linear range of media without cells for %s', met_list_norm{i}.met_name));
        end
        if(best_dilution_with_cells>100)
            error(sprintf('ERROR - outside linear range of media with cells for %s', met_list_norm{i}.met_name));
        end    
    end

    dilution_values=[];
    for i=1:length(sample_list_short)
        spilt_sample_list_str=strsplit(sample_list_short{i},'_');
        last_part_of_sample_list_str = spilt_sample_list_str{end};
        dilution_str = strsplit(last_part_of_sample_list_str,'mMIS');
        dilution_str = dilution_str{1};
        if(isempty(str2num(dilution_str)))
            % do nothing 
        else
            dilution_values(end+1) = str2num(dilution_str);
        end
    end
    dilution_values=[dilution_values(1:length(dilution_values)/2);dilution_values(length(dilution_values)/2+1:length(dilution_values))];
    uptake_secretion_rate_array=[];
    R1_no_cells_array = [];
    R1_with_cells_array = [];
    R2_no_cells_array = [];
    R2_with_cells_array = [];
    for i=1:length(met_list_norm)
        R1_no_cells = R{i}.R1(1);
        R1_with_cells = R{i}.R1(2);
        R2_no_cells = R{i}.R2(1,R{i}.best_dilution_position(1));
        R2_with_cells = R{i}.R2(2,R{i}.best_dilution_position(2));
        if(exist_in_bottle_media(i))
            if(fully_labled_in_bottle_media(i))
                metabolite_bottle_media_13c_umole  = R2_no_cells*dilution_values(1,R{i}.best_dilution_position(1))*MEDIA_FINAL_SAMPLE_SIZE_OF_13C_IN_mL;
                metabolite_with_cells_13c_umole    = R2_with_cells*dilution_values(2,R{i}.best_dilution_position(2))*MEDIA_FINAL_SAMPLE_SIZE_OF_13C_IN_mL;
            else            
                metabolite_bottle_media_13c_umole  = ((R2_no_cells+R2_no_cells/R1_no_cells)/(1-R2_no_cells/R1_no_cells))*dilution_values(1,R{i}.best_dilution_position(1))*MEDIA_FINAL_SAMPLE_SIZE_OF_13C_IN_mL;
                metabolite_with_cells_13c_umole    = ((R2_with_cells+R2_with_cells/R1_with_cells)/(1-R2_with_cells/R1_with_cells))*dilution_values(2,R{i}.best_dilution_position(2))*MEDIA_FINAL_SAMPLE_SIZE_OF_13C_IN_mL;    
            end
        else
            metabolite_bottle_media_13c_umole  = 0;
            metabolite_with_cells_13c_umole    = ((R2_with_cells+R2_with_cells/R1_with_cells)/(1-R2_with_cells/R1_with_cells))*dilution_values(2,R{i}.best_dilution_position(2))*MEDIA_FINAL_SAMPLE_SIZE_OF_13C_IN_mL;    
        end

        uptake_of_sample_in_umole       = metabolite_bottle_media_13c_umole-metabolite_with_cells_13c_umole;
        uptake_of_plate_in_umole        = uptake_of_sample_in_umole*CONVERT_13C_SAMPLE_SIZE_TO_13C_IN_PLATE_FACTOR;
        final_uptake_in_mMole_per_hour  = uptake_of_plate_in_umole/EXPERIMENT_TIME_IN_HOURS/PCV*1000;
%         fprintf('%s uptake = %s\n', met_list_norm{i}.met_name, num2str(final_uptake_in_mMole_per_hour));   
        uptake_secretion_rate_array(i,1) = final_uptake_in_mMole_per_hour;
        R2_with_cells = R{i}.R2(2,R{i}.best_dilution_position(2));
        
        R1_no_cells_array(i,1)   = R1_no_cells;
        R1_with_cells_array(i,1) = R1_with_cells;
        R2_no_cells_array(i,1)   = R2_no_cells;
        R2_with_cells_array(i,1) = R2_with_cells;
    end
end

