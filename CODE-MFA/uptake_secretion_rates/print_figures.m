function print_figures(met_list_norm, sample_list_short)
%PRINT_FIGURES Summary of this function goes here
%   Detailed explanation goes here
    legend_arr = cell(0);
    for x=1:50
        legend_arr{x} = sprintf('m+%d', x-1);
    end
    for i=1:length(met_list_norm)
        mat_bar = met_list_norm{i}.data(met_list_norm{i}.mass_isotopomer_list, :)';
        figure;
        sgtitle(sprintf('%s',met_list_norm{i}.met_name), 'FontSize', 26);
        for(j=1:size(mat_bar,1))
            
            current_bar = zeros(2,size(mat_bar,2));
            current_bar(1,:)=mat_bar(j,:);


            mat_bar_total =  mat_bar;
            subplot(2,4,j);
            bar_label( current_bar, current_bar);


            set(gca, 'Xlim', [0.5 1.5]);

            title(sprintf('%s',strrep(sample_list_short{j},'_',' ')), 'FontSize', 16);
            legend(legend_arr(met_list_norm{i}.mass_isotopomer_list), 'Location', 'NorthEastOutside');

            set(gca, 'XTick', [1:1:length(sample_list_short)]); 
            set(gca, 'XTickLabel', '');
        %     ylabel('% labeling','fontsize',20);
            set(gca, 'Ylim', [0 1]);
            set(gca, 'FontSize', 14);
        end
    end
end

