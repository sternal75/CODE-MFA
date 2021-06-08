% file_name = 'input';

[met_list_norm_neg met_name_arr_neg sample_list_short] = ProcessMavenIsotopicLabel(sprintf('data/%s.xls', file_name), 3);
[met_list_norm_pos met_name_arr_pos sample_list_short] = ProcessMavenIsotopicLabel(sprintf('data/%s.xls', file_name), 3);    
[dM,tM] = xlsread(sprintf('data/%s.xls', file_name),'Contaminated m0');
contaminated_m0_mets = tM;


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
%         met_list_norm{end}.met_name = sprintf('%s (neg)', met_list_norm{end}.met_name);
    end

    for i=1:length(pos_only)
        pos_index = strmatch(pos_only{i}, met_name_arr_pos, 'exact');    
        met_list_norm{end+1} = met_list_norm_pos{pos_index};
%         met_list_norm{end}.met_name = sprintf('%s (pos)', met_list_norm{end}.met_name);
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
mkdir(file_name);

for i=1:length(met_list_norm)
figure;    
%    hbar = bar(met_list_norm{i}.data', 'stacked');
    mat_bar = met_list_norm{i}.data(met_list_norm{i}.mass_isotopomer_list, :)';
%     mat_bar(2,:)=0;
    
%     last_bar = zeros(size(mat_bar(1:2,:)));
last_bar = zeros(2,size(mat_bar,2));
    last_bar(1,:)=mat_bar(end,:);

    
    mat_bar_total =  mat_bar;
%     mat_bar_total = mat_bar .* repmat(met_list_norm{i}.data_total',1,size(mat_bar,2));

    i
%     figure;
    
    if(i>12)
        if(i==13)
            figure;
        end
%         subplot(3,4,i-12);
    else
%         subplot(3,4,i);
    end
    bar_to_plot = zeros(2,50);
    bar_to_plot(1,met_list_norm{i}.mass_isotopomer_list) = last_bar(1,:);
    bar_label( last_bar, last_bar, met_list_norm{i}.mass_isotopomer_list);

    
    set(gca, 'Xlim', [0.5 1.5]);
    
    met_name = strrep(met_list_norm{i}.met_name,'_',' ');
    title(sprintf('%s',met_name), 'FontSize', 18);
    legend(legend_arr(met_list_norm{i}.mass_isotopomer_list), 'Location', 'NorthEastOutside');

    set(gca, 'XTick', [1:1:length(sample_list_short)]); 
    set(gca, 'XTickLabel', '');
%     ylabel('% labeling','fontsize',20);
    set(gca, 'Ylim', [0 1]);
    set(gca, 'FontSize', 17);
    
    s = sprintf('./%s/%s.jpg', file_name, met_list_norm{i}.met_name);
    saveas(gcf, s);
end

% add a flag for metabolites that are m+0 contaminated
% the contaminated m+0 will be considered later on in the optimization and
% the m+0 will not be taken into account in the labeling fitting process
for i=1:length(met_list_norm)
    index = strmatch(met_list_norm{i}.met_name, contaminated_m0_mets, 'exact');    
    if(isempty(index))
        met_list_norm{i}.contaminated_m0=0;
    else
        met_list_norm{i}.contaminated_m0=1;
    end
end



