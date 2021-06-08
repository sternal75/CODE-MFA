% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% The fit been simulated total cellular metabolite isotopic labeling (x-axis) 
% and experimental measurements (y-axis)  
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
load('../mat_files/sensitiviy_analysis_concentration.mat', 'sensitiviy_analysis_concentration');


addpath('../functions/emu') 
addpath('../functions/general') 
addpath('../') 
run ../load_constants;

%     find the index of the best score among all directionalities
best_score = min(directionalities.errors);
index_best_score = find(directionalities.errors==min(directionalities.errors));
best_score_predicted_concentrations     = directionalities.predicted_concentrations(:,index_best_score);    
best_score_predicted_fluxes_fb          = directionalities.predicted_fb_fluxes(:,index_best_score);    
best_score_predicted_net_fluxes         = directionalities.predicted_net_fluxes(:,index_best_score);    


num_of_experiments          = length(WC_known_metabolites_idv);
idv_opt = cell(num_of_experiments,1);
idv_d   = cell(num_of_experiments,1);
idv_known_arr = cell(num_of_experiments,1);
idv_known_mat = cell(num_of_experiments,1);
output_for_excel = cell(num_of_experiments,1);
for (i=1:num_of_experiments)
    idv_known_arr{i} = zeros(0,1);
    idv_known_mat{i}=[];
    for j=1:length(WC_known_metabolites_idv{i})
        idv_known_arr{i} = [idv_known_arr{i};WC_known_metabolites_idv{i}{j}.index_CY;WC_known_metabolites_idv{i}{j}.index_MT];
        idv_known_mat{i} = [idv_known_mat{i};WC_known_metabolites_idv{i}{j}.index_CY WC_known_metabolites_idv{i}{j}.index_MT];
        for(k=1:length(WC_known_metabolites_idv{i}{j}.idv))
            output_for_excel{i}{end+1} = WC_known_metabolites_idv{i}{j}.met_name;
        end
    end    
end
for(i=1:num_of_experiments)
    fluxes_net_after_force_zero{i} = best_score_predicted_net_fluxes;
    fluxes_net_after_force_zero{i}(model_net_fluxes.force_zero_flux{i}==1)=0.1;
    one_direction_flux_after_force_zero{i} = ComputeOneDirectionFluxes(fluxes_net_after_force_zero{i}, best_score_predicted_concentrations, best_score_predicted_fluxes_fb, model_net_fluxes, model_thermodynamics);    
    fcn_name = ['ComputeEmuIDV_Opt' int2str(1)];
    [idv_opt{i} idv_d{i}] = feval(fcn_name,idv{i},one_direction_flux_after_force_zero{i});        
end

measured_labeling_forms_lb      = cell(num_of_experiments,1);
measured_labeling_forms_ub      = cell(num_of_experiments,1);
measured_labeling_forms         = cell(num_of_experiments,1);
WC_convoluted_labeling_forms    = cell(num_of_experiments,1);
for(experiment_index=1:length(idv_known_mat))
    for i=1:size(idv_known_mat{experiment_index},1)
        idv_known_vs_opt=[];
        mass_isotopomer_legend=cell(0);        
        x_CY = idv_known_mat{experiment_index}(i,1);
        x_MT = idv_known_mat{experiment_index}(i,2);
        EMU_indices_CY = find(EMU{experiment_index}.list(:,1)==x_CY);
        EMU_indices_MT = find(EMU{experiment_index}.list(:,1)==x_MT);
        for(j=1:length(WC_known_metabolites_idv{experiment_index}{i}.idv))
            mass_isotopomer_legend{end+1}=sprintf('m+%s',num2str(j-1));
        end

        % calculate cy/mt ratio for metabolie
        x_cy_for_concentrations = WC_known_metabolites_idv{experiment_index}{i}.index_CY_for_concentrations;
        x_mt_for_concentrations = WC_known_metabolites_idv{experiment_index}{i}.index_MT_for_concentrations;
        cy_mt_ratio=CY_WC_VOLUME*exp(best_score_predicted_concentrations(x_cy_for_concentrations))/(CY_WC_VOLUME*exp(best_score_predicted_concentrations(x_cy_for_concentrations))+MT_WC_VOLUME*exp(best_score_predicted_concentrations(x_mt_for_concentrations)));
        
        % output bars for measured and simulated
        idv_known_vs_opt(end+1,:)=WC_known_metabolites_idv{experiment_index}{i}.idv;
        idv_std = sqrt(WC_known_metabolites_idv{experiment_index}{i}.idv_variance)
        if(sum(idv_std>=0.05))
%         if((sum(max(idv_std)>=0.02))&& (sum(max(idv_std)<0.05)))
            continue;
        end
        idv_std(idv_std<0.01) = 0.01;
        if((~isnan(x_CY))&&(~isnan(x_MT)))
            idv_known_vs_opt(end+1,1:length(idv_opt{experiment_index}{EMU_indices_CY(1)}))=cy_mt_ratio*idv_opt{experiment_index}{EMU_indices_CY(1)}+(1-cy_mt_ratio)*idv_opt{experiment_index}{EMU_indices_MT(1)};
            if(WC_known_metabolites_idv{experiment_index}{i}.contaminated_m0)
                m0_diff_from_sim = idv_known_vs_opt(1,1)-idv_known_vs_opt(2,1);
                idv_known_vs_opt(1,1)=idv_known_vs_opt(2,1);
                sum_2_to_end = sum(idv_known_vs_opt(1,2:end));
                idv_known_vs_opt(1,2:end) = idv_known_vs_opt(1,2:end)+(idv_known_vs_opt(1,2:end)/sum_2_to_end)*m0_diff_from_sim;                
                idv_std(1:length(idv_std))=0.01;
            end
        elseif(~isnan(x_CY))
            idv_known_vs_opt(end+1,1:length(idv_opt{experiment_index}{EMU_indices_CY(1)}))=idv_opt{experiment_index}{EMU_indices_CY(1)};           
%            WC_convoluted_labeling_forms{experiment_index} = [WC_convoluted_labeling_forms{experiment_index} idv_opt{experiment_index}{EMU_indices_CY(1)}];
        elseif(~isnan(x_MT))
            idv_known_vs_opt(end+1,1:length(idv_opt{experiment_index}{EMU_indices_MT(1)}))=idv_opt{experiment_index}{EMU_indices_MT(1)};
%            WC_convoluted_labeling_forms{experiment_index} = [WC_convoluted_labeling_forms{experiment_index} idv_opt{experiment_index}{EMU_indices_MT(1)}];
        end        
        measured_labeling_forms_lb{experiment_index} = [measured_labeling_forms_lb{experiment_index} idv_known_vs_opt(1,:)-2*idv_std];
        measured_labeling_forms_lb{experiment_index}(measured_labeling_forms_lb{experiment_index}<0)=0;
        measured_labeling_forms_ub{experiment_index} = [measured_labeling_forms_ub{experiment_index} idv_known_vs_opt(1,:)+2*idv_std];
        measured_labeling_forms_ub{experiment_index}(measured_labeling_forms_ub{experiment_index}>1)=1;
        measured_labeling_forms{experiment_index} = [measured_labeling_forms{experiment_index} idv_known_vs_opt(1,:)];
        WC_convoluted_labeling_forms{experiment_index} = [WC_convoluted_labeling_forms{experiment_index} idv_known_vs_opt(2,:)];                        
    end
end

figure('Position',[1 1 1100 950]);
hold on;
h1=plot(([WC_convoluted_labeling_forms{1};WC_convoluted_labeling_forms{1}]), [measured_labeling_forms_lb{1};measured_labeling_forms_ub{1}],'Color',BLUE_COLOR,'LineWidth',3);
h2=plot(([WC_convoluted_labeling_forms{2};WC_convoluted_labeling_forms{2}]), [measured_labeling_forms_lb{2};measured_labeling_forms_ub{2}],'Color',GREEN_COLOR,'LineWidth',3);
hold off;
xlabel('Simulated labeling form fraction');
ylabel('Measured labeling form fraction');
grid on
set(gcf,'color','w');
set(gca, 'FontSize', 28);
ax=gca;
ax.Box='on';
line([0 1], [0 1], 'color',RED_COLOR,'LineStyle',':','linewidth',2);
legend([h1(1) h2(1)],'U-13C Gln labeling','U-13C Glc labeling  .', 'Location','northwest');

% print gln labeling to excel
xlswrite('temp.xlsx',{'met name'},'Gln labeling','A1');  
xlswrite('temp.xlsx',{'measured'},'Gln labeling','b1');  
xlswrite('temp.xlsx',{'CODE MFA - best fit'},'Gln labeling','c1');  
xlswrite('temp.xlsx',output_for_excel{1}','Gln labeling','A2');  
xlswrite('temp.xlsx',measured_labeling_forms{1}','Gln labeling','b2');  
xlswrite('temp.xlsx',WC_convoluted_labeling_forms{1}','Gln labeling','c2');  

% print glc labeling to excel
xlswrite('temp.xlsx',{'met name'},'Glc labeling','A1');  
xlswrite('temp.xlsx',{'measured'},'Glc labeling','b1');  
xlswrite('temp.xlsx',{'CODE MFA - best fit'},'Glc labeling','c1');  
xlswrite('temp.xlsx',output_for_excel{2}','Glc labeling','A2');  
xlswrite('temp.xlsx',measured_labeling_forms{2}','Glc labeling','b2');  
xlswrite('temp.xlsx',WC_convoluted_labeling_forms{2}','Glc labeling','c2');  

set(gcf,'PaperSize',[50 30]);
s = sprintf('./output_images/3c.pdf');
saveas(gcf, s);
    

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% compute forward/backward fluxes
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
function one_direction_flux = ComputeOneDirectionFluxes(fluxes_net, concentrations, fluxes_fb, model_net_fluxes, model_thermodynamics)
    run ../load_constants;
    one_direction_flux = [];
    one_direction_flux_ind = 1;
    for(i=1:length(model_net_fluxes.is_net_flux))
        if(model_net_fluxes.is_net_flux(i))
            if(model_thermodynamics.thermodynamics_of_reaction_defined(i))
                numerator   = sum(concentrations(model_thermodynamics.product_indexes{i}));
                denominator = sum(concentrations(model_thermodynamics.reactant_indexes{i}));        

                fb_ratio = [((exp(numerator)/exp(denominator))^-1)*exp(-model_thermodynamics.delta_G0(i)/RT)];
                one_direction_flux(one_direction_flux_ind,1)    = fluxes_net(i)*fb_ratio/(fb_ratio-1);
                % backward flux
                one_direction_flux(one_direction_flux_ind+1,1)  = fluxes_net(i)/(fb_ratio-1);                
            else
                one_direction_flux(one_direction_flux_ind,1)    = fluxes_fb(one_direction_flux_ind);
                % backward flux
                one_direction_flux(one_direction_flux_ind+1,1)  = fluxes_fb(one_direction_flux_ind+1);
            end        
            one_direction_flux_ind = one_direction_flux_ind+2;
        else
            one_direction_flux(one_direction_flux_ind,1) = fluxes_net(i);
            one_direction_flux_ind = one_direction_flux_ind+1;
        end
    end     
end

