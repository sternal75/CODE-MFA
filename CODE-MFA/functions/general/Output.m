function Output(final_predicted_flux, deltaG, min_error, model, the_title)
    %OUTPUT Summary of this function goes here
    %   Detailed explanation goes here
    fluxModelImg = imread('modelpic.jpg');
    text_str    = cell(length(final_predicted_flux),1);
    text_str_dG = cell(0);
    [dR,tR] = xlsread('xls_input_files/model_pic_parameters.xlsx');
    position = [];
    position_dG = [];
    arrow_params = [];
    index_in_model_pix_parameters_file = 1;
    for ind=1:length(final_predicted_flux)
        arrow_params=[arrow_params;dR(index_in_model_pix_parameters_file,4:7)];
        text_str{ind} = sprintf('%s=%s', model.rxns{ind}, num2str(final_predicted_flux(ind),'%0.2f'));
        if(~isempty(findstr(model.rxns{ind},'f')))
            position=[position;dR(index_in_model_pix_parameters_file,1:2)];
            text_str_dG{end+1} = sprintf('%s=%s', 'dG', num2str(deltaG(index_in_model_pix_parameters_file),'%0.1f'));
            position_dG=[position_dG;dR(index_in_model_pix_parameters_file,1) dR(index_in_model_pix_parameters_file,2)-11];
        elseif(~isempty(findstr(model.rxns{ind},'b')))
            position=[position;dR(index_in_model_pix_parameters_file,1) dR(index_in_model_pix_parameters_file,2)+11];
            index_in_model_pix_parameters_file = index_in_model_pix_parameters_file+1;
        else
            position=[position;dR(index_in_model_pix_parameters_file,1:2)];
            index_in_model_pix_parameters_file = index_in_model_pix_parameters_file+1;
            text_str_dG{end+1}='';
            position_dG=[position_dG;[-10 -10]];
        end
    end
    
    RGB = insertText(fluxModelImg,position,text_str,'FontSize',12,'BoxColor','white','BoxOpacity',0,'TextColor','black');
%     figure;
    RGB = insertText(RGB,position_dG,text_str_dG','FontSize',12,'BoxColor','white','BoxOpacity',0,'TextColor','black');
    opengl hardware;
    %rmse_str{1} = sprintf('RMSE=%s',num2str(sqrt(min_error),'%0.2f'));
    rmse_str{1} = sprintf('Score=%s',num2str((min_error),'%0.4f'));
    position = [0 0];
    RGB = insertText(RGB,position,rmse_str,'FontSize',10,'BoxColor','green','BoxOpacity',0.4,'TextColor','black');
    imshow(RGB);
    hold on
    quiver(arrow_params(:,1),arrow_params(:,2),arrow_params(:,3),arrow_params(:,4),'color','blue', 'MaxHeadSize', 1,'AutoScale','off', 'LineWidth',2.2);
%     quiver(700,500,100,100,'color','red');
    hold off    
    clear title xlabel ylabel;
    title(the_title);

end

