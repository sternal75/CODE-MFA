% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% Fix flux confidence intervals to be consistent with mass balance
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
clear all; close all;

load('mat_files/sensitiviy_analysis_dG.mat', 'sensitiviy_analysis_dG');
load('mat_files/directionalities.mat', 'directionalities');
load('mat_files/model_thermodynamics.mat','model_thermodynamics');    
load('mat_files/model_net_fluxes.mat','model_net_fluxes');
load('mat_files/model.mat','model');
load('mat_files/EMU.mat','EMU');
load('mat_files/idv.mat','idv');
load('mat_files/WC_known_metabolites_idv.mat','WC_known_metabolites_idv');
load('mat_files/WC_known_metabolites_concentration.mat','WC_known_metabolites_concentration');
load('mat_files/MILP_results.mat','MILP_results');
load('mat_files/iterations_result.mat','iterations_result');
load('mat_files/net_fluxes.mat','net_fluxes');


for(i=1:length(net_fluxes))
    low_flux_our_method(i) = net_fluxes{i}.low_flux;
    high_flux_our_method(i) = net_fluxes{i}.high_flux;
end 


% for constraint in equatoion 17a
Aeq1 = [model_net_fluxes.S(model_net_fluxes.met_extra == 0, :)];
Beq1 = zeros(size(Aeq1,1),1);

Aeq=Aeq1;
Beq = Beq1;


net_fluxes_before_post_processing = net_fluxes;
% find min/max vf-vb
for(i=1:length(net_fluxes))
    % check if lower bound of flux i, fits lower and upper bounds of all
    % fluxes
    lower_bound_for_optimization = low_flux_our_method;
    upper_bound_for_optimization = high_flux_our_method;        
    f=zeros(1,size(Aeq,2));
    f(i)=1;
    [optimization_values,fval,exitflag,output] = linprog(f, [], [], Aeq, Beq, lower_bound_for_optimization, upper_bound_for_optimization);
    if(isempty(optimization_values))
        fprintf('intlinprog is empty\n');
    end
    if (((exitflag ~= 1)&&(exitflag ~= 2)&&(exitflag ~= -3))|(output.constrviolation>0.1))
        fprintf('Error in intlinprog (vf-vb) - %d\n',i);
    else
        results.vf_minus_vb.min(i) = fval;  
        net_fluxes{i}.low_flux = fval;
    end                    

    % check if upper bound of flux i, fits lower and upper bounds of all
    % fluxes    
    lower_bound_for_optimization = low_flux_our_method;
    upper_bound_for_optimization = high_flux_our_method;    
    f=zeros(1,size(Aeq,2));
    f(i)=-1;
    [optimization_values,fval,exitflag,output] = linprog(f, [], [], Aeq, Beq, lower_bound_for_optimization, upper_bound_for_optimization);
    if(isempty(optimization_values))
        fprintf('intlinprog is empty\n');
    end
    if (((exitflag ~= 1)&&(exitflag ~= 2)&&(exitflag ~= -3))|(output.constrviolation>0.1))
        fprintf('Error in intlinprog (vf-vb) - %d\n',i);
    else
        results.vf_minus_vb.max(i) = -fval;  
        net_fluxes{i}.high_flux = -fval;
    end                    
end
save('mat_files/net_fluxes_before_post_processing.mat','net_fluxes_before_post_processing');
save('mat_files/net_fluxes.mat','net_fluxes');

