% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% Find all possible flux directionalities with MILP
% See paper equation 17 and 21
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
function [results] = FindAllFluxDirectionalities(model_net_fluxes, model_thermodynamics,known_net_flux_directions)

load_constants

ts=clock;
rand('seed',(ts(6)*10000));


vf_length   = length(model_net_fluxes.rxns);
vb_length   = length(model_net_fluxes.rxns);
mets_length = length(model_thermodynamics.mets);
yf_length   = length(model_net_fluxes.rxns);
yb_length   = length(model_net_fluxes.rxns);

positive_direction_lower_buond = model_net_fluxes.positive_direction_lb;
negative_direction_lower_buond = model_net_fluxes.negative_direction_lb;
vf_lower_buond_zero = positive_direction_lower_buond;
vf_lower_buond_zero ((model_net_fluxes.is_net_flux) & (vf_lower_buond_zero==0.001))=0;
vb_lower_buond_zero = negative_direction_lower_buond;
vb_lower_buond_zero (:)=0;

positive_direction_upper_buond = model_net_fluxes.positive_direction_ub;
negative_direction_upper_buond = model_net_fluxes.negative_direction_ub;
vf_upper_bound = positive_direction_upper_buond;
vb_upper_bound = negative_direction_upper_buond;
vb_upper_bound(model_net_fluxes.is_net_flux==0) = 0;



% for constraint in equatoion 17d
yf_diag = eye(length(model_net_fluxes.rxns));
yb_diag = eye(length(model_net_fluxes.rxns));
yf_lower_bound=ones(length(yf_diag),1);
yf_lower_bound(model_net_fluxes.is_net_flux==1)=0;

yf_upper_bound=ones(length(yf_diag),1);
yb_lower_bound=zeros(length(yb_diag),1);
yb_upper_bound=zeros(length(yb_diag),1);
yb_upper_bound(model_net_fluxes.is_net_flux==1)=1;

Aeq2 = [zeros(size(yf_diag,1),vf_length+vb_length+mets_length) yf_diag yb_diag];
Beq2 = [ones(size(Aeq2,1),1)];


% these directionalities were taken from regular MFA, where only the fowrard
% net flux was within the confidence intervals
for(i=1:length(known_net_flux_directions))
    if(known_net_flux_directions(i)==1)
        yf_lower_bound(i)=1;
        yb_upper_bound(i)=0;
    elseif(known_net_flux_directions(i)==0)
        yf_upper_bound(i)=0;
        yb_lower_bound(i)=1;
    end
end


% for constraint in equatoion 17a
Aeq1 = [model_net_fluxes.S(model_net_fluxes.met_extra == 0, :)....
       -model_net_fluxes.S(model_net_fluxes.met_extra == 0, :)....
       zeros(size(model_net_fluxes.S(model_net_fluxes.met_extra == 0), 1),mets_length+yf_length+yb_length)];
Beq1 = zeros(size(Aeq1,1),1);

% add equality constraints for fluxes
AeqFluxEqualityConstraints = [model_net_fluxes.equality_constraints zeros(size(model_net_fluxes.equality_constraints,1),vb_length+mets_length+yf_length+yb_length)];
BeqFluxEqualityConstraints = [zeros(size(model_net_fluxes.equality_constraints,1),1)];
AeqFluxEqualityConstraints = [AeqFluxEqualityConstraints;zeros(size(model_net_fluxes.equality_constraints,1),vf_length) model_net_fluxes.equality_constraints zeros(size(model_net_fluxes.equality_constraints,1),mets_length+yf_length+yb_length)];
BeqFluxEqualityConstraints = [BeqFluxEqualityConstraints;zeros(size(model_net_fluxes.equality_constraints,1),1)];


Aeq=[Aeq1;Aeq2;AeqFluxEqualityConstraints];
Beq=[Beq1;Beq2;BeqFluxEqualityConstraints];



% for constraint in equatoion 17f, 17g
vf_diag = eye(length(model_net_fluxes.rxns));
vb_diag = eye(length(model_net_fluxes.rxns));
A3right = [vf_diag zeros(size(vb_diag)) zeros(size(vf_diag,2),length(model_thermodynamics.mets)) -diag(vf_upper_bound) zeros(size(yb_diag))];
A4right = [zeros(size(vf_diag)) vb_diag zeros(size(vf_diag,2),length(model_thermodynamics.mets)) zeros(size(yf_diag))  -diag(vb_upper_bound) ];
A3right(model_net_fluxes.is_net_flux==0,:)=0;
A4right(model_net_fluxes.is_net_flux==0,:)=0;
A34right = [A3right;A4right];
b34right = zeros(size(A34right,1),1);

A3left = [-vf_diag zeros(size(vb_diag)) zeros(size(vf_diag,2),length(model_thermodynamics.mets)) diag(positive_direction_lower_buond) zeros(size(yb_diag))];
A4left = [zeros(size(vf_diag)) -vb_diag zeros(size(vf_diag,2),length(model_thermodynamics.mets)) zeros(size(yf_diag))  diag(negative_direction_lower_buond) ];
A3left(model_net_fluxes.is_net_flux==0,:)=0;
A4left(model_net_fluxes.is_net_flux==0,:)=0;
A34left = [A3left;A4left];
b34left = zeros(size(A34left,1),1);

A34 = [A34right;A34left];
b34 = [b34right;b34left];


% for constraint in equatoion 17e
C=1000;
A1 = zeros(length(model_thermodynamics.rxns), vf_length+vb_length+mets_length+yf_length+yb_length);
b1 = zeros(length(model_thermodynamics.rxns),1);
for(i=1:length(model_thermodynamics.rxns))
    A1(i,vf_length+vb_length+model_thermodynamics.product_indexes{i})  = RT;
    A1(i,vf_length+vb_length+model_thermodynamics.reactant_indexes{i}) = -RT;
    if(~isnan(model_thermodynamics.delta_G_high(i)))
        b1(i) = -model_thermodynamics.delta_G0(i)+model_thermodynamics.delta_G_high(i);
    else
        A1(i,vf_length+vb_length+mets_length+i) = C;
        b1(i) = -model_thermodynamics.delta_G0(i)-dG_MIN_DISTANCE_FROM_ZERO+C;   
    end
end

A2 = zeros(length(model_thermodynamics.rxns), vf_length+vb_length+mets_length+yf_length+yb_length);
b2 = zeros(length(model_thermodynamics.rxns),1);
for(i=1:length(model_thermodynamics.rxns))
    A2(i,vf_length+vb_length+model_thermodynamics.product_indexes{i})  = -RT;
    A2(i,vf_length+vb_length+model_thermodynamics.reactant_indexes{i}) = RT;
    if(~isnan(model_thermodynamics.delta_G_low(i)))
        b2(i) = model_thermodynamics.delta_G0(i)-model_thermodynamics.delta_G_low(i);             
    else
        A2(i,vf_length+vb_length+mets_length+yf_length+i) = C;
        b2(i) = model_thermodynamics.delta_G0(i)-dG_MIN_DISTANCE_FROM_ZERO+C; 
    end        
end
A = [A1;A2;A34];
b = [b1;b2;b34];



% add lb/ub constratins for co-factor ratios
A5 = zeros(length(model_thermodynamics.co_factors), vf_length+vb_length+mets_length+yf_length+yb_length);
A6 = zeros(length(model_thermodynamics.co_factors), vf_length+vb_length+mets_length+yf_length+yb_length);
for(i=1:length(model_thermodynamics.co_factors))
    A5(i,vf_length+vb_length+model_thermodynamics.co_factors{i}.indices(1)) = 1;
    A5(i,vf_length+vb_length+model_thermodynamics.co_factors{i}.indices(2)) = -1;
    b5(i,1)=log(model_thermodynamics.co_factors{i}.ratio_ub);
    A6(i,vf_length+vb_length+model_thermodynamics.co_factors{i}.indices(1)) = -1;
    A6(i,vf_length+vb_length+model_thermodynamics.co_factors{i}.indices(2)) = 1;    
    b6(i,1)=-log(model_thermodynamics.co_factors{i}.ratio_lb);
end
A=[A;A5;A6];
b=[b;b5;b6];
 

milp_first_integer_parameter_in = vf_length+vb_length+mets_length+1;
Intcon = [milp_first_integer_parameter_in:vf_length+vb_length+mets_length+yf_length+yb_length];
f=ones(1,size(Aeq,2));
lower_bound_for_milp = ([vf_lower_buond_zero;vb_lower_buond_zero;model_thermodynamics.mets_lb;yf_lower_bound;yb_lower_bound]);
upper_bound_for_milp = ([vf_upper_bound;vb_upper_bound;model_thermodynamics.mets_ub;yf_upper_bound;yb_upper_bound]);
yf_all_prev_result_array=[];
optimization_values_prev_result_array=[];
counter=0;


A_without_found_directionalities = A;
b_without_found_directionalities = b;
% lower and upper bounds for fluxes and concentrations
lower_bound_for_milp_without_found_directionalities = lower_bound_for_milp;
upper_bound_for_milp_without_found_directionalities = upper_bound_for_milp;

% Find all possible vectors of flux directionalities
while(1)                                
    [optimization_values,fval,exitflag,output] = intlinprog(f, Intcon, A, b, Aeq, Beq, lower_bound_for_milp, upper_bound_for_milp);
    if(isempty(optimization_values))
        break;
    end
    yf_current_result   = optimization_values(vf_length+vb_length+mets_length+1:vf_length+vb_length+mets_length+yf_length)';
    yf_current_result   = round(yf_current_result);
    A9 = [zeros(1,vf_length+vb_length+mets_length) -1+yf_current_result -yf_current_result];
    A(end+1,:)  = A9;
    b(end+1)    = -0.5;   
    counter = counter+1;
    if ((exitflag ~= 1)&&(exitflag ~= 2)&&(exitflag ~= -3))
        fprintf('Error in fmincon\n');
        continue;
    end
    
    lower_bound_for_milp_without_found_directionalities(vf_length+vb_length+mets_length+1:vf_length+vb_length+mets_length+yf_length) = yf_current_result;
    upper_bound_for_milp_without_found_directionalities(vf_length+vb_length+mets_length+1:vf_length+vb_length+mets_length+yf_length) = yf_current_result;
    [optimization_values,fval,exitflag,output] = intlinprog(f, Intcon, A_without_found_directionalities, b_without_found_directionalities, Aeq, Beq, lower_bound_for_milp_without_found_directionalities, upper_bound_for_milp_without_found_directionalities);
    if(isempty(optimization_values))
        continue;
    end    
    yf_all_prev_result_array    = [yf_all_prev_result_array yf_current_result'];
    optimization_values_prev_result_array = [optimization_values_prev_result_array optimization_values];
    counter
end
results.all_directionalities = yf_all_prev_result_array;
