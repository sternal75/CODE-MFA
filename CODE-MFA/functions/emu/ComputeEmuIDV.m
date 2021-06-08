function [idv result_idv_d cycle_error] = ComputeEmuIDV(EMU, idv, idv_d_emu_ind, flux, fcn_name)

output_file = [];
output_file{end+1} = 'function [idv result_idv_d cycle_error] = ComputeEmuIDV(idv, flux)\n';


tic_begin = tic;
ticks = 0;

n = length(EMU.list);
cycle_error = 0;
rxn_num = length(flux);

result_idv_d_arr = zeros(n, 1);
for i=1:length(idv_d_emu_ind)
    result_idv_d_arr(find(EMU.list(:,1)==idv_d_emu_ind(i)))=idv_d_emu_ind(i);
end
idv_d_emu_num = length(idv_d_emu_ind);
% result_idv_d_arr = zeros(n, 1);
% result_idv_d_arr(idv_d_emu_ind) = [1:idv_d_emu_num];

result_idv_d = cell(idv_d_emu_num, 1);
for i=idv_d_emu_ind(:)'
%     result_idv_d{i} = sparse(rxn_num, EMU.size(i)+1);
end

idv_d_initial = cell(n, 1);
output_file{end+1} = sprintf('idv_d_initial = cell(%d, 1);\n', n);

for i=1:n
%    if (~isempty(idv{i}))
%        idv_d_initial{i} = zeros(1, length(idv{i}));
%   end
    output_file{end+1} = sprintf('idv_d_initial{%d} = zeros(1, %d);\n', i, EMU.size(i)+1);
    idv_d_initial{i} = zeros(1, EMU.size(i)+1);
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%v_sym = sym('v_sym', [1 length(flux)]);
%idv_sym = [];
%for i=1:length(idv)
%    t = sym(['idv_sym_' num2str(i) '_'], [1 7]);
%    idv_sym = [idv_sym t];
%end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cluster_num = length(EMU.cluster);
left_mat = cell(cluster_num,1);
right_mat1 = cell(cluster_num,1);
right_mat_idv1 = cell(cluster_num,1);
right_mat2 = cell(cluster_num,1);
right_mat_idv2 = cell(cluster_num,1);

output_file{end+1} = sprintf('left_mat = cell(%d,1);\n', cluster_num);
output_file{end+1} = sprintf('right_mat1 = cell(%d,1);\n', cluster_num);
output_file{end+1} = sprintf('right_mat_idv1 = cell(%d,1);\n', cluster_num);
output_file{end+1} = sprintf('right_mat2 = cell(%d,1);\n', cluster_num);
output_file{end+1} = sprintf('right_mat_idv2 = cell(%d,1);\n', cluster_num);


for i=1:cluster_num
    
   clust_ind = EMU.cluster_order(i);
   
   if (isempty(EMU.cluster{clust_ind}.left_mat))
       %fprintf('Cluster #%d: external EMUs\n', clust_ind);
       continue;
   end

   output_file{end+1} = '%%%%***********************\n';
   output_file{end+1} = sprintf('%%%%Cluster %d\n', i);
   
   cluster_EMU_size = EMU.size(EMU.cluster{clust_ind}.EMU_ind(1));
   left_mat{i} = sparse(   EMU.cluster{clust_ind}.left_mat(:,1), ...
                        EMU.cluster{clust_ind}.left_mat(:,2), ...
                        flux(abs(EMU.cluster{clust_ind}.left_mat(:,3))) .* sign(EMU.cluster{clust_ind}.left_mat(:,3)), ...
                        EMU.cluster{clust_ind}.EMU_num, EMU.cluster{clust_ind}.EMU_num);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 %   s_left = cell(size(left_mat{i},1),1);
    
%    for x=1:length(s_left)
%        s_left{x} = '0';
%    end
%    for x=1:size(EMU.cluster{clust_ind}.left_mat,1)
%        s_left{  EMU.cluster{clust_ind}.left_mat(x,1) } = sprintf('%s + idv_sym_%d*v_sym_%d', s_left{  EMU.cluster{clust_ind}.left_mat(x,1)+0 }, EMU.cluster{clust_ind}.left_mat(x,2)+0, abs(EMU.cluster{clust_ind}.left_mat(x,3)+0));
%    end
    output_file{end+1} = sprintf('left_mat{%d} = zeros(%d, %d);\n', i, EMU.cluster{clust_ind}.EMU_num, EMU.cluster{clust_ind}.EMU_num);
    for x=1:size(EMU.cluster{clust_ind}.left_mat,1)
        output_file{end+1} = sprintf('left_mat{%d}(%d, %d) = left_mat{%d}(%d, %d) + %d * flux(%d);\n', i, EMU.cluster{clust_ind}.left_mat(x,1)+0, EMU.cluster{clust_ind}.left_mat(x,2)+0, i, EMU.cluster{clust_ind}.left_mat(x,1)+0, EMU.cluster{clust_ind}.left_mat(x,2)+0, sign(EMU.cluster{clust_ind}.left_mat(x,3))+0, abs(EMU.cluster{clust_ind}.left_mat(x,3))+0);
    end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%                    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%                    

   right_mat1{i} = sparse( EMU.cluster{clust_ind}.right_mat1(:,1), ...
                        EMU.cluster{clust_ind}.right_mat1(:,2), ...
                        flux(abs(EMU.cluster{clust_ind}.right_mat1(:,3))) .* sign(EMU.cluster{clust_ind}.right_mat1(:,3)), ...
                        EMU.cluster{clust_ind}.EMU_num, EMU.cluster{clust_ind}.prev_EMU_num);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%               
    %idv_right1_sym = sym('idv_right1_sym', [1 EMU.cluster{clust_ind}.EMU_num]);
    
    %s_right1 = cell(size(right_mat1{i},1),1);
    
    %for x=1:length(s_right1)
%        s_right1{x} = 'idv_right1_sym(x) = 0';
%    end
%    for x=1:size(EMU.cluster{clust_ind}.left_mat,1)
%        s_right1{  EMU.cluster{clust_ind}.left_mat(x,1) } = sprintf('%s + EMU_sym_%d*v_sym_%d', s_right1{  EMU.cluster{clust_ind}.right_mat1(x,1)+0 }, EMU.cluster{clust_ind}.prev_EMU_ind(EMU.cluster{clust_ind}.right_mat1(x,2))+0, abs(EMU.cluster{clust_ind}.right_mat1(x,3)+0));
%    end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%                    
    output_file{end+1} = sprintf('right_mat1{%d} = zeros(%d, %d);\n', i, EMU.cluster{clust_ind}.EMU_num, EMU.cluster{clust_ind}.prev_EMU_num);
    for x=1:size(EMU.cluster{clust_ind}.right_mat1,1)
        output_file{end+1} = sprintf('right_mat1{%d}(%d, %d) = right_mat1{%d}(%d, %d) + %d * flux(%d);\n', i, EMU.cluster{clust_ind}.right_mat1(x,1)+0, EMU.cluster{clust_ind}.right_mat1(x,2)+0, i, EMU.cluster{clust_ind}.right_mat1(x,1)+0, EMU.cluster{clust_ind}.right_mat1(x,2)+0, sign(EMU.cluster{clust_ind}.right_mat1(x,3))+0, abs(EMU.cluster{clust_ind}.right_mat1(x,3))+0);
    end
    
    right_mat_idv1{i} = sparse(EMU.cluster{clust_ind}.prev_EMU_num, cluster_EMU_size+1);
    output_file{end+1} = sprintf('right_mat_idv1{%d} = zeros(%d, %d);\n', i, EMU.cluster{clust_ind}.prev_EMU_num, cluster_EMU_size+1);
    if ~isempty(right_mat1{i})
        for z=1:EMU.cluster{clust_ind}.prev_EMU_num
            right_mat_idv1{i}(z, :) = idv{EMU.cluster{clust_ind}.prev_EMU_ind(z)}';
            output_file{end+1} = sprintf('right_mat_idv1{%d}(%d, :) = idv{%d};\n', i, z, EMU.cluster{clust_ind}.prev_EMU_ind(z));
        end
    end
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%                    
%    if ~isempty(right_mat1{i})
%        s1 = sprintf('right_mat_idv1{%d} = idv([', i);
%        s2 = sprintf('%d ', EMU.cluster{clust_ind}.prev_EMU_ind');
%        s3 = '], 1:);';
%        
%        output_file{end+1} = [s1 s2 s3];
%    end
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%                    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%                    
   right_mat2{i} = sparse( EMU.cluster{clust_ind}.right_mat2(:,1), ...
                        EMU.cluster{clust_ind}.right_mat2(:,2), ...
                        flux(abs(EMU.cluster{clust_ind}.right_mat2(:,3))) .* sign(EMU.cluster{clust_ind}.right_mat2(:,3)), ...
                        EMU.cluster{clust_ind}.EMU_num, EMU.cluster{clust_ind}.prev_EMU_pair_num);
                    
    output_file{end+1} = sprintf('right_mat2{%d} = zeros(%d, %d);\n', i, EMU.cluster{clust_ind}.EMU_num, EMU.cluster{clust_ind}.prev_EMU_pair_num);
    for x=1:size(EMU.cluster{clust_ind}.right_mat2,1)
        output_file{end+1} = sprintf('right_mat2{%d}(%d, %d) = right_mat2{%d}(%d, %d) + %d * flux(%d);\n', i, EMU.cluster{clust_ind}.right_mat2(x,1)+0, EMU.cluster{clust_ind}.right_mat2(x,2)+0, i, EMU.cluster{clust_ind}.right_mat2(x,1)+0, EMU.cluster{clust_ind}.right_mat2(x,2)+0,sign(EMU.cluster{clust_ind}.right_mat2(x,3))+0, abs(EMU.cluster{clust_ind}.right_mat2(x,3))+0);
    end
                    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%                    
   right_mat_idv2{i} = sparse(EMU.cluster{clust_ind}.prev_EMU_pair_num, cluster_EMU_size+1);
   output_file{end+1} = sprintf('right_mat_idv2{%d} = zeros(%d, %d);\n', i, EMU.cluster{clust_ind}.prev_EMU_pair_num, cluster_EMU_size+1);
   
   if ~isempty(right_mat2{i})
       [aa, bb, cc] = find( EMU.cluster{clust_ind}.prev_EMU_pair_cluster_ind );
       for z = 1:length(aa)
           EMU1 = aa(z);
           EMU2 = bb(z);       
           EMU_pair_ind = cc(z);      
           
           m = CreateCauchyMat( EMU.size(EMU1)+1, EMU.size(EMU2)+1);    % Move to pre-processing
           [a,b,c] = find(m);
           m = sparse(a, b, idv{EMU1}(c));

           v = m*idv{EMU2}';
           right_mat_idv2{i}(EMU_pair_ind, :) = v;
           
           %%%
           for x=1:(cluster_EMU_size+1)
               s = '0';
               for k=1:EMU.size(EMU1)+1
                   q = x - k;
                   if (q>=0 & q < EMU.size(EMU2)+1)
                       s = sprintf('%s + idv{%d}(%d) * idv{%d}(%d)', s, EMU1, k, EMU2, q+1);
                   end
               end
               output_file{end+1} = sprintf('right_mat_idv2{%d}(%d, %d) = %s;\n',i, EMU_pair_ind, x, s); 
           end
       end
   end


    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%                    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%                    
   %m1 = inv(left_mat);
   m2 = right_mat1{i}*right_mat_idv1{i};
   m3 = right_mat2{i}*right_mat_idv2{i};
   m4 = m2+m3;

   if(any(any(m3)))
      output_file{end+1} = sprintf('m2 = right_mat1{%d}*right_mat_idv1{%d};\n', i, i);
      output_file{end+1} = sprintf('m3 = right_mat2{%d}*right_mat_idv2{%d};\n', i, i);  
      output_file{end+1} = sprintf('m4 = m2+m3;\n');     
   else
      output_file{end+1} = sprintf('m2 = right_mat1{%d}*right_mat_idv1{%d};\n', i, i);
      output_file{end+1} = sprintf('m4 = m2;\n');            
   end
  
   t1 = tic;
   m = left_mat{i}\m4;   
   output_file{end+1} = sprintf('m = left_mat{%d}\\\\m4;\n', i);     

   t2 = toc(t1);
   ticks = ticks + t2;
   
   %fprintf('Cluster #%d:\t', clust_ind);
   for x=1:EMU.cluster{clust_ind}.EMU_num
       idv{EMU.cluster{clust_ind}.EMU_ind(x)} = m(x,:);

      output_file{end+1} = sprintf('idv{%d} = m(%d,:);\n', EMU.cluster{clust_ind}.EMU_ind(x), x);     
       
    if(isnan(idv{EMU.cluster{clust_ind}.EMU_ind(x)}(1)))
        alon=1
    end       
       %fprintf('%s = %s\t', EMU.name{EMU.cluster{clust_ind}.EMU_ind(x)}, DispIDV(m(x,:)));
   end
   %fprintf('\n');
end


   
   
   
%%%%%%%%%%%%%%
% Compute derivatives


for i=1:length(EMU.cluster)
   output_file{end+1} = sprintf('left_mat_inv{%d} = inv(left_mat{%d});\n', i, i);
end


for i=1:length(EMU.cluster)
   clust_ind = EMU.cluster_order(i);
   cluster_EMU_size = EMU.size(EMU.cluster{clust_ind}.EMU_ind(1));  

   left_mat_idv{i} = sparse(EMU.cluster{clust_ind}.EMU_num, cluster_EMU_size+1);
   output_file{end+1} = sprintf('left_mat_idv{%d} = zeros(%d, %d);\n', i, EMU.cluster{clust_ind}.EMU_num, cluster_EMU_size+1);
   for z=1:EMU.cluster{clust_ind}.EMU_num
       left_mat_idv{i}(z, :) = idv{EMU.cluster{clust_ind}.EMU_ind(z)};
       output_file{end+1} = sprintf('left_mat_idv{%d}(%d, :) = idv{%d};\n', i, z, EMU.cluster{clust_ind}.EMU_ind(z));
   end
end
    

toc(tic_begin)

for x = 1:rxn_num
    idv_d = idv_d_initial;
    output_file{end+1} = sprintf('idv_d = idv_d_initial;\n'); 
 

    for i=1:length(EMU.cluster)
        
       clust_ind = EMU.cluster_order(i);
       
       cluster_EMU_size = EMU.size(EMU.cluster{clust_ind}.EMU_ind(1));  
       
       if (clust_ind == 14)
           if (x==1)
               aa=1;
           end
       end
       
       if (isempty(EMU.cluster{clust_ind}.left_mat))
           continue;
       end

               
       
       output_file{end+1} = '%%%%***********************\n';
       output_file{end+1} = sprintf('%%%%Reaction %d Derivative; Cluster %d\n', x, i);
       
        right_mat_idv_d1 = sparse(EMU.cluster{clust_ind}.prev_EMU_num, cluster_EMU_size+1);
        output_file{end+1} = sprintf('right_mat_idv_d1 = zeros(%d, %d);\n', EMU.cluster{clust_ind}.prev_EMU_num, cluster_EMU_size+1); 

        if ~isempty(right_mat1{i})
            for z=1:EMU.cluster{clust_ind}.prev_EMU_num
                if(any(any(idv_d{EMU.cluster{clust_ind}.prev_EMU_ind(z)})))
                    right_mat_idv_d1(z, :) = idv_d{EMU.cluster{clust_ind}.prev_EMU_ind(z)}';
                    output_file{end+1} = sprintf('right_mat_idv_d1(%d, :) = idv_d{%d};\n', z, EMU.cluster{clust_ind}.prev_EMU_ind(z));
                end
            end
        end

       right_mat_idv_d2 = sparse(EMU.cluster{clust_ind}.prev_EMU_pair_num, cluster_EMU_size+1);
       output_file{end+1} = sprintf('right_mat_idv_d2 = zeros(%d, %d);\n', EMU.cluster{clust_ind}.prev_EMU_pair_num, cluster_EMU_size+1);
       
       if ~isempty(right_mat2{i})
           [aa, bb, cc] = find( EMU.cluster{clust_ind}.prev_EMU_pair_cluster_ind );
           for z = 1:length(aa)
               EMU1 = aa(z);
               EMU2 = bb(z);       
               EMU_pair_ind = cc(z);       
               m = CreateCauchyMat( EMU.size(EMU1)+1, EMU.size(EMU2)+1);
               [a,b,c] = find(m);
               
               m = sparse(a, b, idv_d{EMU1}(c));
               v1 = m*idv{EMU2}';

               m = sparse(a, b, idv{EMU1}(c));
               v2 = m*idv_d{EMU2}';
               v = v1+v2;

               right_mat_idv_d2(EMU_pair_ind, :) = v;
               
               for xx=1:(cluster_EMU_size+1)
                   s = '0';
                   for k=1:EMU.size(EMU1)+1
                        q = xx - k;
                        if (q>=0 & q < EMU.size(EMU2)+1)
                            if(idv_d{EMU1}(k) * idv{EMU2}(q+1))
                                s = sprintf('%s + idv_d{%d}(%d) * idv{%d}(%d)', s, EMU1, k, EMU2, q+1);
                            end
                            if(idv{EMU1}(k) * idv_d{EMU2}(q+1))
                                s = sprintf('%s + idv{%d}(%d) * idv_d{%d}(%d)', s, EMU1, k, EMU2, q+1);  
                            end
                        end
                   end
                   if(~strcmp(s,'0'))
                        output_file{end+1} = sprintf('right_mat_idv_d2(%d, %d) = %s;\n',EMU_pair_ind, xx, s); 
                   end
               end
               
           end
       end
       
 
       
       m1 = EMU.cluster{clust_ind}.left_mat_d{x} * left_mat_idv{i};
       m2 = EMU.cluster{clust_ind}.right_mat_d1{x}*right_mat_idv1{i} + right_mat1{i}*right_mat_idv_d1;
       m3 = EMU.cluster{clust_ind}.right_mat_d2{x}*right_mat_idv2{i} + right_mat2{i}*right_mat_idv_d2;       
       m4 = m2+m3-m1;
       t1 = tic;
       m = left_mat{i}\m4;       
       t2 = toc(t1);
       ticks = ticks + t2;
       

       if((sum(sum(abs(right_mat1{i}*right_mat_idv_d1))))&&(sum(sum(abs(right_mat2{i}*right_mat_idv_d2)))))
            output_file{end+1} = sprintf('m4 = right_mat1{%d}*right_mat_idv_d1 + right_mat2{%d}*right_mat_idv_d2;\n', i, i);
       elseif(sum(sum(abs(right_mat1{i}*right_mat_idv_d1))))
            output_file{end+1} = sprintf('m4 = right_mat1{%d}*right_mat_idv_d1;\n', i);
       else
            output_file{end+1} = sprintf('m4 = right_mat2{%d}*right_mat_idv_d2;\n', i);
       end
       
       
       if (sum(sum(abs(EMU.cluster{clust_ind}.right_mat_d1{x}))) > 0)
           [aa bb cc] = find(EMU.cluster{clust_ind}.right_mat_d1{x});
           for q=1:length(aa)
                output_file{end+1} = sprintf('m4(%d,:) = m4(%d,:) + %f .* right_mat_idv1{%d}(%d,:);\n', aa(q), aa(q), cc(q), i, bb(q));
           end
           %output_file{end+1} = sprintf('m4 = m4 + EMU.cluster{%d}.right_mat_d1{%d}*right_mat_idv1{%d};\n', clust_ind, x, i, i);
       end
       
       if (sum(sum(abs(EMU.cluster{clust_ind}.right_mat_d2{x}))) > 0)
           [aa bb cc] = find(EMU.cluster{clust_ind}.right_mat_d2{x});
           for q=1:length(aa)
                output_file{end+1} = sprintf('m4(%d,:) = m4(%d,:) + %f .* right_mat_idv2{%d}(%d,:);\n', aa(q), aa(q), cc(q), i, bb(q));
           end
%           output_file{end+1} = sprintf('m4 = m4 + EMU.cluster{%d}.right_mat_d2{%d}*right_mat_idv2{%d};\n', clust_ind, x, i, i);       
       end
       
       if (sum(sum(abs(EMU.cluster{clust_ind}.left_mat_d{x}))) > 0)
           [aa bb cc] = find(EMU.cluster{clust_ind}.left_mat_d{x});
           for q=1:length(aa)
                output_file{end+1} = sprintf('m4(%d,:) = m4(%d,:) - %f .* left_mat_idv{%d}(%d,:);\n', aa(q), aa(q), cc(q), i, bb(q));
           end
           
            %output_file{end+1} = sprintf('m4 = m4 - EMU.cluster{%d}.left_mat_d{%d} * left_mat_idv{i};\n', clust_ind, x);
       end
       %output_file{end+1} = sprintf('m = left_mat{%d}\\\\m4;\n', i);       
       output_file{end+1} = sprintf('m = left_mat_inv{%d}*m4;\n', i);       

       
       for z=1:EMU.cluster{clust_ind}.EMU_num
           emu_ind = EMU.cluster{clust_ind}.EMU_ind(z);
           if(~isequal(idv_d{emu_ind}, m(z,:)))
                idv_d{emu_ind} = m(z,:);
                output_file{end+1} = sprintf('idv_d{%d} = m(%d,:);\n', emu_ind, z);       
           end
           
           
           if result_idv_d_arr(emu_ind)
%                 result_idv_d{result_idv_d_arr(emu_ind)}(x, :) = m(z,:);
                if((~exist('result_idv_d{emu_ind}(x, :)')) || (~isequal(result_idv_d{emu_ind}(x, :), m(z,:))))
                    result_idv_d{emu_ind}(x, :) = m(z,:);
                    output_file{end+1} = sprintf('result_idv_d{%d}(%d, :) = m(%d,:);\n', emu_ind, x, z);       
                end
           end
       end

    end
   
end
%ticks
tic_end = toc(tic_begin);
fprintf('Fraction of time in inv = %f (%f)\n', ticks/tic_end, tic_end);

   f = fopen([fcn_name '.m'], 'wt');
   for i=1:length(output_file)
       fprintf(f, output_file{i});
   end
   fclose(f);

   a=1;
%toc