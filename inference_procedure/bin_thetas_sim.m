function [mean_var_angle,bins_mid_angles,angle_bins,mm_angle_all,angle_K,tri_K] = bin_thetas_sim(data,angles_T1,binsA,pn,sym,s,new,weighted,N_min,N_max)

% data(~posK_mu) = NaN;

if pn == 1
data = pre_process(data,1,size(data,2),s,new);
if new == 1
angles_T1 = angles_T1(s);
else
end
% data2 = pre_process(data2,1,size(data2,2),s2,new);
else
end
% angles_T2 = angles_T2(s2);

if sym == 1
angle_bins = linspace(N_min,N_max,(binsA*2)-1);

angle_bins(angle_bins ==10 | angle_bins == -10) = [];
else
% 

if binsA == 12
angle_bins = [0,20,30,40,50,60,70,75,80,82.5,85,87.5,90];
elseif binsA == 11
angle_bins = [0,20,30,40,50,60,70,80,82.5,85,87.5,90];
elseif binsA == 9
angle_bins = [0,20,30,40,50,60,70,80,85,90];
elseif binsA == 10
angle_bins = [0,20,30,40,50,60,70,75,80,85,90];
else
angle_bins = linspace(N_min,N_max,(N_max-N_min)/10+1);
end

angle_bins(angle_bins ==10) = [];
% angle_bins = [0,20,30,40,50,60,70,80,90,100,110];

end
% bins_mid_angles = [0,15,25,35,45,55,65,75,90];
% angle_bins = [-10,10,20,30,40,50,60,70,80,100];

angle_K = cell(1,size(angle_bins,2)-1);
tri_K = cell(1,size(angle_bins,2)-1);
weights = cell(1,size(angle_bins,2)-1);
bins_mid_angles = angle_bins(2:size(angle_bins,2))-(0.5*diff(angle_bins));


    
mm_angle_all = zeros(2,size(data,2));
% mm_angle_all2 = zeros(2,size(data2,2));


    for i = 1:size(angle_bins,2)-1
        for n=1:size(data,2)
             min_angle = min(angles_T1{n});
             max_angle = max(angles_T1{n});
             num_angles = size(find(angles_T1{n}<= angle_bins(i+1) & angles_T1{n} >= angle_bins(i)),1);
             if isempty(min_angle) == 1 || isempty(max_angle) == 1 
                 mm_angle_all(:,n) = [NaN;NaN];
             else
             mm_angle_all(:,n) = [min_angle;max_angle];
            if angle_bins(i) <= min_angle && angle_bins(i+1) >= min_angle
                    angle_K{i} = [angle_K{i} data(n)];
                    tri_K{i} = [tri_K{i} n];
                    
                    % weights{i} = [weights{i} angle_bins(i+1)-min_angle];
                    weights{i} = [weights{i} num_angles];
            elseif angle_bins(i) <= max_angle && angle_bins(i+1) >= max_angle
                    angle_K{i} = [angle_K{i} data(n)];
                    tri_K{i} = [tri_K{i} n];
                    weights{i} = [weights{i} num_angles];
            elseif angle_bins(i) >= min_angle && angle_bins(i+1) <= max_angle
                    angle_K{i} = [angle_K{i} data(n)]; 
                    tri_K{i} = [tri_K{i} n];
                    weights{i} = [weights{i} num_angles];
            end
             end
        end
        % for n=1:size(data2,2)
        %     min_angle2 = min(angles_T2{n});
        %     max_angle2 = max(angles_T2{n});
        %     mm_angle_all2(:,n) = [min_angle2;max_angle2];
        %     if angle_bins(i) <= min_angle2 && angle_bins(i+1) >= min_angle2
        %             angle_K{i} = [angle_K{i} data2(n)];
        %             weights{i} = [weights{i} angle_bins(i+1)-min_angle2];
        %     elseif angle_bins(i) <= max_angle2 && angle_bins(i+1) >= max_angle2
        %             angle_K{i} = [angle_K{i} data2(n)];
        %             weights{i} = [weights{i} max_angle2-angle_bins(i)];
        %     elseif angle_bins(i) >= min_angle2 && angle_bins(i+1) <= max_angle2
        %             angle_K{i} = [angle_K{i} data2(n)];  
        %             weights{i} = [weights{i} angle_bins(i+1)-angle_bins(i)];
        %     end
        % end
    end



mean_var_angle = zeros(2,size(angle_bins,2)-1);
for i = 1:size(angle_bins,2)-1
    if weighted == 1
    mean_var_angle(1,i) = sum(angle_K{i}(~isnan(angle_K{i})).*weights{i}(~isnan(angle_K{i})))/sum(weights{i}(~isnan(angle_K{i})));
    mean_var_angle(2,i) = std(angle_K{i},weights{i},'omitnan');
    else
    mean_var_angle(1,i) = mean(angle_K{i},'omitnan');
    mean_var_angle(2,i) = std(angle_K{i},'omitnan');
    end
    mean_var_angle(3,i) = size(angle_K{i}(~isnan(angle_K{i})),2);
    mean_var_angle(4,i) = sum(weights{i}./1000);
    if mean_var_angle(3,i) == 0
        mean_var_angle(1,i) = NaN;
    else
    end
end
