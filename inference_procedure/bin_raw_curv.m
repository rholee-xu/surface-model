function [mean_var_angle,bins_mid_angles,angle_bins,mm_angle_all,angle_K] = bin_raw_curv(data,angle_data,N_min,N_max,N)

% scatter(angle_data,data);
% [angle_data,s] = sort(angle_data);
% data = pre_process(data,1,size(data,2),s,1);
% hold on; scatter(angle_data,data)
angle_bins = linspace(N_min,N_max,N);
bins_mid_angles = angle_bins(2:size(angle_bins,2))-(0.5*diff(angle_bins));


angle_K = cell(1,size(angle_bins,2)-1);

% data(~(data >= mean(data,'omitnan')-2*std(data,'omitnan') & data <= mean(data,'omitnan')+2*std(data,'omitnan'))) = NaN;
    
mm_angle_all = zeros(2,size(data,2));
% mm_angle_all2 = zeros(2,size(data2,2));


    for i = 1:size(angle_bins,2)-1
        for n=1:size(data,2)
            if angle_bins(i) <= angle_data(n) && angle_bins(i+1) >= angle_data(n)
                angle_K{i} = [angle_K{i} data(n)];
            else

            end
         end
     end

for i = 1:size(angle_K,2)
    if isempty(angle_K{i})
    else
    angle_K{i} = pre_process(angle_K{i},1,size(angle_K{i},2),0,0);
    end

end



mean_var_angle = zeros(2,size(angle_bins,2)-1);
for i = 1:size(angle_bins,2)-1
    mean_var_angle(1,i) = mean(angle_K{i},'omitnan');
    mean_var_angle(2,i) = std(angle_K{i},'omitnan');
    mean_var_angle(3,i) = size(angle_K{i}(~isnan(angle_K{i})),2);
    if mean_var_angle(3,i) == 0
        mean_var_angle(1,i) = NaN;
    else
    end
end
end