function [weighted_std] = w_std(data,weights)

weighted_std = zeros(1,size(data,2));
for n = 1:size(data,2)
weighted_std(:,n) = std(data(:,n),weights(:,n),'omitnan');
end