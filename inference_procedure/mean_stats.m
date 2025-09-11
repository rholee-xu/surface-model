function [mm_final] = mean_stats(table_means,weights,scale_all)

% weights = weights./scale_all;

mm_final = zeros(3,size(table_means,2));

mm_final(1,:) = sum(table_means.*weights,1,'omitnan')./sum(weights,1);
mm_final(2,:) = w_std(table_means,weights);
mm_final(3,:) = sum(~isnan(table_means),1);
mm_final(4,:) = sum(weights);

mm_final(mm_final(:) == 0) = NaN;

end