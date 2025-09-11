function [min_triangles,max_triangles] = max_min_t(tri_sim_before,subset_sim_before,dim)
max_triangles = zeros(1,size(tri_sim_before.ConnectivityList,2));
min_triangles = zeros(1,size(tri_sim_before.ConnectivityList,2));

for n = 1:size(tri_sim_before.ConnectivityList,1)
    con_n = tri_sim_before.ConnectivityList(n,:);
    points = subset_sim_before(dim,con_n);
    min_triangles(1,n) = min(points);
    max_triangles(1,n) = max(points);
end
% max_min_triangles = max_min_triangles(:,tri_i);
end