function [final_points,tri_sim_before,triangle_areas_before] = rand_sim(final_r,init_r,N,dist_n)

randp = randi(size(final_r,2),1,N);
points = final_r(:,randp);
point_init = init_r(:,randp);

dists = squareform(pdist(points'));
dists_init = squareform(pdist(point_init'));

final_points = points;
del_all = [];
for i = 1:N
    del = find((dists(:,i)<=dist_n & dists(:,i)>0) | (dists_init(:,i)<=dist_n & dists_init(:,i)>0))';
    if isempty(del) == 1
    else
        i_g = del > i;

        del_all = [del_all del(i_g)];

    end
end
del_all = unique(del_all);
final_points(:,del_all) = [];
% 
final_points = final_points(:,final_points(3,:)>=-0.7);

minz = min(final_points(3,:));
tip = max(final_r(1,:));
center = [0 0 minz-0.2];
plane_pt = max(final_points(3,:))+1;
R = 1;

[~,tri_sim_before,~] = project_points_sim(final_points,center,plane_pt,tip,R);
% tri_sim_after = triangulation(scon_sim,subset_sim_after');

[~,~,~,~,~,triangle_areas_before,~] = ...
    eigenvalue_calc_new(tri_sim_before,tri_sim_before,final_points,final_points,0);
triangle_area = triangle_areas_before <= 1.2;
triangle_areas_before = triangle_areas_before(triangle_area);
tri_sim_before = triangulation(tri_sim_before.ConnectivityList(triangle_area,:),final_points');
% tri_sim_before.ConnectivityList(~triangle_area,:) = [];
% trisurf(tri_sim_before.ConnectivityList,final_points(1,:),final_points(2,:),final_points(3,:),'FaceAlpha',0.3);hold on;quickscatter(final_points)


