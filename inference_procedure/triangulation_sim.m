function [subset_sim_before,subset_sim_after,tri_sim_before,tri_sim_after] = triangulation_sim(init_r,final_r,loc,Master_Connectivity)



subset_sim_after = init_r(:,loc);
subset_sim_before= final_r(:,loc);

% minz = min(subset_sim_before(3,:));
% tip = max(final_r(1,:));
% center = [0 0 minz-0.5];
% plane_pt = max(subset_sim_before(3,:)) +1;
% R = 1;
% 
% [~,tri_sim_before,scon_sim] = project_points_sim(subset_sim_before,center,plane_pt,tip,R);
% tri_sim_after = triangulation(scon_sim,subset_sim_after');

tri_sim_before = triangulation(Master_Connectivity,subset_sim_before');
tri_sim_after = triangulation(Master_Connectivity,subset_sim_after');
