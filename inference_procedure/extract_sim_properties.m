function [E_sim,area_sim,tri_i,A_after,A_before] = ...
    extract_sim_properties(tri_sim_before,tri_sim_after,subset_sim_after,subset_sim_before,tip,ridge,beads,spacing,bias)

%%%% Calculate stretch
[E_sim,~,~,area_sim,A_after,A_before,~] = ...
    eigenvalue_calc_new(tri_sim_after,tri_sim_before,subset_sim_after,subset_sim_before,1);
IC_sim = incenter(tri_sim_before);
%%

%%%% Triangle error analysis
% displ = mean(abs(vecnorm(subset_sim_after - subset_sim_before)));
[tri_errorall,~] = tri_error(tri_sim_before,tri_sim_after,subset_sim_before,subset_sim_after,1,1);

%% 
%%%% Subset triangles (error and large area triangles at the back)
tri_i = find((tri_errorall <=mean(tri_errorall)+std(tri_errorall)) & IC_sim(:,1)'>0.3);

%%%% Circumferential and Meridional Stretch directions

% [mer_stretch,circ_stretch,mdirs,cdirs] = stretch_directions(tri_sim_before,subset_sim_before,subset_sim_after,[1;0;0],ridge);


%%%% Tension
% maxx = max(ridge(1,:));
% [curvs_tensions,princp_directions,xq1,yq1,zq1] = curvature_tension_calc_sim(ridge,tri_sim_before,subset_sim_before,IC_sim,spacing,bias,mdirs,cdirs);

% scatter3(xq1,yq1,zq1,5,'filled')





