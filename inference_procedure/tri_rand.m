function [Master_con,Master_con2,Master_tri_i,Master_tri_i2,locpb,locpb2,Master_i_tri,Master_i_tri2] = tri_rand(final_r,init_r,n_init,dist_n)

z_end = round(max(init_r(1,:)));
% generate triangulation and set of beads
[final_points,tri_before] = rand_sim(final_r,init_r,n_init*(z_end/2),dist_n);
[final_points2,tri_before2] = rand_sim(final_r,init_r,n_init*(z_end/2),dist_n);


%%
loc_final = zeros(1,size(final_points,2));
for n = 1:size(final_points,2)
    p = final_points(:,n);
    c = ismember(final_r,p);
   
    loc_final(n) = find(c(1,:) == 1 & c(2,:) == 1,1);% & c(3,:) == 1,1);
end

loc_final2 = zeros(1,size(final_points2,2));
for n = 1:size(final_points2,2)
    p = final_points2(:,n);
    c = ismember(final_r,p);
   
    loc_final2(n) = find(c(1,:) == 1 & c(2,:) == 1,1);% & c(3,:) == 1);
end

final_points_after = init_r(:,loc_final);
tri_after = triangulation(tri_before.ConnectivityList,final_points_after');
final_points_after2 = init_r(:,loc_final2);
tri_after2 = triangulation(tri_before2.ConnectivityList,final_points_after2');
%%
[~,~,~,~,~,A_before,~] = ...
    eigenvalue_calc_new(tri_after,tri_before,final_points_after,final_points,0);
[~,~,~,~,~,A_before2,~] = ...
    eigenvalue_calc_new(tri_after2,tri_before2,final_points_after2,final_points2,0);
%%

% get master tri_i 
IC_temp = incenter(tri_before);
[~,~,~,~,~,~,tri_errorall] = ...
    eigenvalue_calc_new(tri_after,tri_before,final_points_after,final_points,1);
% [tri_errorall,del_errorall] = tri_error(tri_before,tri_after,final_points,final_points_after,1,1);
i_tri = (log(tri_errorall) <=mean(log(tri_errorall))+std(log(tri_errorall))) & IC_temp(:,1)'>0.3 & A_before < 1;
tri_i = find((log(tri_errorall) <=mean(log(tri_errorall))+std(log(tri_errorall))) & IC_temp(:,1)'>0.3 & A_before < 1);
% trisurf(tri_before.ConnectivityList,final_points(1,:),final_points(2,:),final_points(3,:),'FaceAlpha',0);axis equal
% 

%%
IC_temp2 = incenter(tri_before2);
[~,~,~,~,~,~,tri_errorall2] = ...
    eigenvalue_calc_new(tri_after2,tri_before2,final_points_after2,final_points2,1);
% [tri_errorall2,~] = tri_error(tri_before2,tri_after2,final_points2,final_points_after2,1,1);
i_tri2 = (log(tri_errorall2) <=mean(log(tri_errorall2))+std(log(tri_errorall2))) & IC_temp2(:,1)'>0.3 & A_before2 < 1;
tri_i2 = find((log(tri_errorall2) <=mean(log(tri_errorall2))+std(log(tri_errorall2))) & IC_temp2(:,1)'>0.3 & A_before2 < 1);
% trisurf(tri_before2.ConnectivityList,final_points2(1,:),final_points2(2,:),final_points2(3,:),'FaceAlpha',0);axis equal
% 

% tiledlayout(1,2)
% nexttile;trisurf(tri_before.ConnectivityList,final_points(1,:),final_points(2,:),final_points(3,:),'FaceAlpha',0);hold on;quickscatter(final_points)
% hold on;trisurf(tri_before.ConnectivityList(tri_i,:),final_points(1,:),final_points(2,:),final_points(3,:),'FaceAlpha',0.3);axis equal
% nexttile;trisurf(tri_before2.ConnectivityList,final_points2(1,:),final_points2(2,:),final_points2(3,:),'FaceAlpha',0);hold on;quickscatter(final_points2)
% hold on;trisurf(tri_before2.ConnectivityList(tri_i2,:),final_points2(1,:),final_points2(2,:),final_points2(3,:),'FaceAlpha',0.3);axis equal
%%

% Set masters

Master_con = tri_before.ConnectivityList;
Master_tri_i = tri_i;
locpb = loc_final;
Master_con2 = tri_before2.ConnectivityList;
Master_tri_i2 = tri_i2;
locpb2 = loc_final2;
Master_i_tri = i_tri;
Master_i_tri2 = i_tri2;

end