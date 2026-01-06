%% Load Raw Data

% Load the csv output from Ridge Detection with inputs for the lowest and
% highest z-stack value
before_xy_ridge_all = get_ridge_new('ridge_before.csv',15,90); 
after_xy_ridge_all = get_ridge_new('ridge_after.csv',18,90);
% The function will output the point cloud data in 3xn form

% Load the csv output from RS-FISH
raw_beads_before = readtable('beads_before.csv');
beads_before = get_beads(raw_beads_before);
raw_beads_after = readtable('beads_after.csv');
beads_after = get_beads(raw_beads_after);
% The function will output the point cloud data in 3xn form
%% Adjust beads for matching
 quickscatter(beads_before); hold on; quickscatter(beads_after);
 % beads_after(1,:) = beads_after(1,:) - 10;
 % beads_after(2,:) = beads_after(2,:) -10;
 % beads_after(3,:) = beads_after(3,:) +10;
%% Subset beads for better visualization
beads_subset_after = beads_after(:,beads_after(2,:) <70);
beads_subset_before = beads_before(:,beads_before(2,:) <70);
quickscatter(beads_subset_before); hold on; quickscatter(beads_subset_after);
%%
subset = zeros(3,size(cursor_info,2));
for n = 1:size(cursor_info,2)
    subset(:,n) = cursor_info(n).Position;
end

subset_all_after = subset;
subset_all_after2 = [subset_all_after2 subset];

%%
quickscatter(beads_subset_before); hold on; quickscatter(beads_subset_after);
quickscatter(subset_all_before2); %hold on; quickscatter(subset_all_after);

%% bring beads back

% subset_all_after(1,:) = subset_all_after(1,:) - 20;
% subset_all_after(2,:) = subset_all_after(2,:) - 10;
% subset_all_after(3,:) = subset_all_after(3,:) - 10;
%% Sort beads before checking their order
[~,order1] =  sort(subset_all2_after(1,:));
subset_all2_after = subset_all2_after(:,order1);
[~,order2] =  sort(subset_all2_before(1,:));
subset_all2_before = subset_all2_before(:,order2);
[~,order3] =  sort(subset_all_after(1,:));
subset_all_after = subset_all_after(:,order3);
[~,order4] =  sort(subset_all_before(1,:));
subset_all_before = subset_all_before(:,order4);

%%
check_beads_order(subset_all2_before,subset_all2_after)

%% Calculate tip and side radius
subset_tip_sphere = before_xy_ridge_all(:,before_xy_ridge_all(1,:)>=190);
sphere_obj = sphericalFit(subset_tip_sphere);
rad = sphere_obj.radius;

subset_tip_sphere_after = after_xy_ridge_all(:,after_xy_ridge_all(1,:)>=180);
sphere_obj_after = sphericalFit(subset_tip_sphere_after);
rad_after = sphere_obj_after.radius;
%%
tiledlayout(1,2)
nexttile;showfit(sphere_obj,'FaceAlpha',0.2,'FaceColor','g');
hold on; quickplot(subset_tip_sphere)
nexttile;showfit(sphere_obj_after,'FaceAlpha',0.2,'FaceColor','g');
hold on; quickplot(subset_tip_sphere_after);
%%
tip = max(before_xy_ridge_all(1,:));
side_ridge = before_xy_ridge_all(:,before_xy_ridge_all(1,:)<=150 & before_xy_ridge_all(1,:)>=0);
cyl = cylindricalFit(side_ridge);
rad_side = cyl.radius
showfit(cyl,'FaceAlpha',0.2,'FaceColor','g');hold on; quickplot(side_ridge)
%% Automated Triangulation
tiledlayout(1,2);nexttile
[~,tri_all_before,~]=project_points(subset_all_before,before_xy_ridge_all,rad,3);
tri_all_after = triangulation(tri_all_before.ConnectivityList,subset_all_after');
hold on; quickscatter(subset_all_before);
nexttile;
[~,tri_all2_before,~]=project_points(subset_all2_before,before_xy_ridge_all,rad,3);
tri_all2_after = triangulation(tri_all2_before.ConnectivityList,subset_all2_after');
hold on; quickscatter(subset_all2_before);

IC_all = incenter(tri_all_before);
IC_all2 = incenter(tri_all2_before);
%%

% tiledlayout(1,2);nexttile
% trisurf(tri_all_before.ConnectivityList,subset_all_before(1,:),subset_all_before(2,:),subset_all_before(3,:),...
%     'FaceAlpha',0);axis equal
% text(IC_all(:,1),IC_all(:,2),IC_all(:,3),string(1:size(tri_all_before.ConnectivityList,1)))
% nexttile
trisurf(tri_all2_before.ConnectivityList,subset_all2_before(1,:),subset_all2_before(2,:),subset_all2_before(3,:),...
    'FaceAlpha',0.3);axis equal
text(IC_all2(:,1),IC_all2(:,2),IC_all2(:,3),string(1:size(tri_all2_before.ConnectivityList,1)))
%%
delete = tri_all_before.ConnectivityList;
delete([20,21,18,19,11,14,2,10,28,5],:) = [];
delete = [delete;20,24,28;20,21,28];

tri_all_before_new = triangulation(delete,subset_all_before_new');
tri_all_after_new = triangulation(delete,subset_all_after_new');
delete2 = tri_all2_before_new.ConnectivityList;
delete2([22,4,20,1,18,8],:) = [];
tri_all2_before = triangulation(delete2,subset_all2_before');
tri_all2_after = triangulation(delete2,subset_all2_after');
f10 = figure(10);clf()
f10.Position(3) = 1000;
tiledlayout(1,2);nexttile
trisurf(tri_all_before.ConnectivityList,subset_all_before(1,:),subset_all_before(2,:),subset_all_before(3,:),...
    'FaceAlpha',0.3);axis equal
nexttile;trisurf(tri_all2_before.ConnectivityList,subset_all2_before(1,:),subset_all2_before(2,:),subset_all2_before(3,:),...
    'FaceAlpha',0.3);axis equal
set(gcf,'Position',[100 100 1000 500])
%%
tiledlayout(1,2);nexttile
trisurf(tri_all_after.ConnectivityList,subset_all_after(1,:),subset_all_after(2,:),subset_all_after(3,:),...
    'FaceAlpha',0.3);axis equal
nexttile;trisurf(tri_all2_after.ConnectivityList,subset_all2_after(1,:),subset_all2_after(2,:),subset_all2_after(3,:),...
    'FaceAlpha',0.3);axis equal
set(gcf,'Position',[100 100 1000 500])
%% bead projection
[subset_ob_before,subset_ob_after,tri_ob_before,tri_ob_after,iz] = ...
    bead_projection(subset_all_before,subset_all_after,before_xy_ridge_all,after_xy_ridge_all,tri_all_before);
[subset_ob2_before,subset_ob2_after,tri_ob2_before,tri_ob2_after,iz2] = ...
    bead_projection(subset_all2_before,subset_all2_after,before_xy_ridge_all,after_xy_ridge_all,tri_all2_before);

tiledlayout(1,2);nexttile

quickscatter(subset_ob_before(:,iz));hold on; quickscatter(subset_all_before);
hold on; trisurf(tri_ob_before.ConnectivityList,subset_ob_before(1,:),subset_ob_before(2,:),subset_ob_before(3,:),'FaceAlpha',0.3);axis equal
hold on; quickplot(before_xy_ridge_all)
nexttile;quickscatter(subset_ob2_before(:,iz2));hold on; quickscatter(subset_all2_before);
hold on; trisurf(tri_ob2_before.ConnectivityList,subset_ob2_before(1,:),subset_ob2_before(2,:),subset_ob2_before(3,:),'FaceAlpha',0.3);axis equal
set(gcf,'Position',[100 100 1000 500])
%%
tiledlayout(1,2);nexttile

quickscatter(subset_ob_after(:,iz));hold on; quickscatter(subset_all_after_new);
hold on; trisurf(tri_ob_after.ConnectivityList,subset_ob_after(1,:),subset_ob_after(2,:),subset_ob_after(3,:),'FaceAlpha',0.3);axis equal
hold on; quickplot(after_xy_ridge_all)
nexttile;quickscatter(subset_ob2_after(:,iz2));hold on; quickscatter(subset_all2_after_new);
hold on; trisurf(tri_ob2_after.ConnectivityList,subset_ob2_after(1,:),subset_ob2_after(2,:),subset_ob2_after(3,:),'FaceAlpha',0.3);axis equal
hold on; quickplot(after_xy_ridge_all)
set(gcf,'Position',[100 100 1000 500])
%% redo before axis

xz_angle = -0.95;xy_angle = -8.71; 
new_long_axis = get_long_axis(xz_angle,xy_angle);

rot1 = rotateAbtLine(before_xy_ridge_all,-deg2rad(xz_angle),2);
rotated_ridge_before = rotateAbtLine(rot1,deg2rad(xy_angle),3);

xz_angle_af = 3.33;xy_angle_af = -10.51; 
[~,i] = maxk(rotated_ridge_before(1,:),1);
middle = rotated_ridge_before(3,i);
rotated_ridge_flat = rotated_ridge_before(1:2,rotated_ridge_before(3,:)<=middle+1 & rotated_ridge_before(3,:)>=middle-1);


new_long_axis_After = get_long_axis(xz_angle_af,xy_angle_af);
rot2 = rotateAbtLine(after_xy_ridge_all,-deg2rad(xz_angle_af),2);
rotated_ridge_after = rotateAbtLine(rot2,deg2rad(xy_angle_af),3);
[~,i] = maxk(rotated_ridge_after(1,:),1);
middle = rotated_ridge_after(3,i);
rotated_ridge_flat_after = rotated_ridge_after(1:2,rotated_ridge_after(3,:)<=middle+1 & rotated_ridge_after(3,:)>=middle-1);

rotated_ridge_flat(1,:) = rotated_ridge_flat(1,:)-max(rotated_ridge_flat(1,:));
rotated_ridge_flat(2,:) = rotated_ridge_flat(2,:)-max(rotated_ridge_flat(2,:));
rotated_ridge_flat_after(1,:) = rotated_ridge_flat_after(1,:)-max(rotated_ridge_flat_after(1,:));
rotated_ridge_flat_after(2,:) = rotated_ridge_flat_after(2,:)-max(rotated_ridge_flat_after(2,:));


plot(rotated_ridge_flat(1,:),rotated_ridge_flat(2,:),'.');axis equal 
hold on;plot(rotated_ridge_flat_after(1,:),rotated_ridge_flat_after(2,:),'.');axis equal 

%%
max_del = 1;
[tri_errorob,~] = tri_error(tri_ob_before,tri_ob_after,subset_ob_before,subset_ob_after,max_del,max_del);
i_ob = (log(tri_errorob) <=mean(log(tri_errorob))+std(log(tri_errorob)));
not_tri = find(~i_ob);
tri_i_ob = find(i_ob);
[min_ob,max_ob]=max_min_t(tri_ob_before,subset_ob_before,1);
IC_ob = incenter(tri_ob_before);

[tri_errorob2,~] = tri_error(tri_ob2_before,tri_ob2_after,subset_ob2_before,subset_ob2_after,max_del,max_del);
i_ob2 = (log(tri_errorob2) <=mean(log(tri_errorob2))+std(log(tri_errorob2)));
not_tri2 = find(~i_ob2);
tri_i_ob2 = find(i_ob2);
[min_ob2,max_ob2]=max_min_t(tri_ob2_before,subset_ob2_before,1);
IC_ob2 = incenter(tri_ob2_before);


[mer_stretch_ob_test,circ_stretch_ob_test,mdirs_ob_test,cdirs_ob_test,pdir1_ob,pdir2_ob] = stretch_directions(tri_ob_before,subset_ob_before,subset_ob_after,new_long_axis',before_xy_ridge_all);
[mer_stretch_ob2_test,circ_stretch_ob2_test,mdirs_ob2_test,cdirs_ob2_test,pdir1_ob2,pdir2_ob2] = stretch_directions(tri_ob2_before,subset_ob2_before,subset_ob2_after,new_long_axis',before_xy_ridge_all);
[curvs_tensions_ob3_test,princp_directions_ob3,ref_km,ref_kc,ref_tm,ref_tc,ref_K,ref_MU,angles_sub_all] = ...
    curvature_tension_calc_exp(before_xy_ridge_all,tri_ob_before,subset_ob_before,IC_ob,tip,rad,3,new_long_axis',cdirs_ob_test,mdirs_ob_test,circ_stretch_ob_test,mer_stretch_ob_test,tri_i_ob);
[curvs_tensions_ob23_test,princp_directions_ob23,ref_km2,ref_kc2,ref_tm2,ref_tc2,ref_K2,ref_MU2,angles_sub_all2,norm_comp_all2,xq_all2,yq_all2,zq_all2] = ...
    curvature_tension_calc_exp(before_xy_ridge_all,tri_ob2_before,subset_ob2_before,IC_ob2,tip,rad,3,new_long_axis',cdirs_ob2_test,mdirs_ob2_test,circ_stretch_ob2_test,mer_stretch_ob2_test,tri_i_ob2);
[~,~,~,~,~,A_ob_before,~] = ...
    eigenvalue_calc_new(tri_ob_after,tri_ob_before,subset_ob_after,subset_ob_before,1);
[~,~,~,~,~,A_ob2_before,~] = ...
    eigenvalue_calc_new(tri_ob2_after,tri_ob2_before,subset_ob2_after,subset_ob2_before,1);


%%


[K_ob_test, mu_ob_test] = compute_K_mu_indiv(abs(curvs_tensions_ob3_test),[mer_stretch_ob_test;mer_stretch_ob_test;circ_stretch_ob_test]);
[K_ob2_test, mu_ob2_test] = compute_K_mu_indiv(abs(curvs_tensions_ob23_test),[mer_stretch_ob2_test;mer_stretch_ob2_test;circ_stretch_ob2_test]);


new_i_test = find(K_ob_test'>0 & i_ob');
new_i2_test = find(K_ob2_test'>0 & i_ob2');

tip_a = max(after_xy_ridge_all(1,:));

%% Angle calculation and binning

[angles_T1_test,~,mm_angle1_test,deg_ridge,~,deg_all,x_all,y_all,z_all] = get_thetas(tri_ob_after,after_xy_ridge_all,subset_ob_after,rad,tip_a,0,[1;0;0],new_long_axis_After');
[angles_T2_test,~,mm_angle2_test,~,~,deg_all2,x_all2,y_all2,z_all2] = get_thetas(tri_ob2_after,after_xy_ridge_all,subset_ob2_after,rad,tip_a,0,[1;0;0],new_long_axis_After');

scatter3(x_all,y_all,z_all,10,deg_all,'filled');axis equal;colorbar;
climm = [0,100];cmap = parula(10);colormap(cmap);clim(climm)
%%
angles_T1_sub_test = angles_T1_test(new_i_test);
angles_T2_sub_test = angles_T2_test(new_i2_test);
mm_angle1_sub_test = mm_angle1_test(:,new_i_test);mm_angle2_sub_test = mm_angle2_test(:,new_i2_test);
mm_angle_all_test = [mm_angle1_sub_test,mm_angle2_sub_test];
[~,s_all] = sort([mean(mm_angle_all_test)]);

binsA = 12;
[ww_angles_circ_stretch_t_new2,bins_angles_mid,angle_bins,~,ww_cell_circ_stretch] = bin_thetas([circ_stretch_ob_test(new_i_test),circ_stretch_ob2_test(new_i2_test)],[angles_T1_sub_test,angles_T2_sub_test],binsA,1,0,s_all,1,1);
[ww_angles_mer_stretch_t_new2,~,~,~,ww_cell_mer_stretch] = bin_thetas([mer_stretch_ob_test(new_i_test),mer_stretch_ob2_test(new_i2_test)],[angles_T1_sub_test,angles_T2_sub_test],binsA,1,0,s_all,1,1);

[ww_angles_circ_curv_t_new2,~,~,~] = bin_thetas([abs(curvs_tensions_ob3_test(3,new_i_test)),abs(curvs_tensions_ob23_test(3,new_i2_test))],[angles_T1_sub_test,angles_T2_sub_test],binsA,1,0,s_all,1,1);
[ww_angles_mer_curv_t_new2,~,~,~,angle_curv,tri_curv] = bin_thetas([abs(curvs_tensions_ob3_test(4,new_i_test)),abs(curvs_tensions_ob23_test(4,new_i2_test))],[angles_T1_sub_test,angles_T2_sub_test],binsA,1,0,s_all,1,1);

[ww_angles_curv1_t_new2,~,~,~] = bin_thetas([abs(curvs_tensions_ob3_test(1,new_i_test)),abs(curvs_tensions_ob23_test(1,new_i2_test))],[angles_T1_sub_test,angles_T2_sub_test],binsA,1,0,s_all,1,1);
[ww_angles_curv2_t_new2,~,~,~] = bin_thetas([abs(curvs_tensions_ob3_test(2,new_i_test)),abs(curvs_tensions_ob23_test(2,new_i2_test))],[angles_T1_sub_test,angles_T2_sub_test],binsA,1,0,s_all,1,1);
[ww_angles_tens1_t_new2,~,~,~] = bin_thetas([abs(curvs_tensions_ob3_test(6,new_i_test)),abs(curvs_tensions_ob23_test(6,new_i2_test))],[angles_T1_sub_test,angles_T2_sub_test],binsA,1,0,s_all,1,1);
[ww_angles_tens2_t_new2,~,~,~,angle_curv,tri_curv] = bin_thetas([abs(curvs_tensions_ob3_test(5,new_i_test)),abs(curvs_tensions_ob23_test(5,new_i2_test))],[angles_T1_sub_test,angles_T2_sub_test],binsA,1,0,s_all,1,1);

[ww_angles_circ_tens_t_new2,~,~,~,ww_cell_circ_tens] = bin_thetas([abs(curvs_tensions_ob3_test(8,new_i_test)),abs(curvs_tensions_ob23_test(8,new_i2_test))],[angles_T1_sub_test,angles_T2_sub_test],binsA,1,0,s_all,1,1);
[ww_angles_mer_tens_t_new2,~,~,~,ww_cell_mer_tens] = bin_thetas([abs(curvs_tensions_ob3_test(7,new_i_test)),abs(curvs_tensions_ob23_test(7,new_i2_test))],[angles_T1_sub_test,angles_T2_sub_test],binsA,1,0,s_all,1,1);

[ww_angles_K_t_new2,~,~,~,ww_cell_K,ww_tri] = bin_thetas([K_ob_test(new_i_test),K_ob2_test(new_i2_test)],[angles_T1_sub_test,angles_T2_sub_test],binsA,1,0,s_all,1,1);
[ww_angles_mu_t_new2,~,~,~,ww_cell_mu] = bin_thetas([mu_ob_test(new_i_test),mu_ob2_test(new_i2_test)],[angles_T1_sub_test,angles_T2_sub_test],binsA,1,0,s_all,1,1);

%%
tiledlayout(1,4)
nexttile;plot_bins_error(angle_bins,bins_angles_mid,ww_angles_circ_stretch_t_new2,ww_angles_mer_stretch_t_new2);ylim([0.95,1.3]);title('Stretch');xlabel('Local Angle')
nexttile;plot_bins_error(angle_bins,bins_angles_mid,ww_angles_circ_curv_t_new2,ww_angles_mer_curv_t_new2);ylim([-0.05,0.3]);title('Curvature');xlabel('Local Angle')
nexttile;plot_bins_error(angle_bins,bins_angles_mid,ww_angles_circ_tens_t_new2,ww_angles_mer_tens_t_new2);ylim([-5,35]);title('Tension');xlabel('Local Angle')
nexttile;plot_bins_error_single(angle_bins,bins_angles_mid,ww_angles_K_t_new2,'m');ylim([0,150]);title('K');xlabel('Local Angle')
set(gcf,'Position',[100 100  1098 318])
