

% Use main_material_function.m to generate structs with cell information -
% turgid/unturgid outlines, K/MU dsitributions, etc.

% Load the structs you want to run the inference on 
% ex. load('Run_difgradients_database','Run_CANew_Con10','Run_CANew_Con8','Run_CANew_Con7','Run_CANew_Con9')

% Code is setup to run 4 unique cell cases at the three triangulation levels
% triangulation size, noise, and trial number can be adjusted
%% small
warning('off','all')

long_axis = [1;0;0];
sym = 0;
N_min = 0;
N_max = 90;
% noise parameters as standard deviations
magni = 0.002;
magni_init = 0.012;
magni_z = 0.03;
trialNum = 100; 

% triangulation parameters
init_n = 130;
dist_n = 0.05;

% binning
binsA = 10;

[structAll_Con10CAXS_nNormL_bin10_Kfilt,structAll_Con10CAXS_nNormL_bin10_nofilt]...
= simulation_pol_mv(Run_CANew_Con10,long_axis,sym,N_min,N_max,magni,trialNum,magni_init,init_n,dist_n,magni_z,binsA);
[structAll_Con7CAXS_nNormL_bin10_Kfilt,structAll_Con7CAXS_nNormL_bin10_nofilt]...
= simulation_pol_mv(Run_CANew_Con7,long_axis,sym,N_min,N_max,magni,trialNum,magni_init,init_n,dist_n,magni_z,binsA);
[structAll_Con8CAXS_nNormL_bin10_Kfilt,structAll_Con8CAXS_nNormL_bin10_nofilt]...
= simulation_pol_mv(Run_CANew_Con8,long_axis,sym,N_min,N_max,magni,trialNum,magni_init,init_n,dist_n,magni_z,binsA);
[structAll_Con9CAXS_nNormL_bin10_Kfilt,structAll_Con9CAXS_nNormL_bin10_nofilt]...
= simulation_pol_mv(Run_CANew_Con9,long_axis,sym,N_min,N_max,magni,trialNum,magni_init,init_n,dist_n,magni_z,binsA);

%% medium

init_n = 80;
dist_n = 0.25;

[structAll_Con10CA_nNormL_bin10_Kfilt,structAll_Con10CA_nNormL_bin10_nofilt]= simulation_pol_mv(Run_CANew_Con10,long_axis,sym,N_min,N_max,magni,trialNum,magni_init,init_n,dist_n,magni_z,binsA);
[structAll_Con7CA_nNormL_bin10_Kfilt,structAll_Con7CA_nNormL_bin10_nofilt]...
= simulation_pol_mv(Run_CANew_Con7,long_axis,sym,N_min,N_max,magni,trialNum,magni_init,init_n,dist_n,magni_z,binsA);
[structAll_Con8CA_nNormL_bin10_Kfilt,structAll_Con8CA_nNormL_bin10_nofilt]...
= simulation_pol_mv(Run_CANew_Con8,long_axis,sym,N_min,N_max,magni,trialNum,magni_init,init_n,dist_n,magni_z,binsA);
[structAll_Con9CA_nNormL_bin10_Kfilt,structAll_Con9CA_nNormL_bin10_nofilt]...
= simulation_pol_mv(Run_CANew_Con9,long_axis,sym,N_min,N_max,magni,trialNum,magni_init,init_n,dist_n,magni_z,binsA);

%% large

init_n = 70;
dist_n = 0.3;

[structAll_Con10CAL_nNormL_bin10_Kfilt,structAll_Con10CAL_nNormL_bin10_nofilt]= simulation_pol_mv(Run_CANew_Con10,long_axis,sym,N_min,N_max,magni,trialNum,magni_init,init_n,dist_n,magni_z,binsA);
[structAll_Con7CAL_nNormL_bin10_Kfilt,structAll_Con7CAL_nNormL_bin10_nofilt]...
= simulation_pol_mv(Run_CANew_Con7,long_axis,sym,N_min,N_max,magni,trialNum,magni_init,init_n,dist_n,magni_z,binsA);
[structAll_Con8CAL_nNormL_bin10_Kfilt,structAll_Con8CAL_nNormL_bin10_nofilt]...
= simulation_pol_mv(Run_CANew_Con8,long_axis,sym,N_min,N_max,magni,trialNum,magni_init,init_n,dist_n,magni_z,binsA);
[structAll_Con9CAL_nNormL_bin10_Kfilt,structAll_Con9CAL_nNormL_bin10_nofilt]...
= simulation_pol_mv(Run_CANew_Con9,long_axis,sym,N_min,N_max,magni,trialNum,magni_init,init_n,dist_n,magni_z,binsA);


save('trial100_AllCA_nNormL_NEW','-append')
%%
function [structAll_Kfilt,structAll_nofilt] = ...
    simulation_pol_mv(All_runs,long_axis,sym,N_min,N_max,magni,trialNum,magni_init,init_n,dist_n,magni_z,binsA)
% Set trials and perturbation
% trialNum = 1;
spacing = 0.05;

% saving information (generally not used)
infTensMs=cell(1,trialNum); % tension m
infTensCs=cell(1,trialNum); % tension c
infCurvMs=cell(1,trialNum); % curvature m
infCurvCs=cell(1,trialNum); % curvature c
infStrainMs=cell(1,trialNum); % strain m
infStrainCs=cell(1,trialNum); % strain c
infKs_124=cell(1,trialNum); % K
infMUs_124=cell(1,trialNum); % mu
mm_angles_all=cell(1,trialNum);
infStrainArea = cell(1,trialNum);
angles_T_all=cell(1,trialNum);

K128_trials = zeros(trialNum,128);
MU128_trials = zeros(trialNum,128);
synCurvM128_trials = zeros(trialNum,128);
synCurvC128_trials = zeros(trialNum,128);
synTensM128_trials = zeros(trialNum,128);
synTensC128_trials = zeros(trialNum,128);
synStrainM128_trials = zeros(trialNum,128);
synStrainC128_trials = zeros(trialNum,128);
rv_final128_1_trials = zeros(trialNum,129);
rv_final128_2_trials = zeros(trialNum,129);
scale_all = zeros(trialNum,1);
% Ws=zeros(trialNum,N_T);

for trialCount=1:trialNum
disp(trialCount)

long_axis = [1;0;0];
N_min = 0;
N_max = 90;
    % generating cell of different width (within a range) each trial number
    Run = All_runs;
    scale = (1+0.4*rand(1,1));
    rv_final128 = Run.rv_final128.*scale;
    rv_init128 = Run.rv_init128.*scale;
    rv_init = Run.rv_init.*scale;
    rv_final = Run.rv_final.*scale;


    K128 = Run.K128;
    MU128 = Run.MU128;
    synCurvM128 = Run.synCurvM128;
    synCurvC128 = Run.synCurvC128;
    synTensM128 = Run.synTensM128;
    synTensC128 = Run.synTensC128;
    synStrainM128 = Run.synStrainM128;
    synStrainC128 = Run.synStrainC128;

    %     rv_final128,rv_init128,rv_init,rv_final
    initpb_r = [];
    finalpb_r = [];
    initpb128_r = [];
    finalpb128_r = [];

    % rotation of synthetic cell into 3D marker points and outline
    for n = -pi:0.1:pi
	
        % adding noise
        [rvpb_init128] = perturb_new_norm(rv_init128,magni_init,magni_z);  
        [rvpb_final128] = perturb_new_norm(rv_final128,magni,0.005);
        %   and then extract the perturbed markers again
        [rvpb_init]=extract_markers(rvpb_init128);
        [rvpb_final]=extract_markers(rvpb_final128);
    
        rv_finalpb_3D = rvpb_final;%zeros(1,size(rvpb_final,2))];
        rv_initpb_3D = rvpb_init;%zeros(1,size(rvpb_init,2))];
     
        finalpb_beads = rotateAbtLine(rv_finalpb_3D,n,1);
        initpb_beads = rotateAbtLine(rv_initpb_3D,n,1);
        
        finalpb_r = [finalpb_r finalpb_beads];
        initpb_r = [initpb_r initpb_beads];
        
    end
    for n = -pi:0.05:pi
        rv_finalpb1210_6D = [rv_final128;zeros(1,size(rv_final128,2))];
        rv_initpb1210_6D = [rv_init128;zeros(1,size(rv_init128,2))];
        initpb128_beads = rotateAbtLine(rv_initpb1210_6D,n,1);
        finalpb128_beads = rotateAbtLine(rv_finalpb1210_6D,n,1);
        initpb128_r = [initpb128_r initpb128_beads];
        finalpb128_r = [finalpb128_r finalpb128_beads];
    end
        % rand_a = 0.05*rand(1,2);
        % ry1 = rotateAbtLine(finalpb_r,rand_a(1),2);
        % finalpb_r = rotateAbtLine(ry1,rand_a(2),3);
        % ry2 = rotateAbtLine(finalpb128_r,rand_a(1),2);
        % finalpb128_r = rotateAbtLine(ry2,rand_a(2),3);
        % ry3 = rotateAbtLine(initpb_r,rand_a(1),2);
        % initpb_r = rotateAbtLine(ry3,rand_a(2),3);

% save curvature info on outline for both triangulations
[structCurv_info] = get_curv_raw_sim(finalpb128_r,spacing,magni);
              [Master_con,Master_con1,Master_tri_i,Master_tri_i2,locpb,locpb2,Master_i_tri,Master_i_tri2] = tri_rand(finalpb_r,initpb_r,init_n,dist_n);


    % get triangulation and marker point information
    [subset_simpb_before,subset_simpb_after,tri_simpb_before,tri_simpb_after] = triangulation_sim(initpb_r,finalpb_r,locpb,Master_con);
    tip = max(finalpb_r(1,:));
    IC = incenter(tri_simpb_before);

     % compute tensions and curvatures and strain
    [E_sim_pb,~,~,A_after_pb,A_before_pb] = ...
        extract_sim_properties(tri_simpb_before,tri_simpb_after,subset_simpb_after,subset_simpb_before,tip,finalpb128_r,finalpb_r,spacing,magni);
    [~,~,sdir_input_after,cdir_input_after] = stretch_directions_unturgid(tri_simpb_after,subset_simpb_before,subset_simpb_after,[1;0;0],initpb_r);
    [mer_stretch_pb,circ_stretch_pb,sdir_input,cdir_input,~,~,sdir_after_test] = stretch_directions(tri_simpb_before,subset_simpb_before,subset_simpb_after,[1;0;0],finalpb_r,sdir_input_after,cdir_input_after);
        
    [curvs_tensions_pb,~,ref_km,ref_kc,ref_tm,ref_tc,ref_K,ref_MU,angles_sub_all,norm_comp_all,...
    xq_all,yq_all,zq_all] = curvature_tension_calc_sim(finalpb128_r,tri_simpb_before,subset_simpb_before,IC,spacing,magni,...
    sdir_input,cdir_input,circ_stretch_pb,mer_stretch_pb,Master_tri_i,structCurv_info);
    
    % location marker - not used
    [mm_min,mm_max] = max_min_t(tri_simpb_before,subset_simpb_before,1);
    % angle information
    [angles_T,~,mm_angle_sim] = get_thetas_sim_new(tri_simpb_after,initpb128_r,subset_simpb_after,tip,1,sym,long_axis);
    
   
    % compute K and mu
    [K_pb, mu_pb] = compute_K_mu_indiv(curvs_tensions_pb,[E_sim_pb(1,:);mer_stretch_pb;circ_stretch_pb]);
    [K_new_pb] = (curvs_tensions_pb(7,:)+curvs_tensions_pb(8,:))./(2*(A_before_pb./A_after_pb-1));

    % second triangulation
    [subset_simpb_before2,subset_simpb_after2,tri_simpb_before2,tri_simpb_after2] = triangulation_sim(initpb_r,finalpb_r,locpb2,Master_con1);
    tip = max(finalpb_r(1,:));
    IC2 = incenter(tri_simpb_before2);
    
    [E_sim_pb2,~,~,A_after_pb2,A_before_pb2] = ...
        extract_sim_properties(tri_simpb_before2,tri_simpb_after2,subset_simpb_after2,subset_simpb_before2,tip,finalpb128_r,finalpb_r,spacing,magni);
    [~,~,sdir_input_after2,cdir_input_after2] = stretch_directions_unturgid(tri_simpb_after2,subset_simpb_before2,subset_simpb_after2,[1;0;0],initpb_r);
    [mer_stretch_pb2,circ_stretch_pb2,sdir_input2,cdir_input2,~,~,sdir_after_test2] = stretch_directions(tri_simpb_before2,subset_simpb_before2,subset_simpb_after2,[1;0;0],finalpb_r,sdir_input_after2,cdir_input_after2);
        
    [curvs_tensions_pb2,~,ref_km2,ref_kc2,ref_tm2,ref_tc2,ref_K2,ref_MU2,angles_sub_all2,norm_comp_all2,...
    xq_all2,yq_all2,zq_all2] = curvature_tension_calc_sim(finalpb128_r,tri_simpb_before2,subset_simpb_before2,IC2,spacing,magni,...
    sdir_input2,cdir_input2,circ_stretch_pb2,mer_stretch_pb2,Master_tri_i2,structCurv_info);

    [mm_min1,mm_max2] = max_min_t(tri_simpb_before2,subset_simpb_before2,1);
    [angles_T2,~,mm_angle_sim2] = get_thetas_sim_new(tri_simpb_after2,initpb128_r,subset_simpb_after2,tip,1,sym,long_axis);
    
    [K_pb2, mu_pb2] = compute_K_mu_indiv(curvs_tensions_pb2,[E_sim_pb2(1,:);mer_stretch_pb2;circ_stretch_pb2]);
    [K_new_pb2] = (curvs_tensions_pb2(7,:)+curvs_tensions_pb2(8,:))./(2*(A_before_pb2./A_after_pb2-1));

    % subsetting negative K and triangles tagged by sensitivity study
        posK_mu = find(K_pb>0 & Master_i_tri);
        posK_mu2 = find(K_pb2>0 & Master_i_tri2);
        % posK_mu = find(Master_i_tri);
        % posK_mu2 = find(Master_i_tri2);
         
         mm_angle_sim_sub = mm_angle_sim(:,posK_mu);
         mm_angle_sim_sub2 = mm_angle_sim2(:,posK_mu2);

     [~,s_sim_sub] = sort([mean(mm_angle_sim_sub) mean(mm_angle_sim_sub2)]);

    mm_angle_sim_sub = mm_angle_sim(:,Master_tri_i);
    mm_angle_sim_sub2 = mm_angle_sim2(:,Master_tri_i2);

     [~,s_sim] = sort([mean(mm_angle_sim_sub) mean(mm_angle_sim_sub2)]);


    [~,~,sdir_input_after] = stretch_directions_unturgid(tri_simpb_after,subset_simpb_before,subset_simpb_after,[1;0;0],initpb_r);
    ang_b_a = acosd(dot(sdir_input_after,sdir_after_test));
    ang_b_a(ang_b_a>90) = 180-ang_b_a(ang_b_a>90);

    [~,~,sdir_input_after2] = stretch_directions_unturgid(tri_simpb_after2,subset_simpb_before2,subset_simpb_after2,[1;0;0],initpb_r);
    ang_b_a2 = acosd(dot(sdir_input_after2,sdir_after_test2));
    ang_b_a2(ang_b_a2>90) = 180-ang_b_a2(ang_b_a2>90);

% binning data
[mv_angles_circ_curv,bins_mid_5,angle_bins_5] = bin_thetas_sim([curvs_tensions_pb(1,posK_mu) curvs_tensions_pb2(1,posK_mu2)],[angles_T(posK_mu) angles_T2(posK_mu2)],binsA,1,sym,s_sim_sub,1,1,N_min,N_max);
[mv_angles_mer_curv,~] = bin_thetas_sim([curvs_tensions_pb(2,posK_mu) curvs_tensions_pb2(2,posK_mu2)],[angles_T(posK_mu) angles_T2(posK_mu2)],binsA,1,sym,s_sim_sub,1,1,N_min,N_max);
[mv_angles_circ_tens] = bin_thetas_sim([curvs_tensions_pb(8,posK_mu) curvs_tensions_pb2(8,posK_mu2)],[angles_T(posK_mu) angles_T2(posK_mu2)],binsA,1,sym,s_sim_sub,1,1,N_min,N_max);
[mv_angles_mer_tens,~] = bin_thetas_sim([curvs_tensions_pb(7,posK_mu) curvs_tensions_pb2(7,posK_mu2)],[angles_T(posK_mu) angles_T2(posK_mu2)],binsA,1,sym,s_sim_sub,1,1,N_min,N_max);

[mv_angles_circ_stretch,~,~,~,angle_circ_stretch,tri_circ_stretch] = bin_thetas_sim([circ_stretch_pb(posK_mu) circ_stretch_pb2(posK_mu2)],[angles_T(posK_mu) angles_T2(posK_mu2)],binsA,1,sym,s_sim_sub,1,1,N_min,N_max);
[mv_angles_mer_stretch,~] = bin_thetas_sim([mer_stretch_pb(posK_mu) mer_stretch_pb2(posK_mu2)],[angles_T(posK_mu) angles_T2(posK_mu2)],binsA,1,sym,s_sim_sub,1,1,N_min,N_max);
[mv_angles_K,~] = bin_thetas_sim([K_pb(posK_mu) K_pb2(posK_mu2)],[angles_T(posK_mu) angles_T2(posK_mu2)],binsA,1,sym,s_sim_sub,1,1,N_min,N_max);
[mv_angles_area_K,~,~,~] = bin_thetas_sim([K_new_pb(posK_mu) K_new_pb2(posK_mu2)],[angles_T(posK_mu) angles_T2(posK_mu2)],binsA,1,sym,s_sim_sub,1,1,N_min,N_max);
[mv_angles_ang_b_a,~] = bin_thetas_sim([ang_b_a(posK_mu) ang_b_a2(posK_mu2)],[angles_T(posK_mu) angles_T2(posK_mu2)],binsA,0,sym,s_sim_sub,1,1,N_min,N_max);


        posK_mu = find(mu_pb>0 & Master_i_tri);
        posK_mu2 = find(mu_pb2>0 & Master_i_tri2);
        % posK_mu = find(Master_i_tri);
        % posK_mu2 = find(Master_i_tri2);
        mm_angle_sim_sub = mm_angle_sim(:,posK_mu);
        mm_angle_sim_sub2 = mm_angle_sim2(:,posK_mu2);
        [~,s_sim_sub] = sort([mean(mm_angle_sim_sub) mean(mm_angle_sim_sub2)]);

[mv_angles_mu,~] = bin_thetas_sim([mu_pb(posK_mu) mu_pb2(posK_mu2)],[angles_T(posK_mu) angles_T2(posK_mu2)],binsA,1,sym,s_sim_sub,1,1,N_min,N_max);

% [mv_angles_stretch1_Kfilt,~] = bin_thetas_sim([E_sim_pb(3,posK_mu) E_sim_pb2(3,posK_mu2)],[angles_T(posK_mu) angles_T2(posK_mu2)],binsA,1,sym,s_sim_sub,1,1,N_min,N_max);
% [mv_angles_stretch2_Kfilt,~] = bin_thetas_sim([E_sim_pb(2,posK_mu) E_sim_pb2(2,posK_mu2)],[angles_T(posK_mu) angles_T2(posK_mu2)],binsA,1,sym,s_sim_sub,1,1,N_min,N_max);


[point_bin_circ_curv,mid_point,angle_point,~,~]=bin_raw_curv([ref_kc,ref_kc2],[angles_sub_all,angles_sub_all2],0,100,11);
[point_bin_mer_curv]=bin_raw_curv([ref_km,ref_km2],[angles_sub_all,angles_sub_all2],0,100,11);
[point_bin_circ_tens,~,~,~,~]=bin_raw_curv([ref_tc,ref_tc2],[angles_sub_all,angles_sub_all2],0,100,11);
[point_bin_mer_tens]=bin_raw_curv([ref_tm,ref_tm2],[angles_sub_all,angles_sub_all2],0,100,11);
point_posK_mu = ref_K >0 & ref_MU >0;
point_posK_mu2 = ref_K2 >0 & ref_MU2 >0;
[point_bin_K,~,~,~,~]=bin_raw_curv([ref_K(point_posK_mu),ref_K2(point_posK_mu2)],[angles_sub_all(point_posK_mu),angles_sub_all2(point_posK_mu2)],0,100,11);
[point_bin_MU]=bin_raw_curv([ref_MU(point_posK_mu),ref_MU2(point_posK_mu2)],[angles_sub_all(point_posK_mu),angles_sub_all2(point_posK_mu2)],0,100,11);
% [point_bin_K,~,~,~,~]=bin_raw_curv([ref_K,ref_K2],[angles_sub_all,angles_sub_all2],0,100,11);
% [point_bin_MU]=bin_raw_curv([ref_MU,ref_MU2],[angles_sub_all,angles_sub_all2],0,100,11);




        posK_mu = find(K_pb>0 & Master_i_tri & ang_b_a<7);
        posK_mu2 = find(K_pb2>0 & Master_i_tri2 & ang_b_a2<7);
         
         mm_angle_sim_sub = mm_angle_sim(:,posK_mu);
         mm_angle_sim_sub2 = mm_angle_sim2(:,posK_mu2);

[~,s_sim_sub] = sort([mean(mm_angle_sim_sub) mean(mm_angle_sim_sub2)]);
[mv_angles_K_angsub,~] = bin_thetas_sim([K_pb(posK_mu) K_pb2(posK_mu2)],[angles_T(posK_mu) angles_T2(posK_mu2)],binsA,1,sym,s_sim_sub,1,1,N_min,N_max);
[mv_angles_ang_dev_sub,~] = bin_thetas_sim([ang_b_a(posK_mu) ang_b_a2(posK_mu2)],[angles_T(posK_mu) angles_T2(posK_mu2)],binsA,0,sym,s_sim_sub,1,1,N_min,N_max);
[mv_angles_circ_stretch_angsub,~] = bin_thetas_sim([circ_stretch_pb(posK_mu) circ_stretch_pb2(posK_mu2)],[angles_T(posK_mu) angles_T2(posK_mu2)],binsA,1,sym,s_sim_sub,1,1,N_min,N_max);
[mv_angles_mer_stretch_angsub,~] = bin_thetas_sim([mer_stretch_pb(posK_mu) mer_stretch_pb2(posK_mu2)],[angles_T(posK_mu) angles_T2(posK_mu2)],binsA,1,sym,s_sim_sub,1,1,N_min,N_max);

[point_bin_circ_curv,mid_point,angle_point,~,~]=bin_raw_curv([ref_kc,ref_kc2],[angles_sub_all,angles_sub_all2],0,100,11);
[point_bin_mer_curv]=bin_raw_curv([ref_km,ref_km2],[angles_sub_all,angles_sub_all2],0,100,11);
[point_bin_circ_tens,~,~,~,~]=bin_raw_curv([ref_tc,ref_tc2],[angles_sub_all,angles_sub_all2],0,100,11);
[point_bin_mer_tens]=bin_raw_curv([ref_tm,ref_tm2],[angles_sub_all,angles_sub_all2],0,100,11);
point_posK_mu = ref_K >0 & ref_MU >0;
point_posK_mu2 = ref_K2 >0 & ref_MU2 >0;
[point_bin_K,~,~,~,~]=bin_raw_curv([ref_K(point_posK_mu),ref_K2(point_posK_mu2)],[angles_sub_all(point_posK_mu),angles_sub_all2(point_posK_mu2)],0,100,11);
[point_bin_MU]=bin_raw_curv([ref_MU(point_posK_mu),ref_MU2(point_posK_mu2)],[angles_sub_all(point_posK_mu),angles_sub_all2(point_posK_mu2)],0,100,11);


    %record data
    K128_trials(trialCount,:) = K128';
    MU128_trials(trialCount,:) =  MU128';
    synCurvM128_trials(trialCount,:) = synCurvM128';
    synCurvC128_trials(trialCount,:) = synCurvC128';
    synTensM128_trials(trialCount,:) = synTensM128';
    synTensC128_trials(trialCount,:) = synTensC128';
    synStrainM128_trials(trialCount,:) = synStrainM128';
    synStrainC128_trials(trialCount,:) = synStrainC128';
    rv_final128_1_trials(trialCount,:) = rv_final128(1,:);
    rv_final128_2_trials(trialCount,:) = rv_final128(2,:);
    scale_all(trialCount,:) = scale;

    infTensMs{trialCount}=[infTensMs{trialCount} abs(curvs_tensions_pb(7,Master_tri_i)) abs(curvs_tensions_pb2(7,Master_tri_i2))];
    infTensCs{trialCount}=[infTensCs{trialCount} abs(curvs_tensions_pb(8,Master_tri_i)) abs(curvs_tensions_pb2(8,Master_tri_i2))];
    infCurvMs{trialCount} = [infCurvMs{trialCount} abs(curvs_tensions_pb(4,Master_tri_i)) abs(curvs_tensions_pb2(4,Master_tri_i2))];
    infCurvCs{trialCount} = [infCurvCs{trialCount} abs(curvs_tensions_pb(3,Master_tri_i)) abs(curvs_tensions_pb2(3,Master_tri_i2))];

    infStrainMs{trialCount}=[infStrainMs{trialCount} mer_stretch_pb(Master_tri_i) mer_stretch_pb2(Master_tri_i2)];
    infStrainCs{trialCount}=[infStrainCs{trialCount} circ_stretch_pb(Master_tri_i) circ_stretch_pb2(Master_tri_i2)];
    infStrainArea{trialCount}= E_sim_pb(2,Master_tri_i).*E_sim_pb(3,Master_tri_i);

    K_pb_final = [K_pb(Master_tri_i) K_pb2(Master_tri_i2)];
    K_pb_final(~posK_mu) = NaN;

    mu_pb_final = [mu_pb(Master_tri_i) mu_pb2(Master_tri_i2)];
    mu_pb_final(~posK_mu) = NaN;

    infKs_124{trialCount}=[infKs_124{trialCount} pre_process(K_pb_final,1,size(K_pb_final,2),s_sim,1)];
    infMUs_124{trialCount}=[infMUs_124{trialCount} pre_process(mu_pb_final,1,size(mu_pb_final,2),s_sim,1)];

    final_angles_T=[angles_T(Master_tri_i) angles_T2(Master_tri_i2)];
    angles_T_all{trialCount}=[angles_T_all{trialCount} final_angles_T(:,s_sim)];
    mm_angles_all{trialCount}=[mm_angles_all{trialCount} mean(mm_angle_sim_sub) mean(mm_angle_sim_sub2)];
%}




%
point_bin_circ_curv_nboth_t100(trialCount,:) = point_bin_circ_curv(1,:)*max(rv_final128(2,:));
point_bin_mer_curv_nboth_t100(trialCount,:) = point_bin_mer_curv(1,:)*max(rv_final128(2,:));
point_bin3_circ_curv_nboth_t100(trialCount,:) = point_bin_circ_curv(3,:);
point_bin3_mer_curv_nboth_t100(trialCount,:) = point_bin_mer_curv(3,:);

point_bin_circ_tens_nboth_t100(trialCount,:) = point_bin_circ_tens(1,:)/max(rv_final128(2,:));
point_bin_mer_tens_nboth_t100(trialCount,:) = point_bin_mer_tens(1,:)/max(rv_final128(2,:));
point_bin3_circ_tens_nboth_t100(trialCount,:) = point_bin_circ_tens(3,:);
point_bin3_mer_tens_nboth_t100(trialCount,:) = point_bin_mer_tens(3,:);

point_bin_K_nboth_t100(trialCount,:) = point_bin_K(1,:)/max(rv_init128(2,:));
point_bin_mu_nboth_t100(trialCount,:) = point_bin_MU(1,:)/max(rv_init128(2,:));
point_bin3_K_nboth_t100(trialCount,:) = point_bin_K(3,:);
point_bin3_mu_nboth_t100(trialCount,:) = point_bin_MU(3,:);

        mv_angles_ang_b_a_nboth_t100(trialCount,:) = mv_angles_ang_b_a(1,:);
        mv_angles_K_nboth_t100(trialCount,:) = mv_angles_K(1,:)/max(rv_init128(2,:)); % normalize
        mv_angles_area_K_nboth_t100(trialCount,:) = mv_angles_area_K(1,:)/max(rv_init128(2,:)); % normalize
        mv_angles2_K_nboth_t100(trialCount,:) = mv_angles_K(2,:);
        mv_angles3_K_nboth_t100(trialCount,:) = mv_angles_K(3,:);
        mv_angles4_K_nboth_t100(trialCount,:) = mv_angles_K(4,:);
        mv_angles_mu_nboth_t100(trialCount,:) = mv_angles_mu(1,:);
        mv_angles2_mu_nboth_t100(trialCount,:) = mv_angles_mu(2,:);
        mv_angles3_mu_nboth_t100(trialCount,:) = mv_angles_mu(3,:);
    % [mv_angles_circ_stretch_nboth,~] = bin_thetas_sim(circ_stretch_pb(Master_tri_i),angles_T(Master_tri_i),binsA,1,sym,s_sim,1,1,N_min,N_max);
    % [mv_angles_mer_stretch_nboth,~] = bin_thetas_sim(mer_stretch_pb(Master_tri_i),angles_T(Master_tri_i),binsA,1,sym,s_sim,1,1,N_min,N_max);
    % [mv_angles_area_nboth,~] = bin_thetas_sim(E_sim_pb(2,Master_tri_i).*E_sim_pb(3,Master_tri_i),angles_T(Master_tri_i),binsA,1,sym,s_sim,1,1,N_min,N_max);
        mv_angles_circ_stretch_nboth_t100(trialCount,:) = mv_angles_circ_stretch(1,:);
        mv_angles3_circ_stretch_nboth_t100(trialCount,:) = mv_angles_circ_stretch(3,:);
        mv_angles_mer_stretch_nboth_t100(trialCount,:) = mv_angles_mer_stretch(1,:);
        mv_angles3_mer_stretch_nboth_t100(trialCount,:) = mv_angles_mer_stretch(3,:);
        % mv_angles_area_nboth_t100(trialCount,:) = mv_angles_area_nboth(1,:);
    % [mv_angles_circ_tens_nboth,~] = bin_thetas_sim(curvs_tensions_pb(8,Master_tri_i),angles_T(Master_tri_i),binsA,1,sym,s_sim,1,1,N_min,N_max);
    % [mv_angles_mer_tens_nboth,~] = bin_thetas_sim(curvs_tensions_pb(7,Master_tri_i),angles_T(Master_tri_i),binsA,1,sym,s_sim,1,1,N_min,N_max);
        mv_angles_circ_tens_nboth_t100(trialCount,:) = mv_angles_circ_tens(1,:)/max(rv_final128(2,:));
        mv_angles3_circ_tens_nboth_t100(trialCount,:) = mv_angles_circ_tens(3,:);
        mv_angles_mer_tens_nboth_t100(trialCount,:) = mv_angles_mer_tens(1,:)/max(rv_final128(2,:));
        mv_angles3_mer_tens_nboth_t100(trialCount,:) = mv_angles_mer_tens(3,:);
    % [mv_angles_circ_curv_nboth,~] = bin_thetas_sim(abs(curvs_tensions_pb(3,Master_tri_i)),angles_T(Master_tri_i),binsA,1,sym,s_sim,1,1,N_min,N_max);
    % [mv_angles_mer_curv_nboth,~] = bin_thetas_sim(abs(curvs_tensions_pb(4,Master_tri_i)),angles_T(Master_tri_i),binsA,1,sym,s_sim,1,1,N_min,N_max);
        mv_angles_circ_curv_nboth_t100(trialCount,:) = mv_angles_circ_curv(1,:)*max(rv_final128(2,:));
        mv_angles3_circ_curv_nboth_t100(trialCount,:) = mv_angles_circ_curv(3,:);
        mv_angles_mer_curv_nboth_t100(trialCount,:) = mv_angles_mer_curv(1,:)*max(rv_final128(2,:));
        mv_angles3_mer_curv_nboth_t100(trialCount,:) = mv_angles_mer_curv(3,:);
    % opti_K_nboth = (mv_angles_circ_tens(1,:)+mv_angles_mer_tens(1,:))./(2.*((mv_angles_circ_stretch(1,:).*mv_angles_mer_stretch(1,:))-1));
    % opti_mu_nboth = (mv_angles_mer_tens(1,:)-mv_angles_circ_tens(1,:))./(1./(mv_angles_circ_stretch(1,:).^2)-(1./(mv_angles_mer_stretch(1,:).^2)));
  

         nofilt= find(Master_i_tri);
         nofilt2 = find(Master_i_tri2);
         
         mm_angle_sim_sub = mm_angle_sim(:,nofilt);
         mm_angle_sim_sub2 = mm_angle_sim2(:,nofilt2);

     [~,s_sim_sub] = sort([mean(mm_angle_sim_sub) mean(mm_angle_sim_sub2)]);

    mm_angle_sim_sub = mm_angle_sim(:,Master_tri_i);
    mm_angle_sim_sub2 = mm_angle_sim2(:,Master_tri_i2);

     [~,s_sim] = sort([mean(mm_angle_sim_sub) mean(mm_angle_sim_sub2)]);
        % mm_angle_sim_sub = mm_angle_sim(:,posK_mu);
        % mm_angle_sim_sub2 = mm_angle_sim2(:,posK_mu2);
        % [~,s_sim] = sort([mean(mm_angle_sim_sub) mean(mm_angle_sim_sub2)]);
    % end
mm_angle_all_nboth = [mean(mm_angle_sim_sub) mean(mm_angle_sim_sub2)];
    
[mvNF_angles_circ_curv,bins_angles_mid,angle_bins] = bin_thetas_sim([curvs_tensions_pb(1,nofilt) curvs_tensions_pb2(1,nofilt2)],[angles_T(nofilt) angles_T2(nofilt2)],binsA,1,sym,s_sim_sub,1,1,N_min,N_max);
[mvNF_angles_mer_curv,~] = bin_thetas_sim([curvs_tensions_pb(2,nofilt) curvs_tensions_pb2(2,nofilt2)],[angles_T(nofilt) angles_T2(nofilt2)],binsA,1,sym,s_sim_sub,1,1,N_min,N_max);
[mvNF_angles_circ_tens] = bin_thetas_sim([curvs_tensions_pb(8,nofilt) curvs_tensions_pb2(8,nofilt2)],[angles_T(nofilt) angles_T2(nofilt2)],binsA,1,sym,s_sim_sub,1,1,N_min,N_max);
[mvNF_angles_mer_tens,~] = bin_thetas_sim([curvs_tensions_pb(7,nofilt) curvs_tensions_pb2(7,nofilt2)],[angles_T(nofilt) angles_T2(nofilt2)],binsA,1,sym,s_sim_sub,1,1,N_min,N_max);

[mvNF_angles_circ_stretch,~] = bin_thetas_sim([circ_stretch_pb(nofilt) circ_stretch_pb2(nofilt2)],[angles_T(nofilt) angles_T2(nofilt2)],binsA,1,sym,s_sim_sub,1,1,N_min,N_max);
[mvNF_angles_mer_stretch,~] = bin_thetas_sim([mer_stretch_pb(nofilt) mer_stretch_pb2(nofilt2)],[angles_T(nofilt) angles_T2(nofilt2)],binsA,1,sym,s_sim_sub,1,1,N_min,N_max);
[mvNF_angles_K,~] = bin_thetas_sim([K_pb(nofilt) K_pb2(nofilt2)],[angles_T(nofilt) angles_T2(nofilt2)],binsA,1,sym,s_sim_sub,1,1,N_min,N_max);
[mvNF_angles_mu,~] = bin_thetas_sim([mu_pb(nofilt) mu_pb2(nofilt2)],[angles_T(nofilt) angles_T2(nofilt2)],binsA,1,sym,s_sim_sub,1,1,N_min,N_max);
[mvNF_angles_area_K,~,~,~] = bin_thetas_sim([K_new_pb(nofilt) K_new_pb2(nofilt2)],[angles_T(nofilt) angles_T2(nofilt2)],binsA,1,sym,s_sim_sub,1,1,N_min,N_max);

[mvNF_angles_ang_b_a,~] = bin_thetas_sim([ang_b_a(nofilt) ang_b_a2(nofilt2)],[angles_T(nofilt) angles_T2(nofilt2)],binsA,0,sym,s_sim_sub,1,1,N_min,N_max);

% [mv_angles_stretch1_nboth,~] = bin_thetas_sim([E_sim_pb(3,posK_mu) E_sim_pb2(3,posK_mu2)],[angles_T(posK_mu) angles_T2(posK_mu2)],binsA,1,sym,s_sim_sub,1,1,N_min,N_max);
% [mv_angles_stretch2_nboth,~] = bin_thetas_sim([E_sim_pb(2,posK_mu) E_sim_pb2(2,posK_mu2)],[angles_T(posK_mu) angles_T2(posK_mu2)],binsA,1,sym,s_sim_sub,1,1,N_min,N_max);


% point_posK_mu = ref_K >0 & ref_MU >0;
% point_posK_mu2 = ref_K2 >0 & ref_MU2 >0;
[pointNF_bin_K,~,~,~,~]=bin_raw_curv([ref_K,ref_K2],[angles_sub_all,angles_sub_all2],0,100,11);
[pointNF_bin_MU]=bin_raw_curv([ref_MU,ref_MU2],[angles_sub_all,angles_sub_all2],0,100,11);

        nofilt = find(Master_i_tri & ang_b_a<7);
        nofilt2 = find(Master_i_tri2 & ang_b_a2<7);
         
         mm_angle_sim_sub = mm_angle_sim(:,nofilt);
         mm_angle_sim_sub2 = mm_angle_sim2(:,nofilt2);

[~,s_sim_sub] = sort([mean(mm_angle_sim_sub) mean(mm_angle_sim_sub2)]);
[mvNF_angles_K_angsub,~] = bin_thetas_sim([K_pb(nofilt) K_pb2(nofilt2)],[angles_T(nofilt) angles_T2(nofilt2)],binsA,1,sym,s_sim_sub,1,1,N_min,N_max);
[mvNF_angles_ang_dev_sub,~] = bin_thetas_sim([ang_b_a(nofilt) ang_b_a2(nofilt2)],[angles_T(nofilt) angles_T2(nofilt2)],binsA,0,sym,s_sim_sub,1,1,N_min,N_max);
[mvNF_angles_circ_stretch_angsub,~] = bin_thetas_sim([circ_stretch_pb(nofilt) circ_stretch_pb2(nofilt2)],[angles_T(nofilt) angles_T2(nofilt2)],binsA,1,sym,s_sim_sub,1,1,N_min,N_max);
[mvNF_angles_mer_stretch_angsub,~] = bin_thetas_sim([mer_stretch_pb(nofilt) mer_stretch_pb2(nofilt2)],[angles_T(nofilt) angles_T2(nofilt2)],binsA,1,sym,s_sim_sub,1,1,N_min,N_max);



    K_pb_final = [K_pb(Master_tri_i) K_pb2(Master_tri_i2)];
    K_pb_final(~nofilt) = NaN;

    mu_pb_final = [mu_pb(Master_tri_i) mu_pb2(Master_tri_i2)];
    mu_pb_final(~nofilt) = NaN;

    infKs_124{trialCount}=[infKs_124{trialCount} pre_process(K_pb_final,1,size(K_pb_final,2),s_sim,1)];
    infMUs_124{trialCount}=[infMUs_124{trialCount} pre_process(mu_pb_final,1,size(mu_pb_final,2),s_sim,1)];

    final_angles_T=[angles_T(Master_tri_i) angles_T2(Master_tri_i2)];
    angles_T_all{trialCount}=[angles_T_all{trialCount} final_angles_T(:,s_sim)];
    mm_angles_all{trialCount}=[mm_angles_all{trialCount} mean(mm_angle_sim_sub) mean(mm_angle_sim_sub2)];
%

pointNF_bin_K_nboth_t100(trialCount,:) = pointNF_bin_K(1,:)/max(rv_init128(2,:));
pointNF_bin_mu_nboth_t100(trialCount,:) = pointNF_bin_MU(1,:)/max(rv_init128(2,:));
pointNF_bin3_K_nboth_t100(trialCount,:) = pointNF_bin_K(3,:);
pointNF_bin3_mu_nboth_t100(trialCount,:) = pointNF_bin_MU(3,:);

        mvNF_angles_ang_b_a_nboth_t100(trialCount,:) = mvNF_angles_ang_b_a(1,:);
        mvNF_angles_K_nboth_t100(trialCount,:) = mvNF_angles_K(1,:)/max(rv_init128(2,:)); % normalize
        mvNF_angles_area_K_nboth_t100(trialCount,:) = mvNF_angles_area_K(1,:)/max(rv_init128(2,:)); % normalize
        mvNF_angles2_K_nboth_t100(trialCount,:) = mvNF_angles_K(2,:);
        mvNF_angles3_K_nboth_t100(trialCount,:) = mvNF_angles_K(3,:);
        mvNF_angles_mu_nboth_t100(trialCount,:) = mvNF_angles_mu(1,:);
        mvNF_angles2_mu_nboth_t100(trialCount,:) = mvNF_angles_mu(2,:);
        mvNF_angles3_mu_nboth_t100(trialCount,:) = mvNF_angles_mu(3,:);
        mvNF_angles_circ_stretch_nboth_t100(trialCount,:) = mvNF_angles_circ_stretch(1,:);
        mvNF_angles3_circ_stretch_nboth_t100(trialCount,:) = mvNF_angles_circ_stretch(3,:);
        mvNF_angles_mer_stretch_nboth_t100(trialCount,:) = mvNF_angles_mer_stretch(1,:);
        mvNF_angles3_mer_stretch_nboth_t100(trialCount,:) = mvNF_angles_mer_stretch(3,:);
        mvNF_angles_circ_tens_nboth_t100(trialCount,:) = mvNF_angles_circ_tens(1,:)/max(rv_final128(2,:));
        mvNF_angles3_circ_tens_nboth_t100(trialCount,:) = mvNF_angles_circ_tens(3,:);
        mvNF_angles_mer_tens_nboth_t100(trialCount,:) = mvNF_angles_mer_tens(1,:)/max(rv_final128(2,:));
        mvNF_angles3_mer_tens_nboth_t100(trialCount,:) = mvNF_angles_mer_tens(3,:);
        mvNF_angles_circ_curv_nboth_t100(trialCount,:) = mvNF_angles_circ_curv(1,:)*max(rv_final128(2,:));
        mvNF_angles3_circ_curv_nboth_t100(trialCount,:) = mvNF_angles_circ_curv(3,:);
        mvNF_angles_mer_curv_nboth_t100(trialCount,:) = mvNF_angles_mer_curv(1,:)*max(rv_final128(2,:));
        mvNF_angles3_mer_curv_nboth_t100(trialCount,:) = mvNF_angles_mer_curv(3,:);

        mvNF_angles4_K_nboth_t100(trialCount,:) = mvNF_angles_K(4,:);

end


cell_Kfilt = {bins_angles_mid,angle_bins,mv_angles_K_nboth_t100,mv_angles2_K_nboth_t100,mv_angles_mu_nboth_t100,...
    mv_angles2_mu_nboth_t100,mv_angles_circ_curv_nboth_t100,mv_angles_circ_stretch_nboth_t100,...
    mv_angles3_circ_stretch_nboth_t100,mv_angles_mer_stretch_nboth_t100,mv_angles3_mer_stretch_nboth_t100,...
    mv_angles3_circ_curv_nboth_t100,mv_angles_mer_curv_nboth_t100,mv_angles3_mer_curv_nboth_t100,...
    mv_angles_K,mm_angle_all_nboth,mv_angles_mu,mv_angles_circ_stretch,mv_angles_mer_stretch,...
    mv_angles_mer_tens,mv_angles_circ_tens,mv_angles_circ_curv,mv_angles_mer_curv,mm_angles_all,mv_angles3_K_nboth_t100,...
    mv_angles3_mu_nboth_t100,mv_angles_circ_tens_nboth_t100,mv_angles3_circ_tens_nboth_t100,...
    mv_angles_mer_tens_nboth_t100,mv_angles3_mer_tens_nboth_t100,...
    point_bin_circ_curv,point_bin_mer_curv,point_bin_circ_tens,point_bin_mer_tens,point_bin_K,point_bin_MU,...
    point_bin_K_nboth_t100,mid_point,angle_point,scale_all,magni,magni_init,magni_z,...
    mv_angles4_K_nboth_t100,mv_angles_area_K_nboth_t100,mv_angles_ang_b_a_nboth_t100};
cell_nofilt = {bins_angles_mid,angle_bins,mvNF_angles_K_nboth_t100,mvNF_angles2_K_nboth_t100,mvNF_angles_mu_nboth_t100,...
    mvNF_angles2_mu_nboth_t100,mvNF_angles_circ_curv_nboth_t100,mvNF_angles_circ_stretch_nboth_t100,...
    mvNF_angles3_circ_stretch_nboth_t100,mvNF_angles_mer_stretch_nboth_t100,mvNF_angles3_mer_stretch_nboth_t100,...
    mvNF_angles3_circ_curv_nboth_t100,mvNF_angles_mer_curv_nboth_t100,mvNF_angles3_mer_curv_nboth_t100,...
    mvNF_angles_K,mm_angle_all_nboth,mvNF_angles_mu,mvNF_angles_circ_stretch,mvNF_angles_mer_stretch,...
    mvNF_angles_mer_tens,mvNF_angles_circ_tens,mvNF_angles_circ_curv,mvNF_angles_mer_curv,mm_angles_all,mvNF_angles3_K_nboth_t100,...
    mvNF_angles3_mu_nboth_t100,mvNF_angles_circ_tens_nboth_t100,mvNF_angles3_circ_tens_nboth_t100,...
    mvNF_angles_mer_tens_nboth_t100,mvNF_angles3_mer_tens_nboth_t100,...
    point_bin_circ_curv,point_bin_mer_curv,point_bin_circ_tens,point_bin_mer_tens,pointNF_bin_K,pointNF_bin_MU,...
    pointNF_bin_K_nboth_t100,mid_point,angle_point,scale_all,magni,magni_init,magni_z,...
    mvNF_angles4_K_nboth_t100,mvNF_angles_area_K_nboth_t100,mvNF_angles_ang_b_a_nboth_t100};
fields = ["bins_angles_mid","angle_bins","mv_angles_K_nboth_t100","mv_angles2_K_nboth_t100","mv_angles_mu_nboth_t100",...
    "mv_angles2_mu_nboth_t100","mv_angles_circ_curv_nboth_t100","mv_angles_circ_stretch_nboth_t100",...
    "mv_angles3_circ_stretch_nboth_t100","mv_angles_mer_stretch_nboth_t100","mv_angles3_mer_stretch_nboth_t100",...
    "mv_angles3_circ_curv_nboth_t100","mv_angles_mer_curv_nboth_t100","mv_angles3_mer_curv_nboth_t100",...
    "mv_angles_K","mm_angle_all_nboth","mv_angles_mu","mv_angles_circ_stretch","mv_angles_mer_stretch",...
    "mv_angles_mer_tens","mv_angles_circ_tens","mv_angles_circ_curv","mv_angles_mer_curv","mm_angles_all","mv_angles3_K_nboth_t100",...
    "mv_angles3_mu_nboth_t100","mv_angles_circ_tens_nboth_t100","mv_angles3_circ_tens_nboth_t100",...
    "mv_angles_mer_tens_nboth_t100","mv_angles3_mer_tens_nboth_t100",...
    "point_bin_circ_curv","point_bin_mer_curv","point_bin_circ_tens","point_bin_mer_tens","point_bin_K","point_bin_MU",...
    "point_bin_K_nboth_t100","mid_point","angle_point",...
    "scale_all","magni","magni_init","magni_z","mv_angles4_K_nboth_t100","mv_angles_area_K_nboth_t100","mv_angles_ang_b_a_nboth_t100"];
structAll_Kfilt = cell2struct(cell_Kfilt',fields);
structAll_nofilt = cell2struct(cell_nofilt',fields);
end
% end
