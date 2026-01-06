function [curvs_tensions,princp_directions,ref_km,ref_kc,ref_tm,ref_tc,ref_K,ref_MU,angles_sub_all,norm_comp_all,xq_all,yq_all,zq_all] = curvature_tension_calc_exp(before_xy_ridge_all,tri_before,subset_before,IC,...
    tip_point,rad,spacing,long_axis,cdir_input,sdir_input,circ_stretch,mer_stretch,tri_i)

curvs_tensions = zeros(8,size(tri_before.ConnectivityList,1));
princp_directions = zeros(24,size(tri_before.ConnectivityList,1));

[~,i] = max(before_xy_ridge_all(1,:));
tip3 = before_xy_ridge_all(:,i);
% km = zeros(1,size(tri_before.ConnectivityList,1));
% kc = zeros(1,size(tri_before.ConnectivityList,1));
P = 2;
% count = 1;

[~,i] = max(before_xy_ridge_all(1,:));
mid_y = before_xy_ridge_all(2,i);
   
ref_km = [];
ref_kc = [];
ref_tc = [];
ref_tm = [];
ref_K = [];
ref_MU = [];
xq_all = [];
yq_all = [];
zq_all = [];
angles_sub_all = [];
norm_comp_all = cell(1,size(tri_before.ConnectivityList,1));




[angles1_all,angles2_all,P1_side1_all,P1_side2_all,anglestip_all,P1_tip_all,P2_tip_all,P2_side1_all,P2_side2_all,km_1_all,...
kc_1_all,km_2_all,kc_2_all,km_t_all,kc_t_all,xq1_all,yq1_all,zq1_all,xq2_all,yq2_all,zq2_all,xqt_all,yqt_all,zqt_all,...
VP1_side1_all,VP2_side1_all,VP1_side2_all,VP2_side2_all,VP1_tip_all,VP2_tip_all,mdirs_1_all,cdirs_1_all,mdirs_2_all,cdirs_2_all,mdirs_t_all,cdirs_t_all,...
normals1_all,normals2_all,normalst_all] = get_curv_raw_exp(before_xy_ridge_all,spacing,rad,subset_before,long_axis);




    % quiver3(xqt_all,yqt_all,zqt_all,mdirs_t_all(:,1),mdirs_t_all(:,2),mdirs_t_all(:,3),3,'r');axis equal
    % hold on; quiver3(xqt_all,yqt_all,zqt_all,cdirs_t_all(:,1),cdirs_t_all(:,2),cdirs_t_all(:,3),3,'b');axis equal
    % scatter3(xqt_all,yqt_all,zqt_all,10,km_t_all,'filled')
     % quiver3(xq1_all,yq1_all,zq1_all,mdirs_1_all(:,1),mdirs_1_all(:,2),mdirs_1_all(:,3),3);axis equal
% scatter(angles2_all,km_1_all,'r');hold on; scatter(angles2_all,km_2_all,'r');hold on; scatter(anglestip_all,km_t_all,'r')
% for i = 1:size(normals1_all,1)
        % t1_side = P ./ (2 * curv_1);
        % t2_side = P ./ (2 * curv_1) .* (2 - curv_2 ./ curv_1);
        tm_1_all = P ./ (2 * kc_1_all);
        tc_1_all = P ./ (2 * kc_1_all) .* (2 - km_1_all ./ kc_1_all);
        tm_2_all = P ./ (2 * kc_2_all);
        tc_2_all = P ./ (2 * kc_2_all) .* (2 - km_2_all ./ kc_2_all);
        tm_t_all = P ./ (2 * kc_t_all);
        tc_t_all = P ./ (2 * kc_t_all) .* (2 - km_t_all ./ kc_t_all);

        tc_1_all = pre_process(tc_1_all',1,0,0,0)';
        tm_1_all = pre_process(tm_1_all',1,0,0,0)';
        tc_2_all = pre_process(tc_2_all',1,0,0,0)';
        tm_2_all = pre_process(tm_2_all',1,0,0,0)';
        tc_t_all = pre_process(tc_t_all',1,0,0,0)';
        tm_t_all = pre_process(tm_t_all',1,0,0,0)';

        kc_1_all(isnan(tc_1_all)) = NaN;
        kc_2_all(isnan(tc_2_all)) = NaN;
        kc_t_all(isnan(tc_t_all)) = NaN;
        km_1_all(isnan(tm_1_all)) = NaN;
        km_2_all(isnan(tm_2_all)) = NaN;
        km_t_all(isnan(tm_t_all)) = NaN;
        % save('exp_test.mat')
for i = 1:size(VP1_side1_all,1)
    VP1_side1_all(i,:) = VP1_side1_all(i,:)/norm(VP1_side1_all(i,:));
    VP2_side1_all(i,:) = VP2_side1_all(i,:)/norm(VP2_side1_all(i,:));
    if VP1_side1_all(i,3) > 0
        VP1_side1_all(i,:) = -VP1_side1_all(i,:);
    else
    end
    % if VP2_side1_all(i,2) > 0 && VP2_side1_all(i,3) > 0
    %     VP2_side1_all(i,:) = -VP2_side1_all(i,:);
    % else
    % end
end
for i = 1:size(VP1_side2_all,1)
    VP1_side2_all(i,:) = VP1_side2_all(i,:)/norm(VP1_side2_all(i,:));
    VP2_side2_all(i,:) = VP2_side2_all(i,:)/norm(VP2_side2_all(i,:));
    if VP1_side2_all(i,3) > 0 
        VP1_side2_all(i,:) = -VP1_side2_all(i,:);
    else
    end
    % if VP2_side2_all(i,2) > 0 && VP2_side2_all(i,3) > 0
    %     VP2_side2_all(i,:) = -VP2_side2_all(i,:);
    % else
    % end
end
for i = 1:size(VP1_tip_all,1)
    VP1_tip_all(i,:) = VP1_tip_all(i,:)/norm(VP1_tip_all(i,:));
    VP2_tip_all(i,:) = VP2_tip_all(i,:)/norm(VP2_tip_all(i,:));
    if VP1_tip_all(i,3) > 0 
        VP1_tip_all(i,:) = -VP1_tip_all(i,:);
    else
    end
    % if VP1_tip_all(i,2) <0 && abs(VP1_tip_all(i,2)) >0.7
    %     VP1_tip_all(i,:) = -VP1_tip_all(i,:);
    % else
    % end

end
%%
norms = tri_before.faceNormal;
        for n = 1:size(tri_before.ConnectivityList,1)
    con_n = tri_before.ConnectivityList(n,:);
    points = tri_before.Points(con_n,:);

    e3 = long_axis;
    N_before = norms(n,:);
    v_b = cross(N_before,e3); c = dot(N_before,e3);
    v_x_b = [0 -v_b(3) v_b(2); v_b(3) 0 -v_b(1); -v_b(2) v_b(1) 0];
    Q_b = eye(3) + v_x_b + ((v_x_b^2)*(1/(1+c)));
    test_ridge = Q_b*before_xy_ridge_all;
    rotated_before = Q_b*subset_before;
    
    side1 = (yq1_all>=min(points(:,2))-5) & (yq1_all<=max(points(:,2))+5) & (xq1_all>=min(points(:,1))-5) & (xq1_all<=max(points(:,1))+5)...
        & (zq1_all>=min(points(:,3))-5) & (zq1_all<=max(points(:,3))+5) ;

    side2 = (yq2_all>=min(points(:,2))-5) & (yq2_all<=max(points(:,2))+5) & (xq2_all>=min(points(:,1))-5) & (xq2_all<=max(points(:,1))+5)...
        & (zq2_all>=min(points(:,3))-5) & (zq2_all<=max(points(:,3))+5) ;

    sidet = (yqt_all>=min(points(:,2))-5) & (yqt_all<=max(points(:,2))+5) & (xqt_all>=min(points(:,1))-5) & (xqt_all<=max(points(:,1))+5)...
        & (zqt_all>=min(points(:,3))-5) & (zqt_all<=max(points(:,3))+5) ;

    xq1_sub = xq1_all(side1); yq1_sub = yq1_all(side1);zq1_sub = zq1_all(side1);

    xqt_sub = xqt_all(sidet); yqt_sub = yqt_all(sidet);zqt_sub = zqt_all(sidet);

    xq2_sub = xq2_all(side2); yq2_sub = yq2_all(side2);zq2_sub = zq2_all(side2);

    side1_rot = Q_b*[xq1_sub';yq1_sub';zq1_sub'];
    side2_rot = Q_b*[xq2_sub';yq2_sub';zq2_sub'];
    sidet_rot = Q_b*[xqt_sub';yqt_sub';zqt_sub'];

    cdirs_1_sub = cdirs_1_all(side1,:);cdirs_2_sub = cdirs_2_all(side2,:);cdirs_t_sub = cdirs_t_all(sidet,:);

    mdirs_1_sub = mdirs_1_all(side1,:);mdirs_2_sub = mdirs_2_all(side2,:);mdirs_t_sub = mdirs_t_all(sidet,:);

    P1_side1_sub = P1_side1_all(side1);P1_side2_sub = P1_side2_all(side2);P1_tip_sub = P1_tip_all(sidet);

    P2_side1_sub = P2_side1_all(side1);P2_side2_sub = P2_side2_all(side2);P2_tip_sub = P2_tip_all(sidet);

    VP1_side1_sub = VP1_side1_all(side1,:);VP1_side2_sub = VP1_side2_all(side2,:);VP1_tip_sub = VP1_tip_all(sidet,:);

    VP2_side1_sub = VP2_side1_all(side1,:);VP2_side2_sub = VP2_side2_all(side2,:);VP2_tip_sub = VP2_tip_all(sidet,:);

    % VP1_side1_sub = VP1_side1_all(side1,:);VP1_side2_sub = VP1_side2_all(side2,:);VP1_tip_sub = VP1_tip_all(sidet,:);

    tc_1_sub = tc_1_all(side1);
    tc_2_sub = tc_2_all(side2);
    tc_t_sub = tc_t_all(sidet);

    tm_1_sub = tm_1_all(side1);
    tm_2_sub = tm_2_all(side2);
    tm_t_sub = tm_t_all(sidet);

    kc_1_sub = kc_1_all(side1);
    kc_2_sub = kc_2_all(side2);
    kc_t_sub = kc_t_all(sidet);

    km_1_sub = km_1_all(side1);
    km_2_sub = km_2_all(side2);
    km_t_sub = km_t_all(sidet);

    xq1_sub = xq1_all(side1);
    xq2_sub = xq2_all(side2);
    xqt_sub = xqt_all(sidet);

    yq1_sub = yq1_all(side1);
    yq2_sub = yq2_all(side2);
    yqt_sub = yqt_all(sidet);
    
    zq1_sub = zq1_all(side1);
    zq2_sub = zq2_all(side2);
    zqt_sub = zqt_all(sidet);

    normals1_sub = normals1_all(side1,:);
    normals2_sub = normals2_all(side2,:);
    normalst_sub = normalst_all(sidet,:);

    angles1_sub = angles1_all(side1);
    angles2_sub = angles2_all(side2);
    anglestip_sub = anglestip_all(sidet);
            % if IC(n,2)<=mid_y
                [~,N_before] = rotation_Q(tri_before,subset_before,n);
                
                in_tri = inpolygon(side1_rot(2,:),side1_rot(3,:),rotated_before(2,con_n),rotated_before(3,con_n)); 
                in_tri2 = inpolygon(side2_rot(2,:),side2_rot(3,:),rotated_before(2,con_n),rotated_before(3,con_n)); 
                in_tri_tip = inpolygon(sidet_rot(2,:),sidet_rot(3,:),rotated_before(2,con_n),rotated_before(3,con_n)); 
                
                % in_tri_tip = inpolygon(yqt_all,zqt_all,subset_before(2,con_n),subset_before(3,con_n));
                % curv_c = mean([mean(kc_1_all(in_tri),'omitnan') mean(kc_t_all(in_tri_tip),'omitnan') ],'omitnan');
                % curv_m = mean([mean(km_1_all(in_tri),'omitnan') mean(km_t_all(in_tri_tip),'omitnan') ],'omitnan');
                curv_c_dir = mean([mean(cdirs_1_sub(in_tri,:),1,'omitnan')',mean(cdirs_2_sub(in_tri2,:),1,'omitnan')' mean(cdirs_t_sub(in_tri_tip,:),1,'omitnan')'],2,'omitnan');
                curv_m_dir = mean([mean(mdirs_1_sub(in_tri,:),1,'omitnan')',mean(mdirs_2_sub(in_tri2,:),1,'omitnan')' mean(mdirs_t_sub(in_tri_tip,:),1,'omitnan')'],2,'omitnan');
            
                
                curv_1 = mean([mean(P1_side1_sub(in_tri),'omitnan'),mean(P1_side2_sub(in_tri2),'omitnan') mean(P1_tip_sub(in_tri_tip),'omitnan')],'omitnan');
                curv_2 = mean([mean(P2_side1_sub(in_tri),'omitnan'),mean(P2_side2_sub(in_tri2),'omitnan') mean(P2_tip_sub(in_tri_tip),'omitnan')],'omitnan');
                curv_1_dir = mean([mean(VP1_side1_sub(in_tri,:),1,'omitnan')',mean(VP1_side2_sub(in_tri2,:),1,'omitnan')' mean(VP1_tip_sub(in_tri_tip,:),1,'omitnan')'],2,'omitnan');
                curv_2_dir = mean([mean(VP2_side1_sub(in_tri,:),1,'omitnan')',mean(VP2_side2_sub(in_tri2,:),1,'omitnan')' mean(VP2_tip_sub(in_tri_tip,:),1,'omitnan')'],2,'omitnan');

                tens_c_comb = [tc_1_sub(in_tri);tc_2_sub(in_tri2);tc_t_sub(in_tri_tip)];
                tens_m_comb = [tm_1_sub(in_tri);tm_2_sub(in_tri2);tm_t_sub(in_tri_tip)];
                curv_c_comb = [kc_1_sub(in_tri);kc_2_sub(in_tri2);kc_t_sub(in_tri_tip)];
                curv_m_comb = [km_1_sub(in_tri);km_2_sub(in_tri2);km_t_sub(in_tri_tip)];
                xq_comb = [xq1_sub(in_tri);xq2_sub(in_tri2);xqt_sub(in_tri_tip)];
                yq_comb = [yq1_sub(in_tri);yq2_sub(in_tri2);yqt_sub(in_tri_tip)];
                zq_comb = [zq1_sub(in_tri);zq2_sub(in_tri2);zqt_sub(in_tri_tip)];
                n_sub = [normals1_sub(in_tri,:)',normals2_sub(in_tri2,:)',normalst_sub(in_tri_tip,:)'];
                angles_comb = [angles1_sub(in_tri);angles2_sub(in_tri2);anglestip_sub(in_tri_tip)];


                

            % else
            %     [~,N_before] = rotation_Q(tri_before,subset_before,n);
            %     in_tri = inpolygon(xq2_all,zq2_all,subset_before(1,con_n),subset_before(3,con_n));
            %     in_tri_tip = inpolygon(yqt_all,zqt_all,subset_before(2,con_n),subset_before(3,con_n));
            % 
            %     % curv_c = mean([mean(kc_2_all(in_tri),'omitnan') mean(kc_t_all(in_tri_tip),'omitnan')],'omitnan');
            %     % curv_m = mean([mean(km_2_all(in_tri),'omitnan') mean(km_t_all(in_tri_tip),'omitnan') ],'omitnan');
            %     curv_c_dir = mean([mean(cdirs_2_all(in_tri,:),1,'omitnan')' mean(cdirs_t_all(in_tri_tip,:),1,'omitnan')' ],2,'omitnan');
            %     curv_m_dir = mean([mean(mdirs_2_all(in_tri,:),1,'omitnan')' mean(mdirs_t_all(in_tri_tip,:),1,'omitnan')' ],2,'omitnan');
            % 
            %     curv_1 = mean([mean(P1_side2_all(in_tri),'omitnan') mean(P1_tip_all(in_tri_tip),'omitnan') ],'omitnan');
            %     curv_2 = mean([mean(P2_side2_all(in_tri),'omitnan') mean(P2_tip_all(in_tri_tip),'omitnan') ],'omitnan');
            %     curv_1_dir = mean([mean(VP1_side2_all(in_tri,:),1,'omitnan')' mean(VP1_tip_all(in_tri_tip,:),1,'omitnan')' ],2,'omitnan');
            %     curv_2_dir = mean([mean(VP2_side2_all(in_tri,:),1,'omitnan')' mean(VP2_tip_all(in_tri_tip,:),1,'omitnan')' ],2,'omitnan');
            % 
            %     tens_c_comb = [tc_2_all(in_tri);tc_t_all(in_tri_tip)];
            %     tens_m_comb = [tm_2_all(in_tri);tm_t_all(in_tri_tip)];
            %     curv_c_comb = [kc_2_all(in_tri);kc_t_all(in_tri_tip)];
            %     curv_m_comb = [km_2_all(in_tri);km_t_all(in_tri_tip)];
            %     angles_comb = [angles2_all(in_tri);anglestip_all(in_tri_tip)];
            % 
            % 
            %     xq_comb = [xq2_all(in_tri);xqt_all(in_tri_tip)];
            %     yq_comb = [yq2_all(in_tri);yqt_all(in_tri_tip)];
            %     zq_comb = [zq2_all(in_tri);zqt_all(in_tri_tip)];
            %     n_sub = [normals2_all(in_tri,:)',normalst_all(in_tri_tip,:)'];
                


    
            % end
                norm_comp = zeros(1,size(n_sub,2));
                for i = 1:size(n_sub,2)
                    norm_comp(:,i) = dot(N_before,n_sub(:,i));
                end
                norm_comp_all{n} = norm_comp;
                sub_norm = abs(norm_comp) >= 0;
                % size(find(sub_norm),2)/size(norm_comp,2);
                % if size(find(sub_norm),2) <=30 || ismember(n,tri_i) == 0
                %     tens_c = NaN;
                %     tens_m = NaN;
                %     curv_c = NaN;
                %     curv_m = NaN;
                %     % K = NaN;
                %     % MU = NaN;
                %     tens_c_rem = NaN(size(tens_c_comb(sub_norm)))';
                %     tens_m_rem = NaN(size(tens_c_comb(sub_norm)))';
                %     curv_c_rem = NaN(size(tens_c_comb(sub_norm)))';
                %     curv_m_rem = NaN(size(tens_c_comb(sub_norm)))';
                % else
                    % disp(n)
                    % tens_c_rem = pre_process(tens_c_comb',1,0,0,0);
                    % tens_m_rem = pre_process(tens_m_comb',1,0,0,0);
                    % curv_c_rem = pre_process(curv_c_comb',1,0,0,0);
                    % curv_m_rem = pre_process(curv_m_comb',1,0,0,0);

                    
                    tens_c_rem = (tens_c_comb)';
                    tens_m_rem = (tens_m_comb)';
                    curv_c_rem = (curv_c_comb)';
                    curv_m_rem = (curv_m_comb)';
                    % hold on; scatter3(xq_comb,yq_comb,zq_comb,10,tens_c_rem,'filled')
                    % text(IC(n,1),IC(n,2),IC(n,3),num2str(n))

                tens_c = mean(tens_c_rem,'omitnan');
                tens_m = mean(tens_m_rem,'omitnan');
                curv_c = mean(curv_c_rem,'omitnan');
                curv_m = mean(curv_m_rem,'omitnan');
                % tens_m = P ./ (2 * curv_c);
                % tens_c = P ./ (2 * curv_c) .* (2 - curv_m ./ curv_c);
               
                % end
                stretch_c = circ_stretch(n);
                stretch_m = mer_stretch(n);
                K = (tens_c_rem +tens_m_rem)./(2*(stretch_c*stretch_m-1));
                MU = (tens_m_rem-tens_c_rem)./((1/stretch_c^2)-(1/stretch_m^2));
                
                ref_tm = [ref_tm,tens_m_rem];
                ref_tc = [ref_tc,tens_c_rem];
                ref_km = [ref_km,curv_m_rem];
                ref_kc = [ref_kc,curv_c_rem];
                ref_K = [ref_K,K];
                ref_MU = [ref_MU,MU];
                angles_sub_all = [angles_sub_all,angles_comb'];
                xq_all = [xq_all,xq_comb'];
                yq_all = [yq_all,yq_comb'];
                zq_all = [zq_all,zq_comb'];
                % hold on; scatter(angles_comb(sub_norm)',K)

% curv_1_final = mean(curv_1,1,'omitnan');
% curv_2_final = mean(curv_2,1,'omitnan');
% curv_1_final = curv_1;
% curv_2_final = curv_2;





        cdir = cdir_input(:,n);%/norm(cdir_final(:,n));
        sdir = sdir_input(:,n);%/norm(sdir_final(:,n));

        if isnan(cdir(1)) || isnan(sdir(1))
            m_dirs_local = [NaN,NaN,NaN];
            c_dirs_local = [NaN,NaN,NaN];
        else
            vec_s = tip3 - IC(n,:)';
            vec_s = vec_s/norm(vec_s);
            options = optimoptions('fsolve','TolX',1e-6,'Display','off');
            theta_c = fsolve(@(theta) dot((cos(theta)*cdir+sin(theta)*sdir),vec_s),0,options);
            theta_m = pi/2 + theta_c;
                    
            c_dirs_local = [cos(theta_c)*cdir(1)+sin(theta_c)*sdir(1),cos(theta_c)*cdir(2)+sin(theta_c)*sdir(2),cos(theta_c)*cdir(3)+sin(theta_c)*sdir(3)];
            m_dirs_local = [cos(theta_m)*cdir(1)+sin(theta_m)*sdir(1),cos(theta_m)*cdir(2)+sin(theta_m)*sdir(2),cos(theta_m)*cdir(3)+sin(theta_m)*sdir(3)];
        
        end

        curv_1_dir = curv_1_dir/norm(curv_1_dir);
        curv_2_dir = curv_2_dir/norm(curv_2_dir);
        [~,N_before] = rotation_Q(tri_before,subset_before,n);

        N_old = cross(curv_1_dir,curv_2_dir)/norm(cross(curv_1_dir,curv_2_dir));


        v_new = cross(N_old,N_before); c = dot(N_old,N_before);
        v_x_new = [0 -v_new(3) v_new(2); v_new(3) 0 -v_new(1); -v_new(2) v_new(1) 0];
        new_Q_b = eye(3) + v_x_new + ((v_x_new^2)*(1/(1+c)));

        pdir1 = new_Q_b*curv_1_dir;
        pdir2 = new_Q_b*curv_2_dir;
        pdir1 = pdir1/norm(pdir1);
        pdir2 = pdir2/norm(pdir2);

        % cm = dot(pdir1,m_dirs_local');
        % sm = dot(pdir2,m_dirs_local');
        % if isnan(pdir1(1))
        % else
        %     rank([pdir2,m_dirs_local',c_dirs_local'])
        % end
        % 
        curv_c_dir = curv_c_dir/norm(curv_c_dir);
        curv_m_dir = curv_m_dir/norm(curv_m_dir);
        [~,N_before] = rotation_Q(tri_before,subset_before,n);

        N_old = cross(curv_c_dir,curv_m_dir)/norm(cross(curv_c_dir,curv_m_dir));

        % dot(N_old,N_before)
        v_new = cross(N_old,N_before); c = dot(N_old,N_before);
        v_x_new = [0 -v_new(3) v_new(2); v_new(3) 0 -v_new(1); -v_new(2) v_new(1) 0];
        new_Q_b = eye(3) + v_x_new + ((v_x_new^2)*(1/(1+c)));

        cdir_n = new_Q_b*curv_c_dir;
        mdir_n = new_Q_b*curv_m_dir;
        cdir_n = cdir_n/norm(cdir_n);
        mdir_n = mdir_n/norm(mdir_n);
        % % 
        cm = dot(cdir_n,m_dirs_local');
        sm = dot(mdir_n,m_dirs_local');
        % 
        % if abs(sm) > abs(cm)
        %     sint = norm(cross(pdir2,m_dirs_local));
        %     km = curv_2*(sm^2) + curv_1*(sint^2);
        %     kc = curv_2*(sint^2) + curv_1*(sm^2);
        % 
        % 
        % else
        %     sint = norm(cross(pdir1,m_dirs_local));
        %     km = curv_1*(cm^2) + curv_2*(sint^2);
        %     kc = curv_1*(sint^2) + curv_2*(cm^2);
        % end

        if abs(sm) > abs(cm)
            sint = norm(cross(mdir_n,m_dirs_local));
            km = curv_m*(sm^2) + curv_c*(sint^2);
            kc = curv_m*(sint^2) + curv_c*(sm^2);
        else
            % disp(n)
            sint = norm(cross(cdir_n,m_dirs_local));
            km = curv_c*(cm^2) + curv_m*(sint^2);
            kc = curv_c*(sint^2) + curv_m*(cm^2);
        end

         
        t1_side = P ./ (2 * curv_1);
        t2_side = P ./ (2 * curv_1) .* (2 - curv_2 ./ curv_1);
        tm_side = P ./ (2 * kc);
        tc_side = P ./ (2 * kc) .* (2 - km ./ kc);
        % end
        tm_side = P ./ (2 * curv_c);
        tc_side = P ./ (2 * curv_c) .* (2 - curv_m ./ curv_c);
        
        princp_directions(:,n) = [curv_c_dir;curv_m_dir;curv_1_dir;curv_2_dir;m_dirs_local';c_dirs_local';cdir_n;mdir_n];
        curvs_tensions(:,n) = [curv_1;curv_2;curv_c;curv_m;tm_side;tc_side;tens_m;tens_c];
        % sdir = Q_b*sdir;
        % cm = dot(cdir,m_dirs_local);
        % sm = dot(sdir,m_dirs_local);
        % 
        % 
        % 
        % princp_directions(:,n) = [cdir;sdir;m_dirs_local';c_dirs_local'];
        % 
        %     sint = norm(cross(sdir,m_dirs_local));
        %     km = curv_2_final*(sm^2) + curv_1_final*(sint^2);
        %     kc = curv_2_final*(sint^2) + curv_1_final*(sm^2);
        % 
        % t1_side = P ./ (2 * curv_1_final);
        % t2_side = P ./ (2 * curv_1_final) .* (2 - curv_2_final ./ curv_1_final);
        % tm_side = P ./ (2 * kc);
        % tc_side = P ./ (2 * kc) .* (2 - km ./ kc);
        % % end
        % 
        % curvs_tensions(:,n) = [curv_1_final;curv_2_final;kc;km;t1_side;t2_side;tm_side;tc_side];

        
        end

        % save('exp_test.mat')