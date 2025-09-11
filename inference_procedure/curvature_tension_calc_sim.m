function [curvs_tensions,princp_directions,ref_km,ref_kc,ref_tm,ref_tc,ref_K,ref_MU,angles_sub_all,norm_comp_all,...
    xq_all,yq_all,zq_all] = curvature_tension_calc_sim(final128_r,tri_before,subset_before,IC,spacing,noise,sdir_input,...
    cdir_input,circ_stretch,mer_stretch,tri_i,struct)

[~,it] = max(final128_r(1,:));
tip_point = final128_r(:,it);

curvs_tensions = zeros(8,size(tri_before.ConnectivityList,1));
princp_directions = zeros(24,size(tri_before.ConnectivityList,1));
% km = zeros(1,size(tri_before.ConnectivityList,1));
% kc = zeros(1,size(tri_before.ConnectivityList,1));
P = 2;
% count = 1;
% rank_x = [];

long_axis = [1;0;0];   

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

% [struct.angles1_all,struct.angles2_all,struct.P1_side1_all,struct.P1_side2_all,struct.anglestip_all,struct.P1_tip_all,struct.P2_tip_all,...
%     struct.P2_side1_all,struct.P2_side2_all,struct.km_1_all,...
% struct.kc_1_all,struct.km_2_all,struct.kc_2_all,struct.km_t_all,struct.kc_t_all,struct.xq1_all,struct.yq1_all,struct.zq1_all,struct.xq2_all,...
% struct.yq2_all,struct.zq2_all,struct.xqt_all,struct.yqt_all,struct.zqt_all,...
% struct.VP1_side1_all,struct.VP2_side1_all,struct.VP1_side2_all,struct.VP2_side2_all,struct.VP1_tip_all,...
% struct.VP2_tip_all,struct.mdirs_1_all,struct.cdirs_1_all,struct.mdirs_2_all,struct.cdirs_2_all,struct.mdirs_t_all,struct.cdirs_t_all,...
% struct.normals1_all,struct.normals2_all,struct.normalst_all] = get_curv_raw_sim(final128_r,spacing,noise,newDirS);

    % quiver3(xqt_all,yqt_all,zqt_all,mdirs_t_all(:,1),mdirs_t_all(:,2),mdirs_t_all(:,3),3,'r');axis equal
    % hold on; quiver3(xqt_all,yqt_all,zqt_all,cdirs_t_all(:,1),cdirs_t_all(:,2),cdirs_t_all(:,3),3,'b');axis equal
    
tm_1_all = P ./ (2 * struct.kc_1_all);
tc_1_all = P ./ (2 * struct.kc_1_all) .* (2 - struct.km_1_all ./ struct.kc_1_all);
tm_2_all = P ./ (2 * struct.kc_2_all);
tc_2_all = P ./ (2 * struct.kc_2_all) .* (2 - struct.km_2_all ./ struct.kc_2_all);
tm_t_all = P ./ (2 * struct.kc_t_all);
tc_t_all = P ./ (2 * struct.kc_t_all) .* (2 - struct.km_t_all ./ struct.kc_t_all);

for i = 1:size(struct.VP1_side1_all,1)
    struct.VP1_side1_all(i,:) = struct.VP1_side1_all(i,:)/norm(struct.VP1_side1_all(i,:));
    struct.VP2_side1_all(i,:) = struct.VP2_side1_all(i,:)/norm(struct.VP2_side1_all(i,:));
    % if VP1_side1_all(i,3) > 0
    %     VP1_side1_all(i,:) = -VP1_side1_all(i,:);
    % else
    % end
end
for i = 1:size(struct.VP1_side2_all,1)
    struct.VP1_side2_all(i,:) = struct.VP1_side2_all(i,:)/norm(struct.VP1_side2_all(i,:));
    struct.VP2_side2_all(i,:) = struct.VP2_side2_all(i,:)/norm(struct.VP2_side2_all(i,:));
    % if VP1_side2_all(i,3) > 0 
    %     VP1_side2_all(i,:) = -VP1_side2_all(i,:);
    % else
    % end
    % if VP2_side2_all(i,2) > 0 && VP2_side2_all(i,3) > 0
    %     VP2_side2_all(i,:) = -VP2_side2_all(i,:);
    % else
    % end
end
for i = 1:size(struct.VP1_tip_all,1)
    struct.VP1_tip_all(i,:) = struct.VP1_tip_all(i,:)/norm(struct.VP1_tip_all(i,:));
    struct.VP2_tip_all(i,:) = struct.VP2_tip_all(i,:)/norm(struct.VP2_tip_all(i,:));
    if struct.VP1_tip_all(i,3) > 0 
        struct.VP1_tip_all(i,:) = -struct.VP1_tip_all(i,:);
    else
    end
    % if VP1_tip_all(i,2) <0 && abs(VP1_tip_all(i,2)) >0.7
    %     VP1_tip_all(i,:) = -VP1_tip_all(i,:);
    % else
    % end

end
% scatter(angles2_all,km_1_all,'r');hold on; scatter(angles2_all,km_2_all,'r');hold on; scatter(anglestip_all,km_t_all,'r')

for n = 1:size(tri_before.ConnectivityList,1)
    con_n = tri_before.ConnectivityList(n,:);




            if IC(n,2)<=0
                [~,N_before] = rotation_Q(tri_before,subset_before,n);
                in_tri = inpolygon(struct.xq1_all,struct.zq1_all,subset_before(1,con_n),subset_before(3,con_n)); 
                in_tri_tip = inpolygon(struct.yqt_all,struct.zqt_all,subset_before(2,con_n),subset_before(3,con_n));
                % curv_c = mean([mean(kc_1_all(in_tri),'omitnan') mean(kc_t_all(in_tri_tip),'omitnan') ],'omitnan');
                % curv_m = mean([mean(km_1_all(in_tri),'omitnan') mean(km_t_all(in_tri_tip),'omitnan') ],'omitnan');
                curv_c_dir = mean([mean(struct.cdirs_1_all(in_tri,:),1,'omitnan')' mean(struct.cdirs_t_all(in_tri_tip,:),1,'omitnan')' ],2,'omitnan');
                curv_m_dir = mean([mean(struct.mdirs_1_all(in_tri,:),1,'omitnan')' mean(struct.mdirs_t_all(in_tri_tip,:),1,'omitnan')' ],2,'omitnan');
            
                
                curv_1 = mean([mean(struct.P1_side1_all(in_tri),'omitnan') mean(struct.P1_tip_all(in_tri_tip),'omitnan') ],'omitnan');
                curv_2 = mean([mean(struct.P2_side1_all(in_tri),'omitnan') mean(struct.P2_tip_all(in_tri_tip),'omitnan') ],'omitnan');
                % 
                % tens_c = mean([mean(tc_1_all(in_tri),'omitnan') mean(tc_t_all(in_tri_tip),'omitnan') ],'omitnan');
                % tens_m = mean([mean(tm_1_all(in_tri),'omitnan') mean(tm_t_all(in_tri_tip),'omitnan') ],'omitnan');
                curv_1_dir = mean([mean(struct.VP1_side1_all(in_tri,:),1,'omitnan')' mean(struct.VP1_tip_all(in_tri_tip,:),1,'omitnan')' ],2,'omitnan');
                curv_2_dir = mean([mean(struct.VP2_side1_all(in_tri,:),1,'omitnan')' mean(struct.VP2_tip_all(in_tri_tip,:),1,'omitnan')' ],2,'omitnan');

                % curv1_comb = [VP1_side1_all(in_tri,:);VP1_tip_all(in_tri_tip,:)];
                % curv1_comb = curv1_comb(~isnan(curv1_comb(:,1)),:);
                % curvc_comb = [cdirs_1_all(in_tri,:);cdirs_t_all(in_tri_tip,:)];
                % curvc_comb = curvc_comb(~isnan(curvc_comb(:,1)),:);
                % [data_signs,~,ic] = unique(sign(curv1_comb),'rows');
                % a_counts = accumarray(ic,1);
                % [~,max_a]=max(a_counts);
                % % value_counts = [data_signs, a_counts]
                % 
                % curv_1_dir = (mean(abs(curv1_comb),1).*data_signs(max_a,:))';
                % curv_c_dir = (mean(abs(curvc_comb),1).*data_signs(max_a,:))';
                % 
                % curv2_comb = [VP2_side1_all(in_tri,:);VP2_tip_all(in_tri_tip,:)];
                % curv2_comb = curv2_comb(~isnan(curv2_comb(:,1)),:);
                % curvm_comb = [mdirs_1_all(in_tri,:);mdirs_t_all(in_tri_tip,:)];
                % curvm_comb = curvm_comb(~isnan(curvm_comb(:,1)),:);
                % [data_signs,~,ic] = unique(sign(curv2_comb),'rows');
                % a_counts = accumarray(ic,1);
                % [~,max_a]=max(a_counts);
                % curv_2_dir = (mean(abs(curv2_comb),1).*data_signs(max_a,:))';
                % curv_m_dir = (mean(abs(curvm_comb),1).*data_signs(max_a,:))';
                % if isempty(data_signs)
                %     curv_1_dir = [NaN;NaN;NaN];
                %     curv_2_dir = [NaN;NaN;NaN];
                %     curv_c_dir = [NaN;NaN;NaN];
                %     curv_m_dir = [NaN;NaN;NaN];
                % else
                % end

               tens_c_comb = [tc_1_all(in_tri);tc_t_all(in_tri_tip)];
                tens_m_comb = [tm_1_all(in_tri);tm_t_all(in_tri_tip)];
                curv_c_comb = [struct.kc_1_all(in_tri);struct.kc_t_all(in_tri_tip)];
                curv_m_comb = [struct.km_1_all(in_tri);struct.km_t_all(in_tri_tip)];
                xq_comb = [struct.xq1_all(in_tri);struct.xqt_all(in_tri_tip)];
                yq_comb = [struct.yq1_all(in_tri);struct.yqt_all(in_tri_tip)];
                zq_comb = [struct.zq1_all(in_tri);struct.zqt_all(in_tri_tip)];
                n_sub = [struct.normals1_all(in_tri,:)',struct.normalst_all(in_tri_tip,:)'];
                angles_comb = [struct.angles1_all(in_tri);struct.anglestip_all(in_tri_tip)];
                % if size(find(sub_norm),2) <= 10
                %     hold on;trisurf(tri_before.ConnectivityList(n,:),subset_before(1,:),subset_before(2,:),subset_before(3,:),'FaceAlpha',0)
                % else
                % hold on; scatter3(xq_comb(sub_norm),yq_comb(sub_norm),zq_comb(sub_norm),10,'filled')
                % hold on;trisurf(tri_before.ConnectivityList(n,:),subset_before(1,:),subset_before(2,:),subset_before(3,:),'FaceAlpha',0.3)
                % hold on; text(IC(n,1),IC(n,2),IC(n,3),num2str(n));
                % end
            else
                [~,N_before] = rotation_Q(tri_before,subset_before,n);
                in_tri = inpolygon(struct.xq2_all,struct.zq2_all,subset_before(1,con_n),subset_before(3,con_n));
                in_tri_tip = inpolygon(struct.yqt_all,struct.zqt_all,subset_before(2,con_n),subset_before(3,con_n));
               
                curv_c = mean([mean(struct.kc_2_all(in_tri),'omitnan') mean(struct.kc_t_all(in_tri_tip),'omitnan')],'omitnan');
                curv_m = mean([mean(struct.km_2_all(in_tri),'omitnan') mean(struct.km_t_all(in_tri_tip),'omitnan') ],'omitnan');
                curv_c_dir = mean([mean(struct.cdirs_2_all(in_tri,:),1,'omitnan')' mean(struct.cdirs_t_all(in_tri_tip,:),1,'omitnan')' ],2,'omitnan');
                curv_m_dir = mean([mean(struct.mdirs_2_all(in_tri,:),1,'omitnan')' mean(struct.mdirs_t_all(in_tri_tip,:),1,'omitnan')' ],2,'omitnan');

                curv_1 = mean([mean(struct.P1_side2_all(in_tri),'omitnan') mean(struct.P1_tip_all(in_tri_tip),'omitnan') ],'omitnan');
                curv_2 = mean([mean(struct.P2_side2_all(in_tri),'omitnan') mean(struct.P2_tip_all(in_tri_tip),'omitnan') ],'omitnan');
                curv_1_dir = mean([mean(struct.VP1_side2_all(in_tri,:),1,'omitnan')' mean(struct.VP1_tip_all(in_tri_tip,:),1,'omitnan')' ],2,'omitnan');
                curv_2_dir = mean([mean(struct.VP2_side2_all(in_tri,:),1,'omitnan')' mean(struct.VP2_tip_all(in_tri_tip,:),1,'omitnan')' ],2,'omitnan');

             

                tens_c_comb = [tc_2_all(in_tri);tc_t_all(in_tri_tip)];
                tens_m_comb = [tm_2_all(in_tri);tm_t_all(in_tri_tip)];
                curv_c_comb = [struct.kc_2_all(in_tri);struct.kc_t_all(in_tri_tip)];
                curv_m_comb = [struct.km_2_all(in_tri);struct.km_t_all(in_tri_tip)];
                angles_comb = [struct.angles2_all(in_tri);struct.anglestip_all(in_tri_tip)];

                
                xq_comb = [struct.xq2_all(in_tri);struct.xqt_all(in_tri_tip)];
                yq_comb = [struct.yq2_all(in_tri);struct.yqt_all(in_tri_tip)];
                zq_comb = [struct.zq2_all(in_tri);struct.zqt_all(in_tri_tip)];
                n_sub = [struct.normals2_all(in_tri,:)',struct.normalst_all(in_tri_tip,:)'];
                % if size(find(sub_norm),2) <= 30
                %     hold on;trisurf(tri_before.ConnectivityList(n,:),subset_before(1,:),subset_before(2,:),subset_before(3,:),'FaceAlpha',0)
                % else
                % hold on; scatter3(xq_comb(sub_norm),yq_comb(sub_norm),zq_comb(sub_norm),10,'filled')
                % hold on;trisurf(tri_before.ConnectivityList(n,:),subset_before(1,:),subset_before(2,:),subset_before(3,:),'FaceAlpha',0.3)
                % hold on; text(IC(n,1),IC(n,2),IC(n,3),num2str(n));
                % end

                % curv1_comb = [VP1_side2_all(in_tri,:);VP1_tip_all(in_tri_tip,:)];
                % curv1_comb = curv1_comb(~isnan(curv1_comb(:,1)),:);
                % curvc_comb = [cdirs_2_all(in_tri,:);cdirs_t_all(in_tri_tip,:)];
                % curvc_comb = curvc_comb(~isnan(curvc_comb(:,1)),:);
                % [data_signs,~,ic] = unique(sign(curv1_comb),'rows');
                % a_counts = accumarray(ic,1);
                % [~,max_a]=max(a_counts);
                % value_counts = [data_signs, a_counts]
                
                % curv_1_dir = (mean(abs(curv1_comb),1).*data_signs(max_a,:))';
                % curv_c_dir = (mean(abs(curvc_comb),1).*data_signs(max_a,:))';
                % 
                % curv2_comb = [VP2_side2_all(in_tri,:);VP2_tip_all(in_tri_tip,:)];
                % curv2_comb = curv2_comb(~isnan(curv2_comb(:,1)),:);
                % curvm_comb = [mdirs_2_all(in_tri,:);mdirs_t_all(in_tri_tip,:)];
                % curvm_comb = curvm_comb(~isnan(curvm_comb(:,1)),:);
                % [data_signs,~,ic] = unique(sign(curv2_comb),'rows');
                % a_counts = accumarray(ic,1);
                % [~,max_a]=max(a_counts);
                % curv_2_dir = (mean(abs(curv2_comb),1).*data_signs(max_a,:))';
                % curv_m_dir = (mean(abs(curvm_comb),1).*data_signs(max_a,:))';
                % if isempty(data_signs)
                %     curv_1_dir = [NaN;NaN;NaN];
                %     curv_2_dir = [NaN;NaN;NaN];
                %     curv_c_dir = [NaN;NaN;NaN];
                %     curv_m_dir = [NaN;NaN;NaN];
                % else
                % end
                % 
                %
               
            end
                norm_comp = zeros(1,size(n_sub,2));
                for i = 1:size(n_sub,2)
                    norm_comp(:,i) = dot(N_before,n_sub(:,i));
                end
                norm_comp_all{n} = norm_comp;
                sub_norm = abs(norm_comp) >= 0;
                % size(find(sub_norm),2)/size(norm_comp,2);
                if size(find(sub_norm),2) <=5 || ismember(n,tri_i) == 0
                    tens_c = NaN;
                    tens_m = NaN;
                    curv_c = NaN;
                    curv_m = NaN;
                    % K = NaN;
                    % MU = NaN;
                    tens_c_rem = NaN(size(tens_c_comb(sub_norm)))';
                    tens_m_rem = NaN(size(tens_c_comb(sub_norm)))';
                    curv_c_rem = NaN(size(tens_c_comb(sub_norm)))';
                    curv_m_rem = NaN(size(tens_c_comb(sub_norm)))';
                else
                    % disp(n)
                    tens_c_rem = pre_process(abs(tens_c_comb(sub_norm))',1,0,0,0);
                    tens_m_rem = pre_process(abs(tens_m_comb(sub_norm))',1,0,0,0);
                    curv_c_rem = pre_process(abs(curv_c_comb(sub_norm))',1,0,0,0);
                    curv_m_rem = pre_process(abs(curv_m_comb(sub_norm))',1,0,0,0);

                tens_c = mean(tens_c_rem,'omitnan');
                tens_m = mean(tens_m_rem,'omitnan');
                curv_c = mean(curv_c_rem,'omitnan');
                curv_m = mean(curv_m_rem,'omitnan');
                % tens_m = P ./ (2 * curv_c);
                % tens_c = P ./ (2 * curv_c) .* (2 - curv_m ./ curv_c);
               
                end
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
                angles_sub_all = [angles_sub_all,angles_comb(sub_norm)'];
                xq_all = [xq_all,xq_comb(sub_norm)'];
                yq_all = [yq_all,yq_comb(sub_norm)'];
                zq_all = [zq_all,zq_comb(sub_norm)'];


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
            options = optimoptions('fsolve','TolX',1e-6,'Display','off');

            % if newDirS == 1
            vec_s = tip_point - IC(n,:)';
            vec_s = vec_s/norm(vec_s);
            theta_c = fsolve(@(theta) dot((cos(theta)*cdir+sin(theta)*sdir),vec_s),0,options);
            % else
            % theta_c = fsolve(@(theta) dot((cos(theta)*cdir+sin(theta)*sdir),long_axis),0,options);
            % end
            theta_m = pi/2 + theta_c;
                    
            c_dirs_local = [cos(theta_c)*cdir(1)+sin(theta_c)*sdir(1),cos(theta_c)*cdir(2)+sin(theta_c)*sdir(2),cos(theta_c)*cdir(3)+sin(theta_c)*sdir(3)];
            m_dirs_local = [cos(theta_m)*cdir(1)+sin(theta_m)*sdir(1),cos(theta_m)*cdir(2)+sin(theta_m)*sdir(2),cos(theta_m)*cdir(3)+sin(theta_m)*sdir(3)];
        
        end
% disp(n)
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
        
        princp_directions(:,n) = [curv_c_dir;curv_m_dir;curv_1_dir;curv_2_dir;m_dirs_local';c_dirs_local';cdir_n;mdir_n];
        curvs_tensions(:,n) = [curv_c;curv_m;curv_c;curv_m;t1_side;t2_side;tens_m;tens_c];

   

end