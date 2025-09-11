function [structAll_info] = get_curv_raw_sim(final128_r,spacing,noise)

[tip,it] = max(final128_r(1,:));
tip_point = final128_r(:,it);
R = max(final128_r(2,:));
% gof_all = zeros(1,size(tri_before.NonnectivityList,1));
% curvs_tensions = zeros(8,size(tri_before.NonnectivityList,1));
% princp_directions = zeros(12,size(tri_before.NonnectivityList,1));



% mid_y = final128_r(2,it);
% R = (max(final128_r(2,:))-min(final128_r(2,:)))/2;



maxz = max(final128_r(3,:));
side1_ridge = final128_r(:,final128_r(2,:)<=0 & final128_r(1,:)<=(tip-0.2*R) & final128_r(3,:)<=0.6*maxz & final128_r(3,:)~=0);
side2_ridge = final128_r(:,final128_r(2,:)>0 & final128_r(1,:)<=(tip-0.2*R) & final128_r(3,:)<=0.6*maxz & final128_r(3,:)~=0);
tip_ridge = final128_r(:,final128_r(1,:)>(tip-R)& final128_r(3,:)~=0& final128_r(2,:)~=0);
% quickscatter(tip_ridge)
top_ridge = final128_r(:,final128_r(3,:)>0.5*maxz & final128_r(1,:)<=(tip-0.2*R));


% k = 0;
% m = 0;
% P1_all = [];
% xq1_all = [];
% hold on;

P1_side1_all = [];
P1_side2_all = [];
P1_tip_all = [];
P2_side1_all = [];
P2_side2_all = [];
P2_tip_all = [];
angles1_all = [];
angles2_all = [];
anglestip_all = [];

km_1_all = [];
kc_1_all = [];
km_2_all = [];
kc_2_all = [];
km_t_all = [];
kc_t_all = [];

xq1_all = [];
yq1_all = [];
zq1_all = [];
xq2_all = [];
yq2_all = [];
zq2_all = [];
xqt_all = [];
yqt_all = [];
zqt_all = [];

normals1_all = [];
normals2_all = [];
normalst_all = [];

VP1_side1_all = [];
VP1_side2_all = [];
VP1_tip_all = [];
VP2_side1_all = [];
VP2_side2_all = [];
VP2_tip_all = [];

mdirs_1_all = [];
cdirs_1_all = [];
mdirs_2_all = [];
cdirs_2_all = [];
mdirs_t_all = [];
cdirs_t_all = [];
% save('surfature_variables_testing.mat','-append')
randnoise = rand(1);
for k = 0:0.025:spacing-0.025
    for m = 0:0.025:spacing-0.025
    [xq1,zq1]=meshgrid(min(side1_ridge(1,:))-k:spacing:max(side1_ridge(1,:)),min(side1_ridge(3,:))-m:spacing:max(side1_ridge(3,:)));
        yq1 = griddata(side1_ridge(1,:),side1_ridge(3,:),side1_ridge(2,:),xq1,zq1,'cubic');
    [xq2,zq2]=meshgrid(min(side2_ridge(1,:))-k:spacing:max(side2_ridge(1,:)),min(side2_ridge(3,:))-m:spacing:max(side2_ridge(3,:)));
        yq2 = griddata(side2_ridge(1,:),side2_ridge(3,:),side2_ridge(2,:),xq2,zq2,'cubic');
    [yqt,zqt]=meshgrid(min(tip_ridge(2,:))-k:spacing:max(tip_ridge(2,:)),min(tip_ridge(3,:))-m:spacing:max(tip_ridge(3,:)));
        xqt = griddata(tip_ridge(2,:),tip_ridge(3,:),tip_ridge(1,:),yqt,zqt,'cubic');
    [xqtop,yqtop]=meshgrid(min(top_ridge(1,:))-k:spacing:max(top_ridge(1,:)),min(top_ridge(2,:))-m:spacing:max(top_ridge(2,:)));
        zqtop = griddata(top_ridge(1,:),top_ridge(2,:),top_ridge(3,:),xqtop,yqtop,'cubic');

    %     [xq1,zq1]=meshgrid(min(side1_ridge(1,:))-0.025-k:spacing:max(side1_ridge(1,:))+0.025-k,min(side1_ridge(3,:))-0.025-m:spacing:max(side1_ridge(3,:))+0.025-m);
    %     yq1 = griddata(side1_ridge(1,:),side1_ridge(3,:),side1_ridge(2,:),xq1,zq1,'cubic');
    % [xq2,zq2]=meshgrid(min(side2_ridge(1,:))-0.025-k:spacing:max(side2_ridge(1,:))+0.025-k,min(side2_ridge(3,:))-0.025-m:spacing:max(side2_ridge(3,:))+0.025-m);
    %     yq2 = griddata(side2_ridge(1,:),side2_ridge(3,:),side2_ridge(2,:),xq2,zq2,'cubic');
    % [yqt,zqt]=meshgrid(min(tip_ridge(2,:))-0.025-k:spacing:max(tip_ridge(2,:))+0.025-k,min(tip_ridge(3,:))-0.025-m:spacing:max(tip_ridge(3,:))+0.025-m);
    %     xqt = griddata(tip_ridge(2,:),tip_ridge(3,:),tip_ridge(1,:),yqt,zqt,'cubic');
    % [xqtop,yqtop]=meshgrid(min(top_ridge(1,:))-0.025-k:spacing:max(top_ridge(1,:))+0.025-k,min(top_ridge(2,:))-0.025-m:spacing:max(top_ridge(2,:))+0.025-m);
    %     zqtop = griddata(top_ridge(1,:),top_ridge(2,:),top_ridge(3,:),xqtop,yqtop,'cubic');

        yqt((abs(yqt)<0.003 & abs(zqt)<0.03) | (abs(yqt)<0.03 & abs(zqt)<0.003)) = NaN;
% scatter3(xqt(t),yqt(t),zqt(t),'filled')
        
        
        
%     xq1_all = [xq1_all xq1];
    % hold on;scatter3(xqt,yqt,zqt,20,'k','filled');
% else
% 
% [xq1,zq1]=meshgrid(min(side1_beads(1,:)):spacing:max(side1_beads(1,:)),min(side1_beads(3,:)):spacing:max(side1_beads(3,:)));
%     yq1 = griddata(side1_beads(1,:),side1_beads(3,:),side1_beads(2,:),xq1,zq1,'cubic');
% [xq2,zq2]=meshgrid(min(side2_beads(1,:)):spacing:max(side2_beads(1,:)),min(side2_beads(3,:)):spacing:max(side2_beads(3,:)));
%     yq2 = griddata(side2_beads(1,:),side2_beads(3,:),side2_beads(2,:),xq2,zq2,'cubic');
% [yqt,zqt]=meshgrid(min(tip_beads(2,:)):spacing:max(tip_beads(2,:)),min(tip_beads(3,:)):spacing:max(tip_beads(3,:)));
%     xqt = griddata(tip_beads(2,:),tip_beads(3,:),tip_beads(1,:),yqt,zqt,'cubic');
% end
    [~,~,P1_side1,P2_side1,VP_max1_1,VP_max2_1,VP_max3_1,VP_min1_1,VP_min2_1,VP_min3_1,normals1] = surfature(xq1,yq1,zq1);
    [~,~,P1_side2,P2_side2,VP_max1_2,VP_max2_2,VP_max3_2,VP_min1_2,VP_min2_2,VP_min3_2,normals2] = surfature(xq2,yq2,zq2);
    [~,~,P1_tip,P2_tip,VP_max1_t,VP_max2_t,VP_max3_t,VP_min1_t,VP_min2_t,VP_min3_t,normalstip] = surfature(xqt,yqt,zqt);
    [~,~,P1_top,P2_top,VP_max1_p,VP_max2_p,VP_max3_p,VP_min1_p,VP_min2_p,VP_min3_p,normalstop] = surfature(xqtop,yqtop,zqtop);
    
        % P1_tip(abs(yqt)<0.003) = NaN;
        % P1_tip(abs(zqt)<0.003) = NaN;
        % P1_tip((abs(yqt)<0.003 & abs(zqt)<0.03) | (abs(yqt)<0.03 & abs(zqt)<0.003)) = NaN;
    
    if noise == 0
  
        
    
    else
        
        noise = -0.1 + (0.2)*randnoise;
        xq1_noise = xq1(:) + noise*normals1(:,1); xq1_noise = reshape(xq1_noise,size(xq1));
        yq1_noise = yq1(:) + noise*normals1(:,2); yq1_noise = reshape(yq1_noise,size(xq1));
        zq1_noise = zq1(:) + noise*normals1(:,3); zq1_noise = reshape(zq1_noise,size(xq1));
        xq2_noise = xq2(:) - noise*normals2(:,1); xq2_noise = reshape(xq2_noise,size(xq2));
        yq2_noise = yq2(:) - noise*normals2(:,2); yq2_noise = reshape(yq2_noise,size(xq2));
        zq2_noise = zq2(:) - noise*normals2(:,3); zq2_noise = reshape(zq2_noise,size(xq2));
        xqt_noise = xqt(:) + noise*normalstip(:,1);xqt_noise = reshape(xqt_noise,size(xqt));
        yqt_noise = yqt(:) + noise*normalstip(:,2);yqt_noise = reshape(yqt_noise,size(xqt));
        zqt_noise = zqt(:) + noise*normalstip(:,3);zqt_noise = reshape(zqt_noise,size(xqt));
        xqtop_noise = xqtop(:) + noise*normalstop(:,1);xqtop_noise = reshape(xqtop_noise,size(xqtop));
        yqtop_noise = yqtop(:) + noise*normalstop(:,2);yqtop_noise = reshape(yqtop_noise,size(xqtop));
        zqtop_noise = zqtop(:) + noise*normalstop(:,3);zqtop_noise = reshape(zqtop_noise,size(xqtop));

        [~,~,P1_side1,P2_side1,VP_max1_1,VP_max2_1,VP_max3_1,VP_min1_1,VP_min2_1,VP_min3_1,normals1] = surfature(xq1_noise,yq1_noise,zq1_noise);
        [~,~,P1_side2,P2_side2,VP_max1_2,VP_max2_2,VP_max3_2,VP_min1_2,VP_min2_2,VP_min3_2,normals2] = surfature(xq2_noise,yq2_noise,zq2_noise);
        [~,~,P1_tip,P2_tip,VP_max1_t,VP_max2_t,VP_max3_t,VP_min1_t,VP_min2_t,VP_min3_t,normalstip] = surfature(xqt_noise,yqt_noise,zqt_noise);
        [~,~,P1_top,P2_top,VP_max1_p,VP_max2_p,VP_max3_p,VP_min1_p,VP_min2_p,VP_min3_p,normalstop] = surfature(xqtop_noise,yqtop_noise,zqtop_noise);

    xq1 = xq1_noise;
    xq2 = xq2_noise;
    xqt = xqt_noise;
    % xqtop = xqtop_noise;
    yq1 = yq1_noise;
    yq2 = yq2_noise;
    yqt = yqt_noise;
    % yqtop = yqtop_noise;
    zq1 = zq1_noise;
    zq2 = zq2_noise;
    zqt = zqt_noise;
    % zqtop = zqtop_noise;

       

    end

    % after new surfature change, no more remove border within surfature.m
    % P1_side2(P1_side2==0)=NaN;
    % P2_side2(P2_side2==0)=NaN;
    % P1_side1(P1_side1==0)=NaN;
    % P2_side1(P2_side1==0)=NaN;
    % P1_tip(P1_tip==0)=NaN;
    % P2_tip(P2_tip==0)=NaN;
    % VP_max1_t(isnan(P1_tip))=NaN;VP_max2_t(isnan(P1_tip))=NaN;VP_max3_t(isnan(P1_tip))=NaN;VP_min1_t(isnan(P2_tip))=NaN;VP_min2_t(isnan(P2_tip))=NaN;VP_min3_t(isnan(P2_tip))=NaN;
    % VP_max1_2(isnan(P1_side2))=NaN;VP_max2_2(isnan(P1_side2))=NaN;VP_max3_2(isnan(P1_side2))=NaN;VP_min1_2(isnan(P2_side2))=NaN;VP_min2_2(isnan(P2_side2))=NaN;VP_min3_2(isnan(P2_side2))=NaN;
    % VP_max1_1(isnan(P1_side1))=NaN;VP_max2_1(isnan(P1_side1))=NaN;VP_max3_1(isnan(P1_side1))=NaN;VP_min1_1(isnan(P2_side1))=NaN;VP_min2_1(isnan(P2_side1))=NaN;VP_min3_1(isnan(P2_side1))=NaN;


  

         % P = 2;
        P1_side1 = -P1_side1(:);
        P2_side1 = -P2_side1(:);

        VP_max1_1 = VP_max1_1(:);
        VP_max2_1 = VP_max2_1(:);
        VP_max3_1 = VP_max3_1(:);
        VP_min1_1 = VP_min1_1(:);
        VP_min2_1 = VP_min2_1(:);
        VP_min3_1 = VP_min3_1(:);

        P1_side2 = P1_side2(:);
        P2_side2 = P2_side2(:);

        P1_tip = -P1_tip(:);
        P2_tip = -P2_tip(:);
        P1_top = -P1_top(:);
        P2_top = -P2_tip(:);

        angles1 = zeros(size(P1_side1));
        angles2 = zeros(size(P1_side2));
        anglestip = zeros(size(P1_tip));
        anglestop = zeros(size(P1_top));
        % 
        for i = 1:size(P1_side1,1)
            angles1(i) = acosd(dot(normals1(i,:)',[1;0;0]));
        end
        for i = 1:size(P1_side2,1)
            angles2(i) = 180-acosd(dot(normals2(i,:)',[1;0;0]));
        end
        for i = 1:size(P1_tip,1)
            anglestip(i) = acosd(dot(normalstip(i,:)',[1;0;0]));
        end
        for i = 1:size(P1_top,1)
            anglestop(i) = acosd(dot(normalstop(i,:)',[1;0;0]));
        end
        % 
        
        % hold on;scatter3(xqt(:),yqt(:),zqt(:),20,P1_tip,'filled');axis equal;colorbar
        % colormap('jet')
        km_1 = zeros(size(P1_side1));
        kc_1 = zeros(size(P2_side1));

        mdirs_1 = zeros(3,size(P1_side1,1));
        cdirs_1 = zeros(3,size(P2_side2,1));
        for i = 1:size(P1_side1,1)
            % N_before = normals1(i,:)';
            cdir_final = [VP_max1_1(i),VP_max2_1(i),VP_max3_1(i)];
            sdir_final = [VP_min1_1(i),VP_min2_1(i),VP_min3_1(i)];

            % v_new = cross(N_before,[0;0;1]); c = dot(N_before,[0;0;1]);
            % v_x_new = [0 -v_new(3) v_new(2); v_new(3) 0 -v_new(1); -v_new(2) v_new(1) 0];
            % new_Q_b = eye(3) + v_x_new + ((v_x_new^2)*(1/(1+c)));

            % m_dirs_local= (new_Q_b\[1;0;0]);
            % c_dirs_local= (new_Q_b\[0;1;0]);
            % hold on;quiver3(xq1(i),yq1(i),zq1(i),c_dirs_local(1),c_dirs_local(2),c_dirs_local(3),'r')
            % hold on;quiver3(xq1(i),yq1(i),zq1(i),VP_max1_1(i),VP_max2_1(i),VP_max3_1(i),'b')

            cdir = cdir_final/norm(cdir_final);
            sdir = sdir_final/norm(sdir_final);
            

            if isnan(cdir_final(1))
                m_dirs_local = [NaN,NaN,NaN];
                c_dirs_local = [NaN,NaN,NaN];
            else
            options = optimoptions('fsolve','TolX',1e-6,'Display','off');
            
            % if newDirS == 1
            vec_s = tip_point - [xq1(i);yq1(i);zq1(i)];
            vec_s = vec_s/norm(vec_s);
            theta_c = fsolve(@(theta) dot((cos(theta)*cdir+sin(theta)*sdir),vec_s),0,options);
            % else
            % % % 
            % theta_c = fsolve(@(theta) dot((cos(theta)*cdir+sin(theta)*sdir),[1;0;0]),0,options);
            % end
            theta_m = pi/2 + theta_c;
            % hold on;quiver3(xq1(i),yq1(i),zq1(i),cdir(1),cdir(2),cdir(3),'b')
            % hold on;quiver3(xq1(i),yq1(i),zq1(i),cos(thetas(n))*cdir(1)+sin(thetas(n))*sdir(1),cos(thetas(n))*cdir(2)+sin(thetas(n))*sdir(2),cos(thetas(n))*cdir(3)+sin(thetas(n))*sdir(3),'c')
            
            c_dirs_local = [cos(theta_c)*cdir(1)+sin(theta_c)*sdir(1),cos(theta_c)*cdir(2)+sin(theta_c)*sdir(2),cos(theta_c)*cdir(3)+sin(theta_c)*sdir(3)];
            m_dirs_local = [cos(theta_m)*cdir(1)+sin(theta_m)*sdir(1),cos(theta_m)*cdir(2)+sin(theta_m)*sdir(2),cos(theta_m)*cdir(3)+sin(theta_m)*sdir(3)];

            end

            cm = dot(cdir,m_dirs_local);
            sm = dot(sdir,m_dirs_local);

            if abs(sm) > abs(cm)
            sint = norm(cross(sdir,m_dirs_local));
            km_1(i) = P2_side1(i)*(sm^2) + P1_side1(i)*(sint^2);
            kc_1(i) = P2_side1(i)*(sint^2) + P1_side1(i)*(sm^2);

           
            else

            sint = norm(cross(cdir,m_dirs_local));
            km_1(i) = P1_side1(i)*(cm^2) + P2_side1(i)*(sint^2);
            kc_1(i) = P1_side1(i)*(sint^2) + P2_side1(i)*(cm^2);
            end

            mdirs_1(:,i) = m_dirs_local';
            cdirs_1(:,i) = c_dirs_local';
        end
        % scatter(angles1,kc_1,'b');%hold on; scatter(angles1,km,'r')

        km_2 = zeros(size(P1_side2));
        kc_2 = zeros(size(P2_side2));

        mdirs_2 = zeros(3,size(P1_side2,1));
        cdirs_2 = zeros(3,size(P2_side1,1));
        for i = 1:size(P1_side2,1)
            N_before = normals2(i,:)';
            cdir_final = [VP_max1_2(i),VP_max2_2(i),VP_max3_2(i)];
            sdir_final = [VP_min1_2(i),VP_min2_2(i),VP_min3_2(i)];

            v_new = cross(N_before,[0;0;1]); c = dot(N_before,[0;0;1]);
            v_x_new = [0 -v_new(3) v_new(2); v_new(3) 0 -v_new(1); -v_new(2) v_new(1) 0];
            new_Q_b = eye(3) + v_x_new + ((v_x_new^2)*(1/(1+c)));

            m_dirs_local= (new_Q_b\[1;0;0]);
            c_dirs_local= (new_Q_b\[0;1;0]);
            % hold on;quiver3(xq1(i),yq1(i),zq1(i),c_dirs_local(1),c_dirs_local(2),c_dirs_local(3),'r')
            % hold on;quiver3(xq1(i),yq1(i),zq1(i),VP_max1_1(i),VP_max2_1(i),VP_max3_1(i),'b')

            cdir = cdir_final/norm(cdir_final);
            sdir = sdir_final/norm(sdir_final);
            if isnan(cdir_final(1))
                m_dirs_local = [NaN,NaN,NaN];
                c_dirs_local = [NaN,NaN,NaN];
            else
            options = optimoptions('fsolve','TolX',1e-6,'Display','off');
            
            % if newDirS == 1
            vec_s = tip_point - [xq2(i);yq2(i);zq2(i)];
            vec_s = vec_s/norm(vec_s);
            theta_c = fsolve(@(theta) dot((cos(theta)*cdir+sin(theta)*sdir),vec_s),0,options);
            % else
            % theta_c = fsolve(@(theta) dot((cos(theta)*cdir+sin(theta)*sdir),[1;0;0]),0,options);
            % end
            theta_m = pi/2 + theta_c;
            % hold on;quiver3(xq1(i),yq1(i),zq1(i),cdir(1),cdir(2),cdir(3),'b')
            % hold on;quiver3(xq1(i),yq1(i),zq1(i),cos(thetas(n))*cdir(1)+sin(thetas(n))*sdir(1),cos(thetas(n))*cdir(2)+sin(thetas(n))*sdir(2),cos(thetas(n))*cdir(3)+sin(thetas(n))*sdir(3),'c')
            
            c_dirs_local = [cos(theta_c)*cdir(1)+sin(theta_c)*sdir(1),cos(theta_c)*cdir(2)+sin(theta_c)*sdir(2),cos(theta_c)*cdir(3)+sin(theta_c)*sdir(3)];
            m_dirs_local = [cos(theta_m)*cdir(1)+sin(theta_m)*sdir(1),cos(theta_m)*cdir(2)+sin(theta_m)*sdir(2),cos(theta_m)*cdir(3)+sin(theta_m)*sdir(3)];

            end
            cm = dot(cdir,m_dirs_local);
            sm = dot(sdir,m_dirs_local);




            if abs(sm) > abs(cm)
            sint = norm(cross(sdir,m_dirs_local));
            km_2(i) = P2_side2(i)*(sm^2) + P1_side2(i)*(sint^2);
            kc_2(i) = P2_side2(i)*(sint^2) + P1_side2(i)*(sm^2);
            else
            sint = norm(cross(cdir,m_dirs_local));
            km_2(i) = P1_side2(i)*(cm^2) + P2_side2(i)*(sint^2);
            kc_2(i) = P1_side2(i)*(sint^2) + P2_side2(i)*(cm^2);
            end

            mdirs_2(:,i) = m_dirs_local';
            cdirs_2(:,i) = c_dirs_local';
            

        end

        km_t = zeros(size(P1_tip));
        kc_t = zeros(size(P2_tip));

        mdirs_t = zeros(3,size(P1_tip,1));
        cdirs_t = zeros(3,size(P2_tip,1));
        for i = 1:size(P1_tip,1)
            cdir_final = [VP_max1_t(i),VP_max2_t(i),VP_max3_t(i)];
            sdir_final = [VP_min1_t(i),VP_min2_t(i),VP_min3_t(i)];

            cdir = cdir_final/norm(cdir_final);
            sdir = sdir_final/norm(sdir_final);
            if isnan(cdir_final(1))
                m_dirs_local = [NaN,NaN,NaN];
                c_dirs_local = [NaN,NaN,NaN];
            else
            options = optimoptions('fsolve','TolX',1e-6,'Display','off');
            
            % if newDirS == 1
            vec_s = tip_point - [xqt(i);yqt(i);zqt(i)];
            vec_s = vec_s/norm(vec_s);
            theta_c = fsolve(@(theta) dot((cos(theta)*cdir+sin(theta)*sdir),vec_s),0,options);
            % else
            % theta_c = fsolve(@(theta) dot((cos(theta)*cdir+sin(theta)*sdir),[1;0;0]),0,options);
            % end
            theta_m = pi/2 + theta_c;
            
            c_dirs_local = [cos(theta_c)*cdir(1)+sin(theta_c)*sdir(1),cos(theta_c)*cdir(2)+sin(theta_c)*sdir(2),cos(theta_c)*cdir(3)+sin(theta_c)*sdir(3)];
            m_dirs_local = [cos(theta_m)*cdir(1)+sin(theta_m)*sdir(1),cos(theta_m)*cdir(2)+sin(theta_m)*sdir(2),cos(theta_m)*cdir(3)+sin(theta_m)*sdir(3)];

            end
            cm = dot(cdir,m_dirs_local);
            sm = dot(sdir,m_dirs_local);


            if abs(sm) > abs(cm)
            sint = norm(cross(sdir,m_dirs_local));
            km_t(i) = P2_tip(i)*(sm^2) + P1_tip(i)*(sint^2);
            kc_t(i) = P2_tip(i)*(sint^2) + P1_tip(i)*(sm^2);
           
            else

            sint = norm(cross(cdir,m_dirs_local));
            km_t(i) = P1_tip(i)*(cm^2) + P2_tip(i)*(sint^2);
            kc_t(i) = P1_tip(i)*(sint^2) + P2_tip(i)*(cm^2);
            end

            mdirs_t(:,i) = m_dirs_local';
            cdirs_t(:,i) = c_dirs_local';
        end
        
        % km_t = zeros(size(P1_tip));
        % kc_t = zeros(size(P2_tip));
        % 
        % mdirs_t = zeros(3,size(P1_tip,1));
        % cdirs_t = zeros(3,size(P2_tip,1));
        % for i = 1:size(P1_tip,1)
        %     cdir_final = [VP_max1_t(i),VP_max2_t(i),VP_max3_t(i)];
        %     sdir_final = [VP_min1_t(i),VP_min2_t(i),VP_min3_t(i)];
        % 
        %     cdir = cdir_final/norm(cdir_final);
        %     sdir = sdir_final/norm(sdir_final);
        %     if isnan(cdir_final(1))
        %         m_dirs_local = [NaN,NaN,NaN];
        %         c_dirs_local = [NaN,NaN,NaN];
        %     else
        %     options = optimoptions('fsolve','TolX',1e-6,'Display','off');
        %     theta_c = fsolve(@(theta) dot((cos(theta)*cdir+sin(theta)*sdir),[1;0;0]),0,options);
        %     theta_m = pi/2 + theta_c;
        % 
        %     c_dirs_local = [cos(theta_c)*cdir(1)+sin(theta_c)*sdir(1),cos(theta_c)*cdir(2)+sin(theta_c)*sdir(2),cos(theta_c)*cdir(3)+sin(theta_c)*sdir(3)];
        %     m_dirs_local = [cos(theta_m)*cdir(1)+sin(theta_m)*sdir(1),cos(theta_m)*cdir(2)+sin(theta_m)*sdir(2),cos(theta_m)*cdir(3)+sin(theta_m)*sdir(3)];
        % 
        %     end
        %     cm = dot(cdir,m_dirs_local);
        %     sm = dot(sdir,m_dirs_local);
        % 
        % 
        %     if abs(sm) > abs(cm)
        %     sint = norm(cross(sdir,m_dirs_local));
        %     km_t(i) = P2_tip(i)*(sm^2) + P1_tip(i)*(sint^2);
        %     kc_t(i) = P2_tip(i)*(sint^2) + P1_tip(i)*(sm^2);
        % 
        %     else
        % 
        %     sint = norm(cross(cdir,m_dirs_local));
        %     km_t(i) = P1_tip(i)*(cm^2) + P2_tip(i)*(sint^2);
        %     kc_t(i) = P1_tip(i)*(sint^2) + P2_tip(i)*(cm^2);
        %     end
        % 
        %     mdirs_t(:,i) = m_dirs_local';
        %     cdirs_t(:,i) = c_dirs_local';
        % end
        % hold on;scatter(anglestip,kc_1,'b');%hold on; scatter(anglestip,km,'r')

        xq1_all = [xq1_all;xq1(:)];
        yq1_all = [yq1_all;yq1(:)];
        zq1_all = [zq1_all;zq1(:)];

        xq2_all = [xq2_all;xq2(:)];
        yq2_all = [yq2_all;yq2(:)];
        zq2_all = [zq2_all;zq2(:)];

        xqt_all = [xqt_all;xqt(:)];
        yqt_all = [yqt_all;yqt(:)];
        zqt_all = [zqt_all;zqt(:)];

        km_1_all = [km_1_all;km_1];
        kc_1_all = [kc_1_all;kc_1];
        km_2_all = [km_2_all;km_2];
        kc_2_all = [kc_2_all;kc_2];
        km_t_all = [km_t_all;km_t];
        kc_t_all = [kc_t_all;kc_t];

        mdirs_1_all = [mdirs_1_all;mdirs_1'];
        cdirs_1_all = [cdirs_1_all;cdirs_1'];
        mdirs_2_all = [mdirs_2_all;mdirs_2'];
        cdirs_2_all = [cdirs_2_all;cdirs_2'];
        mdirs_t_all = [mdirs_t_all;mdirs_t'];
        cdirs_t_all = [cdirs_t_all;cdirs_t'];
        
        P1_side1_all = [P1_side1_all;P1_side1];
        P1_side2_all = [P1_side2_all;P1_side2];
        P1_tip_all = [P1_tip_all;P1_tip];
        P2_side1_all = [P2_side1_all;P2_side1];
        P2_side2_all = [P2_side2_all;P2_side2];
        P2_tip_all = [P2_tip_all;P2_tip];
        angles1_all = [angles1_all;angles1];
        angles2_all = [angles2_all;angles2];
        anglestip_all = [anglestip_all;anglestip];

        VP1_side1_all = [VP1_side1_all;[VP_max1_1(:),VP_max2_1(:),VP_max3_1(:)]];
        VP1_side2_all = [VP1_side2_all;[VP_max1_2(:),VP_max2_2(:),VP_max3_2(:)]];
     
        VP1_tip_all = [VP1_tip_all;[VP_max1_t(:),VP_max2_t(:),VP_max3_t(:)]];
        VP2_side1_all = [VP2_side1_all;[VP_min1_1(:),VP_min2_1(:),VP_min3_1(:)]];
        VP2_side2_all = [VP2_side2_all;[VP_min1_2(:),VP_min2_2(:),VP_min3_2(:)]];
        VP2_tip_all = [VP2_tip_all;[VP_min1_t(:),VP_min2_t(:),VP_min3_t(:)]];
        % scatter3(xq1(:),yq1(:),zq1(:),10,P1_side1,'filled');axis equal;colorbar
        normals1_all = [normals1_all;normals1];
        normals2_all = [normals2_all;normals2];
        normalst_all = [normalst_all;normalstip];
        
    end
end

cell_info = {angles1_all,angles2_all,P1_side1_all,P1_side2_all,anglestip_all,P1_tip_all,P2_tip_all,P2_side1_all,P2_side2_all,km_1_all,...
kc_1_all,km_2_all,kc_2_all,km_t_all,kc_t_all,xq1_all,yq1_all,zq1_all,xq2_all,yq2_all,zq2_all,xqt_all,yqt_all,zqt_all,...
VP1_side1_all,VP2_side1_all,VP1_side2_all,VP2_side2_all,VP1_tip_all,VP2_tip_all,mdirs_1_all,cdirs_1_all,mdirs_2_all,cdirs_2_all,mdirs_t_all,cdirs_t_all,...
normals1_all,normals2_all,normalst_all};

fields = ["angles1_all","angles2_all","P1_side1_all","P1_side2_all","anglestip_all","P1_tip_all","P2_tip_all","P2_side1_all","P2_side2_all","km_1_all",...
"kc_1_all","km_2_all","kc_2_all","km_t_all","kc_t_all","xq1_all","yq1_all","zq1_all","xq2_all","yq2_all","zq2_all","xqt_all","yqt_all","zqt_all",...
"VP1_side1_all","VP2_side1_all","VP1_side2_all","VP2_side2_all","VP1_tip_all","VP2_tip_all","mdirs_1_all","cdirs_1_all","mdirs_2_all","cdirs_2_all","mdirs_t_all","cdirs_t_all",...
"normals1_all","normals2_all","normalst_all"];
structAll_info = cell2struct(cell_info',fields);

end