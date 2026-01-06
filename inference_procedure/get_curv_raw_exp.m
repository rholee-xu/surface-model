function [angles1_all,angles2_all,P1_side1_all,P1_side2_all,anglestip_all,P1_tip_all,P2_tip_all,P2_side1_all,P2_side2_all,km_1_all,...
kc_1_all,km_2_all,kc_2_all,km_t_all,kc_t_all,xq1_all,yq1_all,zq1_all,xq2_all,yq2_all,zq2_all,xqt_all,yqt_all,zqt_all,...
VP1_side1_all,VP2_side1_all,VP1_side2_all,VP2_side2_all,VP1_tip_all,VP2_tip_all,mdirs_1_all,cdirs_1_all,mdirs_2_all,cdirs_2_all,mdirs_t_all,cdirs_t_all,...
normals1_all,normals2_all,normalst_all] = get_curv_raw_exp(before_xy_ridge_all,spacing,rad,subset_before,long_axis)
%%
P = 2;

% maxz = max(before_xy_ridge_all(3,:));
maxz = max(subset_before(3,:));
minz = min(subset_before(3,:));
[tip_point,i] = max(before_xy_ridge_all(1,:));
tip3 = before_xy_ridge_all(:,i);
mid_y = before_xy_ridge_all(2,i);

minx = min(subset_before(1,:));
% side1_ridge = before_xy_ridge_all(:,before_xy_ridge_all(2,:)<=70 & before_xy_ridge_all(1,:)<=(tip_point-0.5*rad));
% side2_ridge = before_xy_ridge_all(:,before_xy_ridge_all(2,:)>70 & before_xy_ridge_all(1,:)<=(tip_point-0.5*rad));
side1_ridge = before_xy_ridge_all(:,before_xy_ridge_all(2,:)<=mid_y & before_xy_ridge_all(1,:)<=(tip_point-0.3*rad)...
    & before_xy_ridge_all(3,:)~=0 & before_xy_ridge_all(1,:)>= minx-30 & before_xy_ridge_all(3,:)<=maxz+20  & before_xy_ridge_all(3,:)>=minz-20);
side2_ridge = before_xy_ridge_all(:,before_xy_ridge_all(2,:)>mid_y & before_xy_ridge_all(1,:)<=(tip_point-0.3*rad) ...
    & before_xy_ridge_all(3,:)~=0 & before_xy_ridge_all(1,:)>= minx-30 & before_xy_ridge_all(3,:)<=maxz+20  & before_xy_ridge_all(3,:)>=minz-20);
tip_ridge = before_xy_ridge_all(:,before_xy_ridge_all(1,:)>(tip_point-1.5*rad));%& before_xy_ridge_all(3,:)<=maxz & before_xy_ridge_all(3,:)>=minz);


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

normals1_all = [];
normals2_all = [];
normalst_all = [];

xq1_all = [];
yq1_all = [];
zq1_all = [];
xq2_all = [];
yq2_all = [];
zq2_all = [];
xqt_all = [];
yqt_all = [];
zqt_all = [];

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

for k = 0:spacing-1
    for m = 0:spacing-1
    [xq1,zq1]=meshgrid(min(side1_ridge(1,:))-k:spacing:max(side1_ridge(1,:)),min(side1_ridge(3,:))-m:spacing:max(side1_ridge(3,:)));
        yq1 = griddata(side1_ridge(1,:),side1_ridge(3,:),side1_ridge(2,:),xq1,zq1,'cubic');
    [xq2,zq2]=meshgrid(min(side2_ridge(1,:))-k:spacing:max(side2_ridge(1,:)),min(side2_ridge(3,:))-m:spacing:max(side2_ridge(3,:)));
        yq2 = griddata(side2_ridge(1,:),side2_ridge(3,:),side2_ridge(2,:),xq2,zq2,'cubic');
    [yqt,zqt]=meshgrid(min(tip_ridge(2,:))-k:spacing:max(tip_ridge(2,:)),min(tip_ridge(3,:))-m:spacing:max(tip_ridge(3,:)));
        xqt = griddata(tip_ridge(2,:),tip_ridge(3,:),tip_ridge(1,:),yqt,zqt,'cubic');

    %     [xq1,zq1]=meshgrid(min(side1_ridge(1,:))-1-k:spacing:max(side1_ridge(1,:))+1-k,min(side1_ridge(3,:))-1-m:spacing:max(side1_ridge(3,:))+1-m);
    %     yq1 = griddata(side1_ridge(1,:),side1_ridge(3,:),side1_ridge(2,:),xq1,zq1,'cubic');
    % [xq2,zq2]=meshgrid(min(side2_ridge(1,:))-1-k:spacing:max(side2_ridge(1,:))+1-k,min(side2_ridge(3,:))-1-m:spacing:max(side2_ridge(3,:))+1-m);
    %     yq2 = griddata(side2_ridge(1,:),side2_ridge(3,:),side2_ridge(2,:),xq2,zq2,'cubic');
    % [yqt,zqt]=meshgrid(min(tip_ridge(2,:))-1-k:spacing:max(tip_ridge(2,:))+1-k,min(tip_ridge(3,:))-1-m:spacing:max(tip_ridge(3,:))+1-m);
    %     xqt = griddata(tip_ridge(2,:),tip_ridge(3,:),tip_ridge(1,:),yqt,zqt,'cubic');
    
%     xq1_all = [xq1_all xq1];
%     scatter3(xq1,yq1,zq1,5,'filled');
% else
% 
% [xq1,zq1]=meshgrid(min(side1_beads(1,:)):spacing:max(side1_beads(1,:)),min(side1_beads(3,:)):spacing:max(side1_beads(3,:)));
%     yq1 = griddata(side1_beads(1,:),side1_beads(3,:),side1_beads(2,:),xq1,zq1,'cubic');
% [xq2,zq2]=meshgrid(min(side2_beads(1,:)):spacing:max(side2_beads(1,:)),min(side2_beads(3,:)):spacing:max(side2_beads(3,:)));
%     yq2 = griddata(side2_beads(1,:),side2_beads(3,:),side2_beads(2,:),xq2,zq2,'cubic');
% [yqt,zqt]=meshgrid(min(tip_beads(2,:)):spacing:max(tip_beads(2,:)),min(tip_beads(3,:)):spacing:max(tip_beads(3,:)));
%     xqt = griddata(tip_beads(2,:),tip_beads(3,:),tip_beads(1,:),yqt,zqt,'cubic');
% end

    X2D1_microns = xq1.*0.16;
    Y2D1_microns = yq1.*0.16;
    Z2D1_microns = zq1.*0.16;

    X2D2_microns = xq2.*0.16;
    Y2D2_microns = yq2.*0.16;
    Z2D2_microns = zq2.*0.16;
    X2Dt_microns = xqt.*0.16;
    Y2Dt_microns = yqt.*0.16;
    Z2Dt_microns = zqt.*0.16;
    
  

    [~,~,P1_side1,P2_side1,VP_max1_1,VP_max2_1,VP_max3_1,VP_min1_1,VP_min2_1,VP_min3_1,normals1] = surfature(X2D1_microns,Y2D1_microns,Z2D1_microns);
    % if size(X2D2_microns,1) <=1
    %     % disp(1)
    %     P1_side2 = zeros(size(P1_side1));
    %     P2_side2 = zeros(size(P1_side1));
    %     VP_max1_2 = zeros(size(VP_max1_1));
    %     VP_max2_2= zeros(size(VP_max1_1));
    %     VP_max3_2 = zeros(size(VP_max1_1));
    %     VP_min1_2 = zeros(size(VP_max1_1));
    %     VP_min2_2 = zeros(size(VP_max1_1));
    %     VP_min3_2 = zeros(size(VP_max1_1));
    %     normals2 = zeros(size(normals1));
    % else
    [~,~,P1_side2,P2_side2,VP_max1_2,VP_max2_2,VP_max3_2,VP_min1_2,VP_min2_2,VP_min3_2,normals2] = surfature(X2D2_microns,Y2D2_microns,Z2D2_microns);
    % end
    [~,~,P1_tip,P2_tip,VP_max1_t,VP_max2_t,VP_max3_t,VP_min1_t,VP_min2_t,VP_min3_t,normalstip] = surfature(X2Dt_microns,Y2Dt_microns,Z2Dt_microns);
  

        
        P1_side1 = -P1_side1(:);
        P2_side1 = -P2_side1(:);

        VP_max1_1 = VP_max1_1(:);
        VP_max2_1 = VP_max2_1(:);
        VP_max3_1 = VP_max3_1(:);
        VP_min1_1 = VP_min1_1(:);
        VP_min2_1 = VP_min2_1(:);
        VP_min3_1 = VP_min3_1(:);

        VP_max1_t = VP_max1_t(:);
        VP_max2_t = VP_max2_t(:);
        VP_max3_t = VP_max3_t(:);
        VP_min1_t = VP_min1_t(:);
        VP_min2_t = VP_min2_t(:);
        VP_min3_t = VP_min3_t(:);

        VP_max1_2 = VP_max1_2(:);
        VP_max2_2 = VP_max2_2(:);
        VP_max3_2 = VP_max3_2(:);
        VP_min1_2 = VP_min1_2(:);
        VP_min2_2 = VP_min2_2(:);
        VP_min3_2 = VP_min3_2(:);

        P1_side2 = P1_side2(:);
        P2_side2 = P2_side2(:);
        

        P1_tip = -P1_tip(:);
        P2_tip = -P2_tip(:);

        % km = zeros(size(P1_side1));
        % kc = zeros(size(P1_side1));
        angles1 = zeros(size(P1_side1));
        angles2 = zeros(size(P1_side2));
        anglestip = zeros(size(P1_tip));
        % 
        for i = 1:size(P1_side1,1)
            angles1(i) = acosd(dot(normals1(i,:)',long_axis));
        end
        for i = 1:size(P1_side2,1)
            angles2(i) = 180-acosd(dot(normals2(i,:)',long_axis));
        end
        for i = 1:size(P1_tip,1)
            anglestip(i) = acosd(dot(normalstip(i,:)',long_axis));
        end
     
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
            % 
            % m_dirs_local= (new_Q_b\[1;0;0]);
            % c_dirs_local= (new_Q_b\[0;1;0]);
            % hold on;quiver3(xq1(i),yq1(i),zq1(i),c_dirs_local(1),c_dirs_local(2),c_dirs_local(3),'r')
            % hold on;quiver3(xq1(i),yq1(i),zq1(i),VP_max1_1(i),VP_max2_1(i),VP_max3_1(i),'b')

            cdir = cdir_final/norm(cdir_final);
            sdir = sdir_final/norm(sdir_final);
            

            if isnan(cdir(1)) || isnan(sdir(1))
                m_dirs_local = [NaN,NaN,NaN];
                c_dirs_local = [NaN,NaN,NaN];
            else
                vec_s = tip3 - [xq1(i);yq1(i);zq1(i)];
                vec_s = vec_s/norm(vec_s);
                options = optimoptions('fsolve','TolX',1e-6,'Display','off');
                theta_c = fsolve(@(theta) dot((cos(theta)*cdir+sin(theta)*sdir),vec_s),0,options);
                theta_m = pi/2 + theta_c;
            
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

        km_2 = zeros(size(P1_side2));
        kc_2 = zeros(size(P2_side2));

        mdirs_2 = zeros(3,size(P1_side2,1));
        cdirs_2 = zeros(3,size(P2_side1,1));
        for i = 1:size(P1_side2,1)
            % N_before = normals2(i,:)';
            cdir_final = [VP_max1_2(i),VP_max2_2(i),VP_max3_2(i)];
            sdir_final = [VP_min1_2(i),VP_min2_2(i),VP_min3_2(i)];

            % v_new = cross(N_before,[0;0;1]); c = dot(N_before,[0;0;1]);
            % v_x_new = [0 -v_new(3) v_new(2); v_new(3) 0 -v_new(1); -v_new(2) v_new(1) 0];
            % new_Q_b = eye(3) + v_x_new + ((v_x_new^2)*(1/(1+c)));
            % 
            % m_dirs_local= (new_Q_b\[1;0;0]);
            % c_dirs_local= (new_Q_b\[0;1;0]);
            % hold on;quiver3(xq1(i),yq1(i),zq1(i),c_dirs_local(1),c_dirs_local(2),c_dirs_local(3),'r')
            % hold on;quiver3(xq1(i),yq1(i),zq1(i),VP_max1_1(i),VP_max2_1(i),VP_max3_1(i),'b')

            cdir = cdir_final/norm(cdir_final);
            sdir = sdir_final/norm(sdir_final);
            if isnan(cdir(1)) || isnan(sdir(1))
                m_dirs_local = [NaN,NaN,NaN];
                c_dirs_local = [NaN,NaN,NaN];
            else
                vec_s = tip3 - [xq2(i);yq2(i);zq2(i)];
                vec_s = vec_s/norm(vec_s);
                options = optimoptions('fsolve','TolX',1e-6,'Display','off');
                theta_c = fsolve(@(theta) dot((cos(theta)*cdir+sin(theta)*sdir),vec_s),0,options);
                theta_m = pi/2 + theta_c;
            
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
            % N_before = normalstip(i,:)';
            cdir_final = [VP_max1_t(i),VP_max2_t(i),VP_max3_t(i)];
            sdir_final = [VP_min1_t(i),VP_min2_t(i),VP_min3_t(i)];

            % v_new = cross(N_before,[0;0;1]); c = dot(N_before,[0;0;1]);
            % v_x_new = [0 -v_new(3) v_new(2); v_new(3) 0 -v_new(1); -v_new(2) v_new(1) 0];
            % new_Q_b = eye(3) + v_x_new + ((v_x_new^2)*(1/(1+c)));

            % m_dirs_local= (new_Q_b\[1;0;0]);
            % c_dirs_local= (new_Q_b\[0;1;0]);
            % hold on;quiver3(xq1(i),yq1(i),zq1(i),c_dirs_local(1),c_dirs_local(2),c_dirs_local(3),'r')
            % hold on;quiver3(xq1(i),yq1(i),zq1(i),VP_max1_1(i),VP_max2_1(i),VP_max3_1(i),'b')

            cdir = cdir_final/norm(cdir_final);
            sdir = sdir_final/norm(sdir_final);

            if isnan(cdir(1)) || isnan(sdir(1))
                m_dirs_local = [NaN,NaN,NaN];
                c_dirs_local = [NaN,NaN,NaN];
            else
                vec_s = tip3 - [xqt(i);yqt(i);zqt(i)];
                vec_s = vec_s/norm(vec_s);
                options = optimoptions('fsolve','TolX',1e-6,'Display','off');
                theta_c = fsolve(@(theta) dot((cos(theta)*cdir+sin(theta)*sdir),vec_s),0,options);
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
        normals1_all = [normals1_all;normals1];
        normals2_all = [normals2_all;normals2];
        normalst_all = [normalst_all;normalstip];
    end
end

end