function [mer_stretch,circ_stretch,mdirs,cdirs,pdir1,pdir2,sdir_after_test,lambda1s,lambda2s,pdir_after_test] = stretch_directions(tri_before,subset_before,subset_after,...
    long_axis,finalpb_r)

% tip_ones = [];
% [mer_stretch,circ_stretch,mdirs,cdirs,pdir1,pdir2] = stretch_directions_unturgid(tri_before,subset_before,subset_after,long_axis,finalpb_r)
face_normals = tri_before.faceNormal;
[~,it] = max(finalpb_r(1,:));
tip_point = finalpb_r(:,it);
mer_stretch = zeros(1,size(tri_before.ConnectivityList,1));
circ_stretch = zeros(1,size(tri_before.ConnectivityList,1));
pdir1 = zeros(3,size(tri_before.ConnectivityList,1));
pdir2 = zeros(3,size(tri_before.ConnectivityList,1));
mdirs = zeros(3,size(tri_before.ConnectivityList,1));
cdirs = zeros(3,size(tri_before.ConnectivityList,1));
inc = incenter(tri_before);
for n = 1:size(tri_before.ConnectivityList,1)
    beads_after_t = subset_after(:,tri_before.ConnectivityList(n,:));
    beads_before_t = subset_before(:,tri_before.ConnectivityList(n,:));

    e3 = [0;0;1];
    % after
    v1_after = beads_after_t(:,1) - beads_after_t(:,2);
    v2_after = beads_after_t(:,1) - beads_after_t(:,3);


    N_after = cross(v1_after,v2_after)/norm(cross(v1_after,v2_after));
    
    full_after = [v1_after v2_after N_after];
    
    % v_a = cross(N_after,e3); c = dot(N_after,e3);
    % v_x = [0 -v_a(3) v_a(2); v_a(3) 0 -v_a(1); -v_a(2) v_a(1) 0];
    % 
    % Q_a = eye(3) + v_x + ((v_x^2)*(1/(1+c)));
    % before
    v1_before = beads_before_t(:,1) - beads_before_t(:,2);
    v2_before = beads_before_t(:,1) - beads_before_t(:,3);

    
    N_before = cross(v1_before,v2_before)/norm(cross(v1_before,v2_before));
    
    full_before = [v1_before v2_before N_before];
    % 
    % v_b = cross(N_before,e3); c = dot(N_before,e3);
    % v_x_b = [0 -v_b(3) v_b(2); v_b(3) 0 -v_b(1); -v_b(2) v_b(1) 0];
    % 
    % Q_b = eye(3) + v_x_b + ((v_x_b^2)*(1/(1+c)));
    
    %
    % after_2D = Q_a * full_after;
    % before_2D = Q_b * full_before;
    
        % F = before_2D(2:3,1:2)/after_2D(2:3,1:2);
        F = full_before/full_after;

        F_t = transpose(F);
        [~,D] = eig(F_t*F);
        [V,~] = eig(F*F_t);
        % [V,~] = eig(F_t*F);
        D = sqrt(D);
        % eigenvs(:,n) = [D(1,1);D(2,2)];
        % V_3D = [0 0;V(1,:);V(2,:)];
        % V_3D = Q_b\V_3D;
        V_3D = V;
        
        check = dot(V_3D,[face_normals(n,:)',face_normals(n,:)',face_normals(n,:)']);
        [~,m]=max(abs(check));
        if  m == 2
            % disp(D)
            cdir = V_3D(:,3);
            sdir = V_3D(:,1);
            lambda1 = D(3,3);
            lambda2 = D(1,1);
        elseif m == 1
            cdir = V_3D(:,3);
            sdir = V_3D(:,2);
            lambda1 = D(3,3);
            lambda2 = D(2,2);
        elseif m == 3
            cdir = V_3D(:,2);
            sdir = V_3D(:,1);
            lambda1 = D(2,2);
            lambda2 = D(1,1);
        end
 

    if isnan(cdir(1)) || isnan(sdir(1))
            mdirs(:,n) = [NaN,NaN,NaN];
            cdirs(:,n) = [NaN,NaN,NaN];
    else
    options = optimoptions('fsolve','TolX',1e-6,'Display','off');
    
    % if newDirS == 1
    vec_s = tip_point - inc(n,:)';
    vec_s = vec_s/norm(vec_s);
    theta_c = fsolve(@(theta) dot((cos(theta)*cdir+sin(theta)*sdir),vec_s),0,options);
    % else
    % theta_c = fsolve(@(theta) dot((cos(theta)*cdir+sin(theta)*sdir),long_axis),0,options);
    % end
    theta_m = pi/2 + theta_c;
    % if theta_c == 0
    %     tip_ones = [tip_ones,n];
    % end

    cdirs(:,n)= [cos(theta_c)*cdir(1)+sin(theta_c)*sdir(1),cos(theta_c)*cdir(2)+sin(theta_c)*sdir(2),cos(theta_c)*cdir(3)+sin(theta_c)*sdir(3)];
    mdirs(:,n) = [cos(theta_m)*cdir(1)+sin(theta_m)*sdir(1),cos(theta_m)*cdir(2)+sin(theta_m)*sdir(2),cos(theta_m)*cdir(3)+sin(theta_m)*sdir(3)];

    % hold on; quiver3(IC(n,1),IC(n,2),IC(n,3),cdir(1),cdir(2),cdir(3),'b')
    % temp = [cdir,sdir];
    % if rank(temp) == 2
    % else
    %     disp(rank(temp))
    % end
    end

    
 

    % maxis_z = Q_b*mdirs(:,n);
    % caxis_z = Q_b*cdirs(:,n);


    pdir1(:,n) = cdir;
    pdir2(:,n) = sdir;
    pdir_after_test(:,n) = (inv(F)*pdir1(:,n))/norm(inv(F)*pdir1(:,n));
    lambda1s(n) = lambda1;
    lambda2s(n) = lambda2;
    
        % mer_stretch(n) = 1/norm((inv(F))*maxis_z(2:3));
        % circ_stretch(n) = 1/norm((inv(F))*caxis_z(2:3));
        % 



    % F = full_before/full_after;

    % mer_stretch(n) = 1/norm(dot(sdir_unturgid(:,n),(inv(F))*mdirs(:,n)));
    % circ_stretch(n) = 1/norm(dot(cdir_unturgid(:,n),(inv(F))*cdirs(:,n)));% newS
    % mer_stretch(n) = norm(dot(mdirs(:,n),F*sdir_unturgid(:,n)));
    % circ_stretch(n) = norm(dot(cdirs(:,n),F*cdir_unturgid(:,n))); % newS2
    mer_stretch(n) = 1/norm((inv(F))*mdirs(:,n));
    circ_stretch(n) = 1/norm((inv(F))*cdirs(:,n));
    % 
    
    sdir_after_test(:,n) = (inv(F)*mdirs(:,n))/norm(inv(F)*mdirs(:,n));

    
    % mer_stretch(n) = 1/norm((sqrt(inv(F*F_t)))*maxis_z(1:2));
    % circ_stretch(n) = 1/norm((sqrt(inv(F*F_t)))*caxis_z(1:2));

    % mer_stretch(n) = 1/norm((inv(F))*maxis_z(1:2));
    % circ_stretch(n) = 1/norm((inv(F))*caxis_z(1:2));
 % save('temp_angle_mat')

end


end
