function [Q_b,N_before] = rotation_Q(tri_before,subset_all_before,n)


    beads_before_t = subset_all_before(:,tri_before.ConnectivityList(n,:));

    e3 = [1;0;0];
 
    % before
    v1_before = beads_before_t(:,1) - beads_before_t(:,2);
    v2_before = beads_before_t(:,1) - beads_before_t(:,3);
    
    N_before = cross(v1_before,v2_before)/norm(cross(v1_before,v2_before));
    
%     full_before = [v1_before v2_before N_before];
    
    v_b = cross(N_before,e3); c = dot(N_before,e3);
    v_x_b = [0 -v_b(3) v_b(2); v_b(3) 0 -v_b(1); -v_b(2) v_b(1) 0];
    
    Q_b = eye(3) + v_x_b + ((v_x_b^2)*(1/(1+c)));





end