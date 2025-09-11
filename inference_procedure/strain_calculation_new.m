function [V_3D,D,triangle_area_after,triangle_area_before,err_max] = strain_calculation_new(beads_after_t,beads_before_t,del)



e3 = [0;0;1];
% after
v1_after = beads_after_t(:,1) - beads_after_t(:,2);
v2_after = beads_after_t(:,1) - beads_after_t(:,3);

N_after = cross(v1_after,v2_after)/norm(cross(v1_after,v2_after));

full_after = [v1_after v2_after N_after];

v_a = cross(N_after,e3); c = dot(N_after,e3);
v_x = [0 -v_a(3) v_a(2); v_a(3) 0 -v_a(1); -v_a(2) v_a(1) 0];

Q_a = eye(3) + v_x + ((v_x^2)*(1/(1+c)));
% before
v1_before = beads_before_t(:,1) - beads_before_t(:,2);
v2_before = beads_before_t(:,1) - beads_before_t(:,3);

N_before = cross(v1_before,v2_before)/norm(cross(v1_before,v2_before));

full_before = [v1_before v2_before N_before];

v_b = cross(N_before,e3); c = dot(N_before,e3);
v_x_b = [0 -v_b(3) v_b(2); v_b(3) 0 -v_b(1); -v_b(2) v_b(1) 0];

Q_b = eye(3) + v_x_b + ((v_x_b^2)*(1/(1+c)));

%
after_2D = Q_a * full_after;
before_2D = Q_b * full_before;

F = before_2D(1:2,1:2)/after_2D(1:2,1:2);

F_t = transpose(F);
[~,D] = eig(F_t*F);
[V,~] = eig(F*F_t);
D = sqrt(D);
V_3D = [V(1,:);V(2,:);0 0];
V_3D = Q_b\V_3D;

triangle_area_after = .5*norm(cross(v1_after,v2_after));
triangle_area_before = .5*norm(cross(v1_before,v2_before));

Mb = before_2D(1:2,1:2);
Ma = after_2D(1:2,1:2);
% F = (Ma/inv(Mb))/det(Mb);
b_e = sum(abs(Ma),'all');
%err_max = ((2*norm(Ma))/det(Mb) + (norm(F)/det(Mb))*(b_e)+2*norm(inv(Mb)))*del;
err_max = ((2*norm(Mb))/(det(Ma)*norm(F)) + (b_e/det(Ma))+((2*norm(inv(Ma)))/norm(F)))*del;
% err_max = (2*norm(inv(Mb))/norm(eye(2)))*del;
% err_max = ((2*norm(Ma))/(det(Mb)*norm(F)))*del;
% err_max = ((2*norm(inv(Mb))/norm(F)))*del;
% err_max = (b_e/det(Mb))*del;

end

