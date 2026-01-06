function [clp,clp_a,tri_ob_before,tri_ob_after,iz] = bead_projection(subset_before,subset_after,before_ridge,after_ridge,tri_before)

clp = zeros(size(subset_before));
dists = zeros(1,size(subset_before,2));
for n = 1:size(subset_before,2)
    bead = subset_before(:,n);
    ic = dsearchn(before_ridge',bead');
    clp(:,n) = before_ridge(:,ic);
    dists(n) = norm(bead - clp(:,n));
end


clp_a = zeros(size(subset_before));
dists_a = zeros(1,size(subset_before,2));
for n = 1:size(subset_before,2)
    bead = subset_after(:,n);
    ic = dsearchn(after_ridge',bead');
    clp_a(:,n) = after_ridge(:,ic);
    dists_a(n) = norm(bead - clp_a(:,n));
end

z_dists = subset_before(3,:) - clp(3,:);
x_dists = subset_before(1,:) - clp(1,:);

iz_b = abs(z_dists)<=5 & abs(x_dists)<=20;
z_dists_a = subset_after(3,:) - clp_a(3,:);
x_dists_a = subset_after(1,:) - clp_a(1,:);

iz_a = abs(z_dists_a)<=5 & abs(x_dists_a)<=20;

iz = iz_a & iz_b;
% if size(find(iz_a)) <= size(find(iz_b))
%     iz = iz_a;
% else
%     iz = iz_b;
% end

iz_ind = find(~iz);
tri_ob_before = triangulation(tri_before.ConnectivityList,clp');
test_i=ismember(tri_ob_before.ConnectivityList,iz_ind);
test_i2 = test_i(:,1)|test_i(:,2)|test_i(:,3);
tri_ob_before_Con = tri_before.ConnectivityList;
tri_ob_before_Con(test_i2,:)=[];
tri_ob_before = triangulation(tri_ob_before_Con,clp');
tri_ob_after = triangulation(tri_ob_before_Con,clp_a');

end