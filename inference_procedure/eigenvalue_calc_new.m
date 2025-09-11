function [eigenvalues,eigenvectors,ratio,area,triangle_areas_after,triangle_areas_before,errors] = ...
    eigenvalue_calc_new(tri_after,tri_before,pts_after,pts_before,del)

IC = incenter(tri_after);
eigenvalues = zeros(3,size(tri_before.ConnectivityList,1));
errors = zeros(1,size(tri_before.ConnectivityList,1));
eigenvectors = zeros(6,size(tri_before.ConnectivityList,1));
triangle_areas_after = zeros(1,size(tri_before.ConnectivityList,1));
triangle_areas_before = zeros(1,size(tri_before.ConnectivityList,1));
for n = 1:size(tri_before.ConnectivityList,1)
    [V,D,A_after,A_before,err_max] = strain_calculation_new(pts_after(:,tri_after.ConnectivityList(n,:)), ...
        pts_before(:,tri_before.ConnectivityList(n,:)),del);
    eigenvalues(:,n) = [IC(n);D(1,1);D(2,2)];
    eigenvectors(:,n) = [V(:,1);V(:,2)];
    errors(n) = err_max;
    triangle_areas_after(n) = A_after;
    triangle_areas_before(n) = A_before;
end

area = eigenvalues(2,:).*eigenvalues(3,:);
ratio = eigenvalues(3,:)./eigenvalues(2,:);

end