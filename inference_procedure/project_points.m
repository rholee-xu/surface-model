function [int_point,tri_proj,scon_proj] = project_points(subset_beads,ridge,R,plane)

minz = min(subset_beads(3,:));
[tip,tipi] = max(ridge(1,:));
tip2 = ridge(2,tipi);
plane_pt = max(subset_beads(3,:)) + 20; % 2D plane is 10 pixels above max z being projected
center = [0;tip2;minz-20]; % center is the y point of the tip, and the minimum z point - 10

int_point = zeros(size(subset_beads));
hold on;
for n = 1:size(subset_beads,2)
    bead = subset_beads(:,n);
    if bead(1) >= tip - (2*R) +20 % For the tip points past a certain point
        p0 = [1;1;plane_pt];
        np = p0 - [1;1;plane_pt-1];
%         p0 = [tip+20;center(2);plane_pt+20];
%         np = p0 - [tip+19;center(2);plane_pt+19];
        center_n = [tip - (2*R)+20;center(2);center(3)];
    else
        p0 = [1;1;plane_pt];
        np = p0 - [1;1;plane_pt-1];
        center_n = [bead(1);center(2);center(3)];
    end
    l0 = bead;
    l = bead - center_n;
   
    d = dot((p0-l0),np)/dot(l,np);
    point = l0 + l*d;
 
    int_point(:,n) = point;
 
    project_line = [center_n bead int_point(:,n)];
    plot3(project_line(1,:),project_line(2,:),project_line(3,:),'LineWidth',1.2,'Color',[0,0,0,0.2])

end
if plane == 3
    scon_proj= delaunay(int_point(1,:),int_point(2,:));
elseif plane == 2
    scon_proj= delaunay(int_point(1,:),int_point(3,:));
else
    scon_proj= delaunay(int_point(2,:),int_point(3,:));
end
tri_proj = triangulation(scon_proj,subset_beads');
end