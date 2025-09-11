function [int_point,tri_proj,scon_proj] = project_points_sim(subset_beads,center,plane_pt,tip,R)


int_point = zeros(size(subset_beads));

for n = 1:size(subset_beads,2)
    bead = subset_beads(:,n);
    if bead(1) >= tip - R
        p0 = [1;1;plane_pt];
        np = p0 - [1;1;plane_pt-1];
%         p0 = [tip+20;center(2);plane_pt+20];
%         np = p0 - [tip+19;center(2);plane_pt+19];
        center_n = [tip - R;center(2);center(3)];
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
    % hold on;plot3(project_line(1,:),project_line(2,:),project_line(3,:),'LineWidth',1.2,'Color',[0,0,0,0.2])

end

scon_proj= delaunay(int_point(1,:),int_point(2,:));
tri_proj = triangulation(scon_proj,subset_beads');
% triplot(scon_proj,int_point(1,:),int_point(2,:))
end