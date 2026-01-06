function [new_long_axis] = get_long_axis(angle,y_tilt_angle)


% tp = pointCloud(ridge');
% t = pcnormals(tp);
% 
% [newt,~]=max(ridge(1,ridge(3,:) == z_stack));
% 
% i = find(ridge(1,:) == newt);
% 
% new_long_axis = t(i,:);


a = cosd(angle);
b = -sind(y_tilt_angle);
if angle>=0
    c = -sqrt(1-(cosd(angle)^2));
else
    c = sqrt(1-(cosd(angle)^2));
end

new_long_axis = [a,b,c];

new_long_axis = new_long_axis/norm(new_long_axis);
