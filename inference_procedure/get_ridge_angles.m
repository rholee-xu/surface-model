function deg = get_ridge_angles(ridge,long_axis,sym)

% ridge = before_xy_ridge_all;
% x = ridge(1,:); y = ridge(2,:); 
% z = ridge(3,:);
[~,tipi] = max(ridge(1,:));
tp = pointCloud(ridge');
t = pcnormals(tp,10);

% quiver3(x',y',z',t(:,1),t(:,2),t(:,3))
% quiver(x',y',t(:,1),t(:,2))

% ab = t(:,2)./t(:,1);
% angle = atan(ab);
% deg = rad2deg(angle);
% 
% ab100 = t(:,3)./t(:,1);
% angle100 = atan(ab100);
% deg100 = rad2deg(angle100);
deg = zeros(size(ridge,2),1);
x_deg = zeros(size(ridge,2),1);
% top_temp = ridge(3,:)>0;

for i = 1:size(ridge,2)
    
    if ridge(2,i) >ridge(2,tipi)
        if t(i,1)<0 && t(i,2)<0 || t(i,1)>0 && t(i,2)<0 %||  t(i,1)<0 && t(i,2)>0
            t(i,1)=-t(i,1);t(i,2)=-t(i,2);
            t(i,3) = -t(i,3);
        else
        end
    else
        if t(i,1)<0 && t(i,2)>0 || t(i,1)>0 && t(i,2)>0
            t(i,1)=-t(i,1);t(i,2)=-t(i,2);
            t(i,3) = -t(i,3);
        else
        end
        % deg(i) = -(acosd(dott(i,:))+180);
    end
    deg(i) = acosd(dot(t(i,:),long_axis));
    x_deg(i) = acosd(dot(t(i,:),[1;0;0]));
    if isequal(long_axis,[0;0;1]) == 1
           if ridge(3,i) >=0 && deg(i) >= 100 || ridge(3,i) <=0 && deg(i) < 100
            deg(i) = acosd(dot(t(i,:),-long_axis));
           else
           end
    end
end
side1_temp = ridge(2,:)<=0;
if isequal(long_axis,[0;0;1])==0 ||isequal(long_axis,[0;1;0])==0
    if sym == 1
        deg(side1_temp) = -deg(side1_temp);
       deg(deg>100) = 180-deg(deg>100);
       % deg(deg>100) = 0;
    else
    % b_deg(i) = deg(i);
       deg(deg>100) = 180-deg(deg>100);
       % deg(deg>100) = 0;
    end
else
end

% deg(tipi) = 0;
tip_i_all = (ridge(1,:)==max(ridge(1,:)));
deg(tip_i_all) = 0;



% after_sid = ridge(1,:)>2.7;
% deg(after_sid & deg'>84) = NaN;
 
end
