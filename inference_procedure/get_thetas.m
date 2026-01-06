function [angles_T,angles_T_mean,mm_angle,deg,b_deg,deg_all,x_all,y_all,z_all] = get_thetas(tri_before,ridge,subset_before,rad,tip,sym,long_axis,new_long_axis)
% ridge = before_xy_ridge_all;
% x = ridge(1,:); y = ridge(2,:); 
% z = ridge(3,:);
[~,tipi] = max(ridge(1,:));
tp = pointCloud(ridge');
t = pcnormals(tp,10);

deg_all = [];
x_all = [];
y_all = [];
z_all = [];
           
% quiver3(x',y',z',t(:,1),t(:,2),t(:,3))
% quiver(x',y',t(:,1),t(:,2))

ab = t(:,2)./t(:,1);
angle = atan(ab);
% deg = rad2deg(angle);

% dott =dot(t,repmat([1,0,0],size(t,1),1),2);
% deg=acosd(dot(t,repmat([1,0,0],size(t,1),1),2));
deg = zeros(size(ridge,2),1);
b_deg = zeros(size(ridge,2),1);
x_deg = zeros(size(ridge,2),1);
% for i = 1:size(ridge,2)
% 
%     if ridge(2,i) > ridge(2,tipi)
%         if t(i,1)<0 && t(i,2)<0 || t(i,1)>0 && t(i,2)<0
%             t(i,1)=-t(i,1);t(i,2)=-t(i,2);
%             t(i,3) = -t(i,3);
%         else
%         end
%     else
%         if t(i,1)<0 && t(i,2)>0 || t(i,1)>0 && t(i,2)>0
%             t(i,1)=-t(i,1);t(i,2)=-t(i,2);
%             t(i,3) = -t(i,3);
%         else
%         end
%         % deg(i) = -(acosd(dott(i,:))+180);
%     end
%     deg(i) = acosd(dot(t(i,:),long_axis'));
%      x_deg(i) = acosd(dot(t(i,:),[1;0;0]));
%     b_deg(i) = deg(i);
%     % if abs(deg(i))>90
%         % deg(i) = acosd(dot(t(i,:),-long_axis));
%         % deg(i) = NaN;
%     % else
% end



for i = 1:size(ridge,2)

    if ridge(2,i) >ridge(2,tipi) && ridge(1,i) < tip-10
        if t(i,1)<0 && t(i,2)<0 || t(i,1)>0 && t(i,2)<0 %||  t(i,1)<0 && t(i,2)>0
            t(i,1)=-t(i,1);t(i,2)=-t(i,2);
            t(i,3) = -t(i,3);
        else
        end
    elseif ridge(2,i) <=ridge(2,tipi) && ridge(1,i) < tip-10
        if t(i,1)<0 && t(i,2)>0 || t(i,1)>0 && t(i,2)>0
            t(i,1)=-t(i,1);t(i,2)=-t(i,2);
            t(i,3) = -t(i,3);
        else
        end
        % deg(i) = -(acosd(dott(i,:))+180);
    end
    deg(i) = acosd(dot(t(i,:),new_long_axis));
    x_deg(i) = acosd(dot(t(i,:),[1;0;0]));
    if isequal(long_axis,[0;0;1]) == 1
           if ridge(3,i) >=0 && deg(i) >= 100 || ridge(3,i) <=0 && deg(i) < 100
            deg(i) = acosd(dot(t(i,:),-long_axis));
           else
           end
    end
end
side1_temp = ridge(2,:)<=ridge(2,tipi);
if isequal(long_axis,[0;0;1])==0 ||isequal(long_axis,[0;1;0])==0
    if sym == 1
        deg(side1_temp) = -deg(side1_temp);
       deg(deg>100) = NaN;%180-deg(deg>100);
       % deg(deg>100) = 0;
    else
    % b_deg(i) = deg(i);
       deg(deg>100) = NaN;%180-deg(deg>100);

       % deg(deg>100) = 0;
    end
else
end

x_deg(abs(x_deg)>100) = NaN;
b_deg = deg;
if isequal(long_axis,[0;0;1])==1 ||isequal(long_axis,[0;1;0])==1 
    side = x_deg < 60 ;
    deg(side) = NaN;
else
% 
end


IC = incenter(tri_before);
norms = tri_before.faceNormal;
angle_tri = dot(norms,repmat(new_long_axis',size(norms,1),1),2);
%%
% hold on;
angles_T = cell(1,size(tri_before.ConnectivityList,1));
mm_angle = zeros(2,size(tri_before.ConnectivityList,1));
for n = 1:size(tri_before.ConnectivityList,1)
    con_n = tri_before.ConnectivityList(n,:);
    points = tri_before.Points(con_n,:);

    e3 = new_long_axis;
    N_before = norms(n,:);
    v_b = cross(N_before,e3); c = dot(N_before,e3);
    v_x_b = [0 -v_b(3) v_b(2); v_b(3) 0 -v_b(1); -v_b(2) v_b(1) 0];
    Q_b = eye(3) + v_x_b + ((v_x_b^2)*(1/(1+c)));
    test_ridge = Q_b*ridge;
    rotated_before = Q_b*subset_before;
    
    side1 = (ridge(2,:)>=min(points(:,2))-5) & (ridge(2,:)<=max(points(:,2))+5) & (ridge(1,:)>=min(points(:,1))-5) & (ridge(1,:)<=max(points(:,1))+5)...
        & (ridge(3,:)>=min(points(:,3))-5) & (ridge(3,:)<=max(points(:,3))+5) ;

    x1_rot = test_ridge(1,side1); y1_rot = test_ridge(2,side1); 
    z1_rot = test_ridge(3,side1);
 
    x1 = ridge(1,side1); y1 = ridge(2,side1); 
    z1 = ridge(3,side1);
    % x2 = ridge(1,side2); y2 = ridge(2,side2); 
    % z2 = ridge(3,side2);


    deg1 = deg(side1);
    % deg2 = deg(side2);

        if abs(norms(n,3)) > 0.85
            
        else
        % if IC(n,2)<=ridge(2,tipi)
            in_tri = inpolygon(y1_rot,z1_rot,rotated_before(2,con_n),rotated_before(3,con_n)); 
            % in_tri_tip = inpolygon(yt,zt,subset_before(2,con_n),subset_before(3,con_n)); 
            if sym == 1
                deg_T = [-abs(deg1(in_tri));-abs(degt(in_tri_tip))];
            else
               
                    deg_T = [abs(deg1(in_tri))];
                    x_all = [x_all,x1(in_tri)];
                    y_all = [y_all,y1(in_tri)];
                    z_all = [z_all,z1(in_tri)];
            end

            deg_all = [deg_all;deg_T];
            
            % hold on;scatter3([x1(in_tri),xt(in_tri_tip)],[y1(in_tri),yt(in_tri_tip)],[z1(in_tri),zt(in_tri_tip)])
            % hold on;scatter3([x1(in_tri)],[y1(in_tri)],[z1(in_tri)])
            % text(IC(n,1),IC(n,2),IC(n,3),num2str(n))
        % else
%             in_tri = inpolygon(y2_rot,z2_rot,rotated_before(2,con_n),rotated_before(3,con_n));
%             % in_tri_tip = inpolygon(yt,zt,subset_before(2,con_n),subset_before(3,con_n));
% 
%             if sym == 1
%                 deg_T = [abs(deg2(in_tri));-abs(degt(in_tri_tip))];
%             else
%                     deg_T = [abs(deg2(in_tri))];
%                     x_all = [x_all,x2(in_tri)];
%                     y_all = [y_all,y2(in_tri)];
%                     z_all = [z_all,z2(in_tri)];
%             end
%             hold on;scatter3([x2(in_tri)],[y2(in_tri)],[z2(in_tri)]);
%             text(IC(n,1),IC(n,2),IC(n,3),num2str(n))
% %             deg_T = deg2(in_tri);
%             deg_all = [deg_all;deg_T];
% 
%         end

        deg_T_sub = deg_T(deg_T >= mean(deg_T,'omitnan')-2*std(deg_T,'omitnan') & deg_T <= mean(deg_T,'omitnan')+2*std(deg_T,'omitnan'));
        if sym == 0
            if max(deg_T_sub)-min(deg_T_sub) >=70
            % max(deg_T_sub)-min(deg_T_sub)
                deg_T_sub = deg_T(deg_T >= mean(deg_T,'omitnan')-1*std(deg_T,'omitnan') & deg_T <= mean(deg_T,'omitnan')+1*std(deg_T,'omitnan'));
            % max(deg_T_sub)-min(deg_T_sub)
            else
            end
        else
            if max(deg_T_sub)-min(deg_T_sub) >=90
            deg_T_sub = deg_T(deg_T >= mean(deg_T,'omitnan')-0.8*std(deg_T,'omitnan') & deg_T <= mean(deg_T,'omitnan')+0.8*std(deg_T,'omitnan'));
           
            else
            end
        end
        % deg_T_sub = deg_T(deg_T >= mean(deg_T,'omitnan')-std(deg_T,'omitnan') & deg_T <= mean(deg_T,'omitnan')+std(deg_T,'omitnan'));
        angles_T{n} = deg_T_sub;
        if isempty(angles_T{n})==1
            mm_angle(:,n) = [NaN;NaN];
        else
        mm_angle(:,n) = [min(angles_T{n});max(angles_T{n})];
        end

        end
end
%%
angles_T_mean = zeros(1,size(tri_before.ConnectivityList,1));
for i = 1:size(tri_before.ConnectivityList,1)
    angles_T_mean(i) = mean(angles_T{i},'omitnan');
    
end

