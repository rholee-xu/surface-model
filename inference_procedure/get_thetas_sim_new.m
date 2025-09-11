function [angles_T,angles_T_mean,mm_angle,deg,b_deg] = get_thetas_sim_new(tri_before,ridge,subset_before,tip,rad,sym,long_axis)
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
% ab90 = t(:,3)./t(:,1);
% angle90 = atan(ab90);
% deg90 = rad2deg(angle90);
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
           if ridge(3,i) >=0 && deg(i) >= 90 || ridge(3,i) <=0 && deg(i) < 90
            deg(i) = acosd(dot(t(i,:),-long_axis));
           else
           end
    end
end
side1_temp = ridge(2,:)<=0;
if isequal(long_axis,[1;0;0]) == 1
    if sym == 1
        deg(side1_temp) = -deg(side1_temp);
       deg(abs(deg)>100) = NaN;
    else
    % b_deg(i) = deg(i);
       deg(abs(deg)>100) = 0;
    end
else
end

x_deg(abs(x_deg)>100) = NaN;
b_deg = deg;
if isequal(long_axis,[0;0;1]) == 1 || isequal(long_axis,[0;1;0]) == 1
    side = x_deg < 50 | ridge(1,:)' > 1.5;
    deg(side) = NaN;
else
% 
end

tip_i_all = (ridge(1,:)==max(ridge(1,:)));
deg(tip_i_all) = 0;

onet = find(tip_i_all,1);

deg = [deg(~tip_i_all);deg(onet)];
ridge = [ridge(:,~tip_i_all),ridge(:,onet)];

% after_sid = ridge(1,:)>2.7;
% deg(after_sid & deg'>84) = NaN;
%%
side1 = ridge(2,:)<=0 & ridge(1,:)<=tip-0.8*rad & ridge(3,:)<0.8;
side2 = ridge(2,:)>0 & ridge(1,:)<=tip-0.8*rad& ridge(3,:)<0.8;
tip_i = ridge(1,:)>tip-0.8*rad & ridge(3,:)<0.8;
top = ridge(3,:)>=0.8;
x1 = ridge(1,side1); y1 = ridge(2,side1); 
z1 = ridge(3,side1);
x2 = ridge(1,side2); y2 = ridge(2,side2); 
z2 = ridge(3,side2);
xt = ridge(1,tip_i); yt = ridge(2,tip_i); 
zt = ridge(3,tip_i);
xtop = ridge(1,top); ytop = ridge(2,top); 
ztop = ridge(3,top);

deg1 = deg(side1);
deg2 = deg(side2);
if isequal(long_axis,[0;0;1]) == 1 || isequal(long_axis,[0;1;0]) == 1
    degt = NaN(size(xt'));
else
    degt = deg(tip_i);
end
% degt90 = deg90(tip_i);
degtop = deg(top);
% scatter3(ridge(1,tip_i),ridge(2,tip_i),ridge(3,tip_i),10,degt,'filled')

IC = incenter(tri_before);
%%
angles_T = cell(1,size(tri_before.ConnectivityList,1));
mm_angle = zeros(2,size(tri_before.ConnectivityList,1));
for n = 1:size(tri_before.ConnectivityList,1)
    con_n = tri_before.ConnectivityList(n,:);
 
        % hold on;
        if IC(n,2)<=0
            in_tri = inpolygon(x1,z1,subset_before(1,con_n),subset_before(3,con_n)); 
            % hold on;scatter3(x1(in_tri),y1(in_tri),z1(in_tri),10,deg1(in_tri),'filled')
            in_tri_tip = inpolygon(yt,zt,subset_before(2,con_n),subset_before(3,con_n)); 
            % hold on;scatter3(xt(in_tri_tip),yt(in_tri_tip),zt(in_tri_tip),10,degt(in_tri_tip),'filled')
            % hold on; trisurf(tri_before.ConnectivityList(n,:),subset_before(1,:),subset_before(2,:),subset_before(3,:),'FaceAlpha',0,'LineWidth',2)
            % text(IC(n,1),IC(n,2),IC(n,3),num2str(n))
            if IC(n,3)>0
            in_tri_top = inpolygon(xtop,ytop,subset_before(1,con_n),subset_before(2,con_n)); 
          
                if isempty(find(in_tri_top,1))==1
                    
                    if sym == 1
                        deg_T = [-abs(deg1(in_tri));degt(in_tri_tip);degtop(in_tri_top)];
                    else
                        deg_T = [abs(deg1(in_tri));abs(degt(in_tri_tip));abs(degtop(in_tri_top))];
                    end
                else
                    deg_T = [abs(deg1(in_tri));abs(degt(in_tri_tip));abs(degtop(in_tri_top))];
                end

            else
                if sym == 1
                    deg_T = [-abs(deg1(in_tri));degt(in_tri_tip)];
                else
                    deg_T = [abs(deg1(in_tri));abs(degt(in_tri_tip))];
                end
       
            end
        else
            in_tri = inpolygon(x2,z2,subset_before(1,con_n),subset_before(3,con_n));
            in_tri_tip = inpolygon(yt,zt,subset_before(2,con_n),subset_before(3,con_n));
            if IC(n,3)>0
            in_tri_top = inpolygon(xtop,ytop,subset_before(1,con_n),subset_before(2,con_n)); 
                if isempty(find(in_tri_top,1))==1
                    if sym == 1
                        deg_T = [abs(deg2(in_tri));degt(in_tri_tip);degtop(in_tri_top)];
                    else
                        deg_T = [abs(deg2(in_tri));abs(degt(in_tri_tip));abs(degtop(in_tri_top))];
                    end
                else
                    deg_T = [abs(deg2(in_tri));abs(degt(in_tri_tip));abs(degtop(in_tri_top))];
                end
    
            else
        
                if sym == 1
                    deg_T = [abs(deg2(in_tri));degt(in_tri_tip)];
                else
                    deg_T = [abs(deg2(in_tri));abs(degt(in_tri_tip))];
                end

            end

            
        end

        deg_T_sub = deg_T(deg_T >= (mean(deg_T,'omitnan')-2*std(deg_T,'omitnan')) & deg_T <= (mean(deg_T,'omitnan')+2*std(deg_T,'omitnan')));
        angles_T{n} = deg_T_sub;


        if isempty(angles_T{n})==1
            mm_angle(:,n) = [NaN;NaN];
        else
        mm_angle(:,n) = [min(angles_T{n});max(angles_T{n})];
        end
end
%%
angles_T_mean = zeros(1,size(tri_before.ConnectivityList,1));
for i = 1:size(tri_before.ConnectivityList,1)
    angles_T_mean(i) = mean(angles_T{i},'omitnan');
    
end

