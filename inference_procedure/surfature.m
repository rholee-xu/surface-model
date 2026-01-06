function [K,H,DP_max,DP_min,VP_max1,VP_max2,VP_max3,VP_min1,VP_min2,VP_min3,n] = surfature(X,Y,Z)
% SURFATURE -  COMPUTE GAUSSIAN AND MEAN CURVATURES OF A SURFACE
%   [K,H] = SURFATURE(X,Y,Z), WHERE X,Y,Z ARE 2D ARRAYS OF POINTS ON THE
%   SURFACE.  K AND H ARE THE GAUSSIAN AND MEAN CURVATURES, RESPECTIVELY.
%   SURFATURE RETURNS 2 ADDITIONAL ARGUEMENTS,
%   [K,H,Pmax,Pmin] = SURFATURE(...), WHERE Pmax AND Pmin ARE THE MINIMUM
%   AND MAXIMUM CURVATURES AT EACH POINT, RESPECTIVELY.


% First Derivatives
[Xu,Xv] = gradient(X);
[Yu,Yv] = gradient(Y);
[Zu,Zv] = gradient(Z);

% Second Derivatives
[Xuu,Xuv] = gradient(Xu);
[Yuu,Yuv] = gradient(Yu);
[Zuu,Zuv] = gradient(Zu);

[Xuv,Xvv] = gradient(Xv);
[Yuv,Yvv] = gradient(Yv);
[Zuv,Zvv] = gradient(Zv);

% Reshape 2D Arrays into Vectors
Xu = Xu(:);   Yu = Yu(:);   Zu = Zu(:); 
Xv = Xv(:);   Yv = Yv(:);   Zv = Zv(:); 
Xuu = Xuu(:); Yuu = Yuu(:); Zuu = Zuu(:); 
Xuv = Xuv(:); Yuv = Yuv(:); Zuv = Zuv(:); 
Xvv = Xvv(:); Yvv = Yvv(:); Zvv = Zvv(:); 

Xu          =   [Xu Yu Zu];
Xv          =   [Xv Yv Zv];
Xuu         =   [Xuu Yuu Zuu];
Xuv         =   [Xuv Yuv Zuv];
Xvv         =   [Xvv Yvv Zvv];

% First fundamental Coeffecients of the surface (E,F,G)
E           =   dot(Xu,Xu,2);
F           =   dot(Xu,Xv,2);
G           =   dot(Xv,Xv,2);

m           =   cross(Xu,Xv,2);
p           =   sqrt(dot(m,m,2));
n           =   m./[p p p];

% Second fundamental Coeffecients of the surface (L,M,N)
L           =   dot(Xuu,n,2);
M           =   dot(Xuv,n,2);
N           =   dot(Xvv,n,2);

[s,t] = size(Z);
% Gaussian Curvature
K = (L.*N - M.^2)./(E.*G - F.^2);
K = reshape(K,s,t);

% Mean Curvature
H = (E.*N + G.*L - 2.*F.*M)./(2*(E.*G - F.^2));
H = reshape(H,s,t);

% Principal Curvatures
Pmax = H + sqrt(H.^2 - K);
Pmin = H - sqrt(H.^2 - K);

VP_max1 = zeros(size(E,1),1);
VP_max2 = zeros(size(E,1),1);
VP_max3 = zeros(size(E,1),1);
VP_min1 = zeros(size(E,1),1);
VP_min2 = zeros(size(E,1),1);
VP_min3 = zeros(size(E,1),1);
DP_min = zeros(size(E,1),1);
DP_max = zeros(size(E,1),1);
for i = 1:size(E,1)
    Wm1 = [E(i,1) F(i,1);F(i,1) G(i,1)];
    Wm2 = [L(i,1) M(i,1);M(i,1) N(i,1)];
    Wm = Wm1\Wm2;
%     Wm = Wm2*inv(Wm1);
    [VP_n,DP_n] = eig(Wm);
    ru = [Xu(i);Yu(i);Zu(i)];
    rv = [Xv(i);Yv(i);Zv(i)];
    alpha1 = VP_n(1,1)*ru+VP_n(2,1)*rv;
    alpha2 = VP_n(1,2)*ru+VP_n(2,2)*rv;
    % VP_min1(i,:) = alpha1(1);
    % VP_min2(i,:) = alpha1(2);
    % VP_min3(i,:) = alpha1(3);
    % VP_max1(i,:) = alpha2(1);
    % VP_max2(i,:) = alpha2(2);
    % VP_max3(i,:) = alpha2(3);
    % DP_min(i,:) = DP_n(1,1);
    % DP_max(i,:) = DP_n(2,2);
    if abs(DP_n(2,2)) > abs(DP_n(1,1))
    DP_max(i,:) = DP_n(2,2);
    DP_min(i,:) = DP_n(1,1);
    VP_min1(i,:) = alpha1(1);
    VP_min2(i,:) = alpha1(2);
    VP_min3(i,:) = alpha1(3);
    VP_max1(i,:) = alpha2(1);
    VP_max2(i,:) = alpha2(2);
    VP_max3(i,:) = alpha2(3);
    else
    DP_max(i,:) = DP_n(1,1);
    DP_min(i,:) = DP_n(2,2);
    VP_min1(i,:) = alpha2(1);
    VP_min2(i,:) = alpha2(2);
    VP_min3(i,:) = alpha2(3);
    VP_max1(i,:) = alpha1(1);
    VP_max2(i,:) = alpha1(2);
    VP_max3(i,:) = alpha1(3);
    end
%     dot(testV(:,1),testV(:,2))
end
DP_max = reshape(DP_max,s,t);
DP_min = reshape(DP_min,s,t);
VP_max1 = reshape(VP_max1,s,t);
VP_max2 = reshape(VP_max2,s,t);
VP_max3 = reshape(VP_max3,s,t);
VP_min1 = reshape(VP_min1,s,t);
VP_min2 = reshape(VP_min2,s,t);
VP_min3 = reshape(VP_min3,s,t);
% % 
DP_max(isnan(DP_max))=0;
DP_min(isnan(DP_min))=0;
VP_max1(isnan(VP_max1))=0;
VP_max2(isnan(VP_max2))=0;
VP_max3(isnan(VP_max3))=0;
VP_min1(isnan(VP_min1))=0;
VP_min2(isnan(VP_min2))=0;
VP_min3(isnan(VP_min3))=0;
% 
DP_max = remove_border(DP_max);
DP_min = remove_border(DP_min);
VP_max1 = remove_border(VP_max1);
VP_max2 = remove_border(VP_max2);
VP_max3 = remove_border(VP_max3);
VP_min1 = remove_border(VP_min1);
VP_min2 = remove_border(VP_min2);
VP_min3 = remove_border(VP_min3);
% % % 
DP_max(isnan(DP_max))=0;
DP_min(isnan(DP_min))=0;
VP_max1(isnan(VP_max1))=0;
VP_max2(isnan(VP_max2))=0;
VP_max3(isnan(VP_max3))=0;
VP_min1(isnan(VP_min1))=0;
VP_min2(isnan(VP_min2))=0;
VP_min3(isnan(VP_min3))=0;
% 
DP_max = remove_border(DP_max);
DP_min = remove_border(DP_min);
VP_max1 = remove_border(VP_max1);
VP_max2 = remove_border(VP_max2);
VP_max3 = remove_border(VP_max3);
VP_min1 = remove_border(VP_min1);
VP_min2 = remove_border(VP_min2);
VP_min3 = remove_border(VP_min3);
% % % % 
DP_max(isnan(DP_max))=0;
DP_min(isnan(DP_min))=0;
VP_max1(isnan(VP_max1))=0;
VP_max2(isnan(VP_max2))=0;
VP_max3(isnan(VP_max3))=0;
VP_min1(isnan(VP_min1))=0;
VP_min2(isnan(VP_min2))=0;
VP_min3(isnan(VP_min3))=0;
% 
DP_max = remove_border(DP_max);
DP_min = remove_border(DP_min);
VP_max1 = remove_border(VP_max1);
VP_max2 = remove_border(VP_max2);
VP_max3 = remove_border(VP_max3);
VP_min1 = remove_border(VP_min1);
VP_min2 = remove_border(VP_min2);
VP_min3 = remove_border(VP_min3);
% save('surfature_variables_testing.mat')
function [P1_new] = remove_border(P1)
    test = bwmorph(P1,'remove');
    P1_new = ~test.*P1;
    P1_new(P1_new==0)=NaN;

% Pnan = isnan(P1);
% test = bwmorph(Pnan,'remove');
% [test2x,testy] = gradient(test);
% ind =abs(test2x)==0.5 | abs(testy)==0.5;
% ind = ~ind;
% test2 = ind.*P1;
% test2(test2==0)=NaN;
% test = isoutlier(P1);
% test = ~test;
% P1_new = test.*P1;
% P1_new(P1_new==0)=NaN;
% 
% P2_new = test.*P2;
% P2_new(P2_new==0)=NaN;
end
end


