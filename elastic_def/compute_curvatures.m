% this part computes the curvatures and the tensions calculated
function [ ks, kphi, t_m, t_c,angle,D] = compute_curvatures(rv,LUT,P)

% from deng et al. 2022

rb = LUT*rv';%boundary vectors
%Compute instantaneous tension vectors
D = sqrt(sum(rb.^2,2)); %Vector of edge lengths.
rb(:,1) = rb(:,1)./D;
rb(:,2) = rb(:,2)./D;%director
LUTP = LUT>0;
LUTM = LUT<0;
rb1 = (LUTP+LUTM)*rv(2,:)';%yL+yR 
rm = 0.5*rb1;
rbn = [rb(:,2) -rb(:,1)];%outward normal
%angle1 = -atan(rb(:,2)./rb(:,1));%tangent 
%angle2 = atan(rbn(:,2)./rbn(:,1));%outward normal
angle = pi/2-atan(rb(:,2)./rb(:,1));%outward normal from r axis
ks = ones(1,size(D,1));
kphi = ones(1,size(D,1));
% i=1;
% %ks(i) = (3*angle2(i)-4*angle2(i+1)+angle2(i+2))./(D(i+1)+0.5*D(i)+0.5*D(i+2));%forward difference, not strictly second order, because D is not uniform
% ks(i) = (angle(i+1)+angle(i)-pi)./(2*D(i));
% kphi(i) = sin(angle(i))/rm(i);

% calculate tensions from curvuratures using Eq.36 & 37
% t_m(i) = P / (2 * kphi(i));
% t_c(i) = P / (2 * kphi(i)) * (2 - ks(i) / kphi(i));

for i=2:(size(D,1)-1)    
    ks(i) = (angle(i+1)-angle(i-1))./(2*D(i));
    kphi(i) = sin(angle(i))/rm(i);
    % calculate tensions from curvuratures using Eq.36 & 37
    t_m(i) = P / (2 * kphi(i));
    t_c(i) = P / (2 * kphi(i)) * (2 - ks(i) / kphi(i));
    % ks(i+1) = (angle(i+2)-angle(i))./(2*D(i+1));
    % kphi(i+1) = sin(angle(i+1))/rm(i+1);
    % % % calculate tensions from curvuratures using Eq.36 & 37
    % t_m(i+1) = P / (2 * kphi(i+1));
    % t_c(i+1) = P / (2 * kphi(i+1)) * (2 - ks(i+1) / kphi(i+1));
end

i=size(D,1);
%ks(i) = -(3*angle2(i)-4*angle2(i-1)+angle2(i-2))./(D(i-1)+0.5*D(i-2)+0.5*D(i));%backward difference, not strictly second order, because D is not uniform
ks(i) = ((2*pi-angle(i))-angle(i-1))./(2*D(i));
kphi(i) = sin(angle(i))/rm(i);



% calculate tensions
t_m(i) = P / (2 * kphi(i));
t_c(i) = P / (2 * kphi(i)) * (2 - ks(i) / kphi(i));

% RX - added extrapolation for side curvature

for i = [3,2,1]
ks(i) = ks(i+2) + (rv(1,i)-rv(1,i+2))*(ks(i+1)-ks(i+2))/(rv(1,i+1)-rv(1,i+2));
kphi(i) = kphi(i+2) + (rv(1,i)-rv(1,i+2))*(kphi(i+1)-kphi(i+2))/(rv(1,i+1)-rv(1,i+2));
t_m(i) = P / (2 * kphi(i));
t_c(i) = P / (2 * kphi(i)) * (2 - ks(i) / kphi(i));
end



end
