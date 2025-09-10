function [K,MU]=setup_K_MU(k,mu,rv,N,slope,endk)

z_end = max(rv(1,:));

% CONSTANT
% both uniform
% {
if slope == 0
K = k*ones(N,1);
MU = mu*ones(N,1);
% }

% LINEAR
else
[zm,~]=compute_mid_points(rv,N);
slope= (k-endk)/(max(rv(1,:))-min(rv(1,:)));

K=zeros(N,1);
for i=1:N
    K(i,1) = k - slope * zm(i);
end
MU=zeros(N,1);
for i=1:N
    MU(i,1)= mu- slope * zm(i);
end

% for longer length - needed to move the location of slope change to be at
% the transition zone

if z_end >= 4
slope= -(k-endk)/2;
rv_sub = rv(1,:)>=(z_end-2);
y_int = endk-slope*z_end;
x = rv(1,rv_sub);
y = slope*x+y_int;


x_con = rv(1,rv(1,:)<(z_end-2));
y_con = repmat(k,1,size(x_con,2));

% scatter([x,x_con],[y,y_con])
K = [y_con,y(1:end-1)]';

endmu = mu/(k/endk);
slope= -(mu-endmu)/2;
y_int = endmu-slope*z_end;
x = rv(1,rv_sub);
y = slope*x+y_int;


x_con = rv(1,rv(1,:)<(z_end-2));
y_con = repmat(mu,1,size(x_con,2));

% scatter([x,x_con],[y,y_con])
MU = [y_con,y(1:end-1)]';
else
end


end
%}



