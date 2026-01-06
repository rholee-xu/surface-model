function [ zv, adj,R0,L0 ] = generate_arc_hyphoid(z_end,N,a)

% INPUTS
% z_end: length of cell
% N: number of marker points
% a: taper parameter

z = @(r) (pi*r)/a.*cot(pi*r)-(1/a) +z_end; % hyphoid equation

curv = @(r) sqrt(1+((pi)/a.*cot(pi*r)-(pi*r*pi)/a.*(csc(pi*r).^2)).^2); %sqrt(1+z'(r)^2) - arclength


deltas = 1e-6;
ds = integral(curv,0,base_height(a,-z_end));%compute total arclength
ds = ds/N;%the arclength covered by each linear segment
%Generate vertex positions and adjacency matrix.
zv = zeros(2,N+1);
adj = zeros(N+1);

zv(2,1)=0;
zv(1,1) = z_end;


zz = 0;
for i=1:N

    di = 0;
    while di<ds         
            di = di + deltas;
            if i == 1
                zz = zz+deltas./1;
            else
            zz = zz + deltas./curv(zz);
            end
    end
    zv(2,i+1) = zz;
    zv(1,i+1)= z(zv(2,i+1));

end
% zv(2,end) = zv(2,end-2)+(-zv(1,end-2)*(zv(2,end-1)-zv(2,end-2))/(zv(1,end-1)-zv(1,end-2)));
% zv(2,end) = 1;
zv(2,end) = base_height(a,-z_end);
zv(1,end) = 0;%zv(1,i+1);

for ii = 1:N+1
    if (ii > 1 && ii< N+1)
        adj(ii,(ii-1)) = 1;
        adj((ii-1),ii) = 1;
        adj(ii,(ii+1)) = 1;
        adj((ii+1),ii) = 1; 
    end
end


%Generate reference configurations.
for i=N:-1:1
   
   L0(i,1) = sqrt((zv(1,i+1)-zv(1,i)).^2+(zv(2,i+1)-zv(2,i)).^2);
   R0(i,1) = (zv(2,i+1)+zv(2,i))/2;
   
end


end

function mid = base_height(a,endx)
    % returns the y s.t. pi*y/a * cot*pi*y) = endx. this is the base of the
    % cell
    f = @(y) pi*y/a .* cot(pi*y) - (1/a); % the profile curve in the form x = f(y)
    yl = 0; % bisection method
    yr = 0.5;
    while f(yr) > endx
        yr = (1 + yr) / 2;
    end
    mid = (yl+yr)/2;
    fmid = f(mid);
    while abs(fmid - endx) > 0.001
        if fmid < endx
            yr = mid;
        else
            yl = mid;
        end
        mid = (yl+yr)/2;
        fmid = f(mid);
    end
    % mid = mid-0.05;
end

