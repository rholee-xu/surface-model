function [ rv, adj,R0,L0 ] = generate_sphere(R,N)
% global L
L = 1;
%GENERATE a circular strip in the first quadrant
theta = pi/2/N;
%Generate vertex positions and adjacency matrix.
rv = zeros(2,N+1);
adj = zeros(N+1);

for ii = 1:N+1
    rv(1,ii) = L.*R.*cos(pi/2-theta*(ii-1));
    rv(2,ii) = L.*(3).*sin(pi/2-theta*(ii-1));
    if (ii > 1 && ii< N+1)
        adj(ii,(ii-1)) = 1;
        adj((ii-1),ii) = 1;
        adj(ii,(ii+1)) = 1;
        adj((ii+1),ii) = 1; 
    end
end
%due to Matlab cos function error at pi/2: set boundary coordinate to 0's
%need to be 0 to tag the boundary
rv(1,1) = 0;
%Generate reference configurations.
for i=N:-1:1
   
   L0(i,1) = sqrt((rv(1,i+1)-rv(1,i)).^2+(rv(2,i+1)-rv(2,i)).^2);
   R0(i,1) = (rv(2,i+1)+rv(2,i))/2;

end
end

