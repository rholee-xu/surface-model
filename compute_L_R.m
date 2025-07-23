function[Lm, Rm] = compute_L_R(rv, N)

Lm = zeros(N,1);
Rm = zeros(N,1);
for i=1:N
   Lm(i,1) = sqrt((rv(1,i+1)-rv(1,i)).^2+(rv(2,i+1)-rv(2,i)).^2);
   Rm(i,1) = (rv(2,i+1)+rv(2,i))/2;
end

