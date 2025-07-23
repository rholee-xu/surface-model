function [strainM,strainC]=compute_strains(rv_init,rv_final,N)

strainM = zeros(1,N);
strainC = zeros(1,N);

[L_i,~]=compute_L_R(rv_init,N);
[L_f,~]=compute_L_R(rv_final,N);

[~,rm_i]=compute_mid_points(rv_init,N);
[~,rm_f]=compute_mid_points(rv_final,N);

for i=1:N
    strainM(i) = L_f(i) / L_i(i);
    strainC(i) = rm_f(i)/rm_i(i);
end

rv(1,:) = deg_to_sym(rv_final);

% RX - extrapolation to avoid side issues
for i = flip(1:30)
strainM(i) = strainM(i+2) + (rv(1,i)-rv(1,i+2))*(strainM(i+1)-strainM(i+2))/(rv(1,i+1)-rv(1,i+2));
strainC(i) = strainC(i+2) + (rv(1,i)-rv(1,i+2))*(strainC(i+1)-strainC(i+2))/(rv(1,i+1)-rv(1,i+2));
end

end
    