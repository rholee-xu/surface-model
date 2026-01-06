function[z_mid, r_mid] = compute_mid_points(rv,N)

z_mid = zeros(1, N);
for i=N:-1:1
    z_mid(1,i) = (rv(1,i+1) +rv(1,i))/2;
end

r_mid = zeros(1, N);
for i=N:-1:1
    r_mid(1,i) = (rv(2,i+1) +rv(2,i))/2;
end

end