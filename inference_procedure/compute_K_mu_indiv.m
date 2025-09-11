function [K, mu] = compute_K_mu_indiv(final_tensions,final_stretch,A_after)

K = zeros(1,size(final_stretch,2));
mu = zeros(1,size(final_stretch,2));

for i = 1:size(final_stretch,2)
    tc = abs(final_tensions(8,i));
    tm = abs(final_tensions(7,i));
    K(i)=(tm+tc)/(2*((final_stretch(2,i).*final_stretch(3,i))-1));
    % mu(i) = (tm-tc)/((1/(mean(final_stretch(3,i))^2))-(1/(mean(final_stretch(2,i))^2)));
    mu(i) = (tm-tc)/((1/(final_stretch(3,i))^2)-(1/(final_stretch(2,i))^2));
end     

% w_K = mean(K.*A_after)./mean(A_after);
% w_mu = mean(mu.*A_after)./mean(A_after);

end