function [rvpb] = perturb_new(rv,magni_std,magni_z_std)
% global magni
% rb = LUT*rv';
% perturbation relative to segment width
% pbm = magni * max(max(abs(rv(2,:))));
% perturbation relative to segment length
%D = sqrt(sum(rb.^2,2));
%pbm = magni*max(max(abs(D)));
% pb = pbm * (rand(size(rv)) - 0.5);
% pb(1,1)=0;%no perturbation at the boundary along x
% pb(2,end)=0;%no perturbation at the tip along x
% rvpb = rv +pb; % the perturbed profile
% pbm = zeros(1,size(rv,2));
% for i = 1:size(rv,2)
%     pbm(i) = magni* (max(abs(rv(2,i)))+0.01);
% end

rv = [rv;zeros(1,size(rv,2))];
% pb = pbm.* (rand(size(rv)) - 0.5);
% pb = magni.*(normrnd(0,0.045,size(rv)));
pb = (normrnd(0,magni_std,2,size(rv,2)));
pbz = (normrnd(0,magni_z_std,1,size(rv,2)));

% dist = [normrnd(0.5/40,0.01,1,500),normrnd(-0.5/40,0.01,1,500)];
% idx = randi(1000,2,size(rv,2));
% dist2 = [normrnd(0.5/40,0.02,1,500),normrnd(-0.5/40,0.02,1,500)];
% idx2 = randi(1000,1,size(rv,2));
% pb = dist(idx);

 % pb1 = dist(idx);
 % pb2 = dist2(idx2);
rvpb = rv(1:2,:) +pb;
rvpb(3,:) = rv(3,:)+pbz;