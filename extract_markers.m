function [markers_extracted]=extract_markers(markers129)

% from deng et al. 2022

% global N N_I
N = 128;
N_I = 16;
indices=ones(1,N_I+1);
markers_extracted=zeros(size(markers129,1),N_I+1);

markers_extracted(:,1)=markers129(:,1); % the first marker

% indices were all one at first, by the loop below:
%  the second element in indices would be 1+n/n_i
%   for instance, with 8 markers, it would be 9, the third will be 17...
% not initialize all to zeros because it depends on the previous index

for i=2:N_I+1
   indices(i)=indices(i-1)+N/N_I;
   
   markers_extracted(:,i)=markers129(:,indices(i));
   
end

end