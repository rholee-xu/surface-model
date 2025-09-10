function [K128,MU128,synCurvM128,synCurvC128,synTensM128,synTensC128,synStrainM128,synStrainC128,...
    rv_final128,rv_init128,rv_init,rv_final,synCurvM_i128,synCurvC_i128,error] = main_material_function(a,endk,scaling,slope,k,mu)

% adapted from yaqi/min code used in deng et al. 2022
% FOR COMPUTING MATERIAL PROPERTIES K AND MU

% INPUTS
% a: controls hyphoid taper
% endk: value of tip K
% scaling: used to adjust width of cell
% slope: 0 for constant, 1 for linear decreasing
% k: side k
% mu: side mu

% UNCOMMENT FOR LINEAR/CONSTANT or NONLINEAR (see below)

% global Tol Rtol TolFun TolX Inc L N_I N trialNum 

L = 1; 
P = 2;

% number of marker points
N_I = 16; % 8,16,32,64,128
N = 128;


% computational parameters
% Rtol = 1e-13;
% Tol =  1e-3; 
% TolFun = 1e-16; 
% TolX = 1e-16; 
% Inc = 1;


%% generate the high resolution profile of cell shape
z_end = 4*L; % controls length of the cell
[rv128,adj,R0,L0] = generate_arc_hyphoid(z_end,N,a); % generates initial hyphoid outline based on parameter a

LUT = zeros(N,size(adj,1)); % oriented adjacency map
ii = 1;
for n = 1:N
    nv = find(adj(n,:));
    nv = nv(nv>n); 
    for nn = nv
        LUT(ii,n) = 1; LUT(ii,nn) = -1;
        ii = ii + 1;
    end
end
nbonds = sum(adj(:))/2;
% initial curvatures and tensions
rv128 = flip(rv128,2);
rv128 = rv128*scaling;
[synCurvM_i128,synCurvC_i128,synTensM_i128,synTensC_i128] = compute_curvatures(rv128,LUT,P);
rv_init128 = rv128;   % initial shape


% for linear varying or constant k and mu
% [K128,MU128] = setup_K_MU(k,mu,rv_init128/scaling,N,slope,endk);

% for nonlinear K and mu
[K128]=tanh_kmu(z_end,k,endk,rv_init128/scaling,N);
endmu = mu/(k/endk);
[MU128]=tanh_kmu(z_end,mu,endmu,rv_init128/scaling,N);



ext_verts = find(sum(adj,1)==1); %boundary vertices

% add turgor and grow
[rv128, adj,~,~ ] = generate_sphere(2*z_end,N);
[Tl,Tr,rv_final128,error] = equilibrium(rv128,LUT,flip(L0),flip(R0),K128,MU128,ext_verts);
rv_final128 = rv_final128*scaling;

% measure turgid parameters

[synCurvM128,synCurvC128,synTensM128,synTensC128] = compute_curvatures_turgid(rv_final128,LUT,P);
[synStrainM128,synStrainC128]=compute_strains(rv_init128,rv_final128,N);
[K_inf128, Mu_inf128]= compute_kmu(synTensM128,synTensC128,synStrainM128,synStrainC128,N);
% quickscatter(rv_init128);hold on; quickscatter(rv_final128);
% max(synStrainC128)
%% extract the markers for comparison with inference
%   as synthetic data "no noise"

% gets the 16 marker points from the n=128 outline
[rv_init]=extract_markers(rv_init128);
[rv_final]=extract_markers(rv_final128);



