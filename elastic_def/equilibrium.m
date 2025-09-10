function [ Tl,Tr,rv,error ] = equilibrium( rv,LUT,L0,R0,K,MU,ext_verts)
    
% from deng et al. 2022 - no changes
% Solver for turgid solution given unturgid outline and elastic moduli
% distributions

    % global  Tol Rtol TolFun TolX Inc 
    Rtol = 1e-13;
    Tol =  1e-3; 
    TolFun = 1e-10; 
    TolX = 1e-10; 
    Inc = 1;

    N = size(rv,2);
    X0 = reshape(rv',2*N,1)';
    options = odeset('RelTol',Rtol);
    options2 = optimoptions('fsolve','TolFun',TolFun,'TolX',TolX,'Display','iter');
    error = 10*Tol; 
    inc = Inc;
    
    %Find initial guess near the solution
    %
    while error>Tol 
        [tX,X] = ode45(@solver,[0 inc],X0,options,LUT,L0,R0,K,MU,ext_verts);
        error = max(abs(X(end,:)-X0));
        X0=X(end,:);
    end
    %}
    %Find solution given the effective initial guess
    X = fsolve(@solver_fast,X0,options2,LUT,L0,R0,K,MU,ext_verts);
    X0=X;
    X = fsolve(@solver_fast,X0,options2,LUT,L0,R0,K,MU,ext_verts);
    error = max(abs(X-X0))
    X0=X;
    rv = reshape(X0',N,2)';  
    rb = LUT*rv';
    D = sqrt(sum(rb.^2,2)); %Vector of edge lengths.
    rm = 0.5*(rv(2,1:end-1)+rv(2,2:end))';
  
    %Compute tensions
    Tl = K.*((D./L0.*rm./R0)-1)+0.5*MU.*((R0.^2./rm.^2)-L0.^2./D.^2);
    Tr = K.*((D./L0.*rm./R0)-1)+0.5*MU.*(L0.^2./D.^2-(R0.^2./rm.^2));
    Tl=Tl.';
    Tr=Tr.';

end
