function [tri_max_all,err_max_all] = tri_error(tri_before,tri_after,subset_beads_before,subset_beads_after,displ,del)

tri_max_all = zeros(1,size(tri_before.ConnectivityList,1));
err_max_all = zeros(1,size(tri_before.ConnectivityList,1));
for n = 1:size(tri_before.ConnectivityList,1)
    beads_before = subset_beads_before(:,tri_before.ConnectivityList(n,:));
    beads_after = subset_beads_after(:,tri_after.ConnectivityList(n,:));
    

    va = beads_before(:,1) - beads_before(:,2);
    vb = beads_before(:,1) - beads_before(:,3);
    va_0 = beads_after(:,1) - beads_after(:,2);
    vb_0 = beads_after(:,1) - beads_after(:,3);
    
    vc = beads_before(:,1) - beads_before(:,2);
    vd = beads_before(:,2) - beads_before(:,3);
    vc_0 = beads_after(:,1) - beads_after(:,2);
    vd_0 = beads_after(:,2) - beads_after(:,3);
    
    ve = beads_before(:,2) - beads_before(:,3);
    vf = beads_before(:,1) - beads_before(:,3);
    ve_0 = beads_after(:,2) - beads_after(:,3);
    vf_0 = beads_after(:,1) - beads_after(:,3);
    
    [e1,b1] = calc_err(va,vb,va_0,vb_0,displ,del);
    [e2,b2] = calc_err(vc,vd,vc_0,vd_0,displ,del);
    [e3,b3] = calc_err(ve,vf,ve_0,vf_0,displ,del);
    tri_max_all(n) = max([e1 e2 e3]);
    err_max_all(n) = max(abs([b1 b2 b3]));
end

    function [err_max,err_beads] = calc_err(v1,v2,va1,va2,displ,del)
        sintheta = norm(cross(v2,v1))/(norm(v2)*norm(v1));
        err_max = (displ/norm(v2)+displ/norm(v1))*(2/sintheta);
        sintheta_init = norm(cross(va2,va1))/(norm(va2)*norm(va1));
        err_beads = ((-2*(norm(v1)^2)*(norm(v2)^2)*sintheta)/((norm(va1)^2)*(norm(va2)^2)*(sintheta_init^2)))*((1/norm(v1))+(1/norm(v2)))*displ;
    end

end