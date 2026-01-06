function ridge_subset = get_ridge_new(raw_ridge,zmin,zmax)

raw_ridge_read = readtable(raw_ridge,'VariableNamingRule','preserve');

xyridge_array = table2array(raw_ridge_read(:,1:7));
ridge = xyridge_array(:,[5,6,2])';
ridge_subset = ridge(:,ridge(3,:) <= zmax & ridge(3,:) >= zmin);
%ridge(3,:) = repmat(z,[size(ridge,2),1]);

end