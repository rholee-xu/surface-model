function coordArr = rotateAbtLine(arr, slope,axis)
% sintheta = slope/sqrt(slope^2+1);
% costheta = 1/sqrt(slope^2+1);
if axis == 1
    costheta = cos(slope);
    sintheta = sin(slope);

    transform = [1 0 0; 0 costheta -sintheta ;
        0 sintheta costheta];

    coordArr = transform*arr;
elseif axis == 2
    costheta = cos(slope);
    sintheta = sin(slope);

    transform = [costheta 0 sintheta ;
        0 1 0;-sintheta 0 costheta];

    coordArr = transform*arr;
else
    costheta = cos(slope);
    sintheta = sin(slope);

    transform = [costheta -sintheta 0;
        sintheta costheta 0; 0 0 1];

    coordArr = transform*arr;
end
end