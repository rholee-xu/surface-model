function [] = quickscatter(data,ptsize,color)
    if nargin < 2
    if size(data,1) == 3
        scatter3(data(1,:),data(2,:),data(3,:),10,'filled'); axis equal;
    else
        scatter(data(1,:),data(2,:),10,'filled'); axis equal;
    end

    elseif nargin < 3
    if size(data,1) == 3
        scatter3(data(1,:),data(2,:),data(3,:),ptsize,'filled'); axis equal;
    else
        scatter(data(1,:),data(2,:),ptsize,'filled'); axis equal;
    end

    else
        if size(data,1) == 3
        scatter3(data(1,:),data(2,:),data(3,:),ptsize,color,'filled'); axis equal;
        else
        scatter(data(1,:),data(2,:),ptsize,color,'filled'); axis equal;
        end
    end
    axis equal; grid on
end