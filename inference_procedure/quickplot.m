function [] = quickplot(data,l)

    

    if nargin > 1
    if size(data,1) == 3
        plot3(data(1,:),data(2,:),data(3,:),'LineWidth',l)
    else
        plot(data(1,:),data(2,:),'LineWidth',l)
    end
    axis equal; grid on
    else
    if size(data,1) == 3
        plot3(data(1,:),data(2,:),data(3,:))
    else
        plot(data(1,:),data(2,:))
    end
    axis equal; grid on
    end
end