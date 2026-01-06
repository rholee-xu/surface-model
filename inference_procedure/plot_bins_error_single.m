function [] = plot_bins_error_single(bins_all,bins_mid,mv_all1,col,line)
if size(bins_mid,2) ~= size(mv_all1,2)
    disp('size of bins and result do not match')
    return
else
% if mv_all1(3,end) < 4 && isnan(mv_all1(1,end)) == 0 
%     mv_all1 = combine_bins(mv_all1,size(mv_all1,2)-1);
%     % bins_all(bins_all==90)=[];
%     bins_all(:,end-1)=[];
%     bins_mid = bins_all(2:size(bins_all,2))-(0.5*diff(bins_all));
% else
% end
% if mv_all1(3,1) < 4 && isnan(mv_all1(1,1)) == 0 
%     mv_all1 = combine_bins(mv_all1,1);
%     bins_all(:,2)=[];
%     bins_mid = bins_all(2:size(bins_all,2))-(0.5*diff(bins_all));
% else
% end
% if mv_all1(3,2) < 4 && isnan(mv_all1(1,2)) == 0 
%     mv_all1 = combine_bins(mv_all1,2);
%     bins_all(bins_all==20)=[];
%     bins_mid = [5,20,bins_mid(4:end)];
% else
% end

if nargin < 5
    bins = size(mv_all1,2)+1;
    hold on;
    for n = 1:bins-1
        plot([bins_all(n) bins_all(n+1)],[mv_all1(1,n) mv_all1(1,n)],'LineWidth',2,'LineStyle','-','Color',col,'Marker','none');%text(bins_all(n),mv_all1(1,n)+mv_all1(2,n),append('n=',string(mv_all1(3,n))),'FontSize',12)
        % if length(col) == 1
        % else
        % text(bins_all(n),mv_all1(1,n)+mv_all1(2,n),append('n=',string(mv_all1(3,n))),'FontSize',12)
        % end
    end
    hold on;l3 = errorbar(bins_mid(1:bins-1),mv_all1(1,1:bins-1),mv_all1(2,1:bins-1),'LineWidth',2,'LineStyle','-','Color',col,'Marker','none');
    

else
    bins = size(mv_all1,2)+1;
    hold on;
    for n = 1:bins-1
        plot([bins_all(n) bins_all(n+1)],[mv_all1(1,n) mv_all1(1,n)],'LineWidth',2,'LineStyle','--','Color',col,'Marker','none');%text(bins_all(n),mv_all1(1,n)+mv_all1(2,n),append('n=',string(mv_all1(3,n))),'FontSize',12)
        % if length(col) == 1
        % else
        % text(bins_all(n),mv_all1(1,n)+mv_all1(2,n),append('n=',string(mv_all1(3,n))),'FontSize',12)
        % end
    end
    hold on;l3 = errorbar(bins_mid(1:bins-1),mv_all1(1,1:bins-1),mv_all1(2,1:bins-1),'LineWidth',2,'LineStyle','--','Color',col,'Marker','none');
    
    
end

end

%%
function [new] = combine_bins(ww_angles,n)

N1 = ww_angles(3,n);
N2 = ww_angles(3,n+1);
s1 = ww_angles(2,n);
s2 = ww_angles(2,n+1);
X1 = ww_angles(1,n);
X2 = ww_angles(1,n+1);
X12 = mean([ww_angles(1,n),ww_angles(1,n+1)]);


% new_sd = sqrt(((N1*(s1^2))+(N2*(s2^2))+(N1*(d1^2))+(N2*(d2^2)))/(N1+N2));
term1 = ((N1-1)*(s1^2)+(N2-1)*(s2^2))/(N1+N2-1);
term2 = (N1*N2*(X1-X2)^2)/((N1+N2)*(N1+N2-1));
new_sd = sqrt(term1+term2);

new = ww_angles;
new(:,n) = [];
new(:,n) = [X12;new_sd;(N1+N2)];
end

end