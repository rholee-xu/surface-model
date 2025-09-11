function[processed_data]=pre_process(originalResult,trialNum,N_I,s,new)
%
% global trialNum N_I
%
    N_I = size(originalResult,2);

if new == 1
originalResult = originalResult(:,s);
else
end
for h=1:trialNum
    if new == 1
    TF=isoutlier(originalResult(h,:),'movmedian',10);
    else
    TF=isoutlier(originalResult(h,:),'median');
    end
    for i=1:N_I
        
        if TF(i)==1
            originalResult(h,i)=NaN;   
        % elseif originalResult(h,i)>=1000
        %     originalResult(h,i)=NaN;   
        end
    end
end
%}

processed_data=originalResult;

end
        