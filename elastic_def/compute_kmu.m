function [k, mu] = compute_kmu(tm,tc,sm,sc,N)



k=zeros(1,N);mu=zeros(1,N);

for i=1:N
    k(i)=(tm(i)+tc(i))/(2*(sm(i)*sc(i)-1));
    mu(i)=(tm(i)-tc(i))/((1/(sc(i).^2))-(1/(sm(i).^2)));

end

end
