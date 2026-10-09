function SNNG=ShareNNG(DecKNNG)
    A=full(DecKNNG>0);
    n=size(A,1);
    SNNG=zeros(n,n);
    for i=1:n
        for j=(i+1):n
            SNNG(i,j)=sum(A(:,i).*A(:,j));
        end
    end
    SNNG=max(SNNG,SNNG');
end