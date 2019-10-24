function [ M ] = rough(X,f)
c=5;n=19;reflux=9.5;
M=zeros(c,n);
for j=2:18
    for i=1:5
    M(i,j)=X(c+1+i,j)+X(i,j)-1*X(c+1+i,j-1)-1*X(i,j+1)-f(i,j);
    end
end 
for i=1:5
    M(i,1)=X(c+1+i,1)*(1+1/reflux)+X(i,1)-1*X(i,2);
    M(i,19)=X(c+1+i,19)+X(i,19)-1*X(c+1+i,18);
end
end