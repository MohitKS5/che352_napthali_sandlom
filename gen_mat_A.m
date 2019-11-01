function [A]= gen_mat_A(X,j)
    %% matrix A  for jth plate calculation
    X1=X;
    X2=X;
    c=5;
    
    A1=zeros(1,c);
    % diff using central difference
    X1(c+1,j-1)=X(c+1,j-1)+(X(c+1,j-1)/(10^4));
    X2(c+1,j-1)=X(c+1,j-1)-(X(c+1,j-1)/(10^4));
    A2=(H_j(X1,j)-H_j(X2,j))/(2*X(6,j-1)/(10^4));
    
    A3=zeros(1,c);
    for i=7:11
        X1=X;
        X2=X;
        X1(i,j-1)=X(i,j-1)+X(i,j-1)/(10^4);
        X2(i,j-1)=X(i,j-1)-X(i,j-1)/(10^4);
        A3(i-6)=(H_j(X1,j)-H_j(X2,j))/(2*X(i,j-1)/(10^4));
    end
    
    %A4 to A9
    A4=zeros(c,c);
    A5=zeros(c,1);
    A6=-eye(c); % Identity cxc
    A7=zeros(c,c);
    A8=zeros(c,1);
    A9=zeros(c,c);
  
    A=[
        A1 A2 A3
        A4 A5 A6
        A7 A8 A9
      ];
end
