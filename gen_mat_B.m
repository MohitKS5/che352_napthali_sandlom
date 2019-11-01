function [B] = gen_mat_B(X,F,j)
%% matrix B for jth plate
    c=5;
    
    %B1
    B1=zeros(1,c);
    for i=1:c
        X1=X;
        X2=X;
        % derivative using central difference
        X1(i,j)=X(i,j)+X(i,j)/(10^4);
        X2(i,j)=X(i,j)-X(i,j)/(10^4);
        B1(i)=(H_j(X1,j)-H_j(X2,j))/(2*X(i,j)/(10^4));
    end

    %B2
    X1=X;
    X2=X;
    X1(6,j)=X(6,j)+(X(6,j)/(10^4));
    X2(6,j)=X(6,j)-(X(6,j)/(10^4));
    B2=(H_j(X1,j)-H_j(X2,j))/(2*X(6,j)/(10^4));
    
    %B3
    B3=zeros(1,c);
    for i=7:11
        X1=X;
        X2=X;
        % derivative using central difference
        X1(i,j)=X(i,j)+X(i,j)/(10^4);
        X2(i,j)=X(i,j)-X(i,j)/(10^4);
        B3(i-6)=(H_j(X1,j)-H_j(X2,j))/(2*X(i,j)/(10^4));
    end
    
    %B4
    B4=zeros(c,c);
    for i=1:c
        X1=X;
        X2=X;
        % derivative using central difference
        X1(i,j)=X(i,j)+X(i,j)/(10^4);
        X2(i,j)=X(i,j)-X(i,j)/(10^4);
        B4(:,i)=(M_j(X1,F,j)-M_j(X2,F,j))./(2*X(i,j)/(10^4));
    end
  
    %B5 & B6
    B5=zeros(c,1);
    B6=zeros(c,c);
    for i=7:11
        X1=X;
        X2=X;
        % derivative using central difference
        X1(i,j)=X(i,j)+X(i,j)/(10^4);
        X2(i,j)=X(i,j)-X(i,j)/(10^4);
        B6(:,i-6)=(M_j(X1,F,j)-M_j(X2,F,j))./(2*X(i,j)/(10^4));
    end
    
    %B7
    B7=zeros(c,c);
    for i=1:c
        X1=X;
        X2=X;
        % derivative using central difference
        X1(i,j)=X(i,j)+X(i,j)/(10^4);
        X2(i,j)=X(i,j)-X(i,j)/(10^4);
        B7(:,i)=(E_j(X1,j)-E_j(X2,j))./(2*X(i,j)/(10^4));
    end
    
    %B8
    B8=zeros(c,1);
    X1=X;
    X2=X;
    % derivative using central difference
    X1(6,j)=X(6,j)+X(6,j)/(10^4);
    X2(6,j)=X(6,j)-X(6,j)/(10^4);
    B8(:,1)=(E_j(X1,j)-E_j(X2,j))./(2*X(6,j)/(10^4));
    
    %B9
    B9=zeros(c,c);
    for i=7:11
        X1=X;
        X2=X;
        % derivative using central difference
        X1(i,j)=X(i,j)+X(i,j)/(10^4);
        X2(i,j)=X(i,j)-X(i,j)/(10^4);
        B9(:,(i-6))=(E_j(X1,j)-E_j(X2,j))./(2*X(i,j)/(10^4));
    end

    B=[
        B1 B2 B3
        B4 B5 B6
        B7 B8 B9
      ];
end