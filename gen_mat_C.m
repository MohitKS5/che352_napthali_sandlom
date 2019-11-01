function [C]=gen_mat_C(X,j)
    %% matrix C for jth plate
    X1=X;
    X2=X;
    c=5;

    C1=zeros(1,c);
    for i=1:c
        X1=X;
        X2=X;
        % derivative using central difference
        X1(i,j+1)=X(i,j+1)+X(i,j+1)/(10^4);
        X2(i,j+1)=X(i,j+1)-X(i,j+1)/(10^4);
        C1(i)=(H_j(X1,j)-H_j(X2,j))/(2*X(i,j+1)/(10^4));
    end

    % derivative using central difference
    X1(6,j+1)=X(6,j+1)+(X(6,j+1)/(10^4));
    X2(6,j+1)=X(6,j+1)-(X(6,j+1)/(10^4));
    C2=(H_j(X1,j)-H_j(X2,j))/(2*X(6,j+1)/(10^4));

    C3=zeros(1,c);
    C4=-eye(c);
    C5=zeros(c,1);
    C6=zeros(c,c);

    C7=zeros(c,c);
    for i=1:c
        X1=X;
        X2=X;
        % derivative using central difference
        X1(i,j+1)=X(i,j+1)+X(i,j+1)/(10^4);
        X2(i,j+1)=X(i,j+1)-X(i,j+1)/(10^4);
        C7(:,i)=(E_j(X1,j)-E_j(X2,j))./(2*X(i,j+1)/(10^4));
    end

    C8=zeros(c,1);
    C9=zeros(c,c);
    C=[
        C1 C2 C3
        C4 C5 C6
        C7 C8 C9
      ];
end








