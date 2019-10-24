% material eq for the jth tray
% takes X , Feed and stage index j and returns jth stage balance
function [ M ] =M_j(X,F,j)
    %% given
    c=5;
    n=19;
    reflux=9.5;
    
    %% Equations M = out - in
    M=zeros(c,1);
    if j==1
        for i=1:c
         M(i)=X(c+1+i,1)*(1+1/reflux)+X(i,1)-1*X(i,2);
        end
    elseif j==n
        for i=1:c
          M(i)=X(c+1+i,19)+X(i,19)-1*X(c+1+i,18);
        end
    else
       for i=1:c 
        M(i)=X(c+1+i,j)+X(i,j)-1*X(c+1+i,j-1)-1*X(i,j+1)-F(i,j);
       end
    end
end
