%derivative of FunV
function [dfdx]= Calc_dfdx(X,F)
     c=5;
     n=19;
     m=2*c+1;
     dfdx=zeros((m)*n,(m)*n);     
     % stage 1
     B=gen_mat_B(X,F,1);
     C=gen_mat_C(X,1);
     dfdx(1:m,1:2*(m))=horzcat(B,C);
     
     % stage n
     A=gen_mat_A(X,n);
     B=gen_mat_B(X,F,n);
     dfdx(1+m*(n-1):m*n,(m)*n+1-2*(m):(m)*n)=horzcat(A,B);
     
     for j=2:n-1
         A=gen_mat_A(X,j);
         B=gen_mat_B(X,F,j);
         C=gen_mat_C(X,j);
         dfdx(1+m*(j-1):m*j,1+m*(j-2):m*(j+1))=horzcat(A,B,C);
     end
end