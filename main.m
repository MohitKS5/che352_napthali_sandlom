%% given variables
n=19; % 17 plates + condenser + reboiler
c=5;
feed_stage = 7;
feed = [10 30 5 20 20];
total_feed = sum(feed);
P = zeros(1,n);

% Pressure variable Pj
P(:) = 1.01325; %in bar

% Feed variable Fij
F = zeros(c,n);
F(:,feed_stage) = feed;

% Composition variable Zij
Z=F./total_feed;

%% initial X (newton raphson variable)
X = zeros(2*c+1,n);

% lowest boiling point among all components
minBP=329.4;
% highest boiling point  among all components
maxBP=353.2;

% Guess temperature values by linear interpolation
for i=1:n
    X(c+1,i)=(minBP*(n-i)+maxBP*(i-1))/(n-1);
end

% calc k using antoine at feed stage (T,P)
[~,k]=antoine(X(c+1,feed_stage));
[L,V,x,y]=rachford_rice(total_feed, k, Z(:,feed_stage));

% molar flows
l=x.*L;
v=y.*V;

% use values obtained as initial guesses
for i=1:c
    X(i,:)=v(i);
    X(i+c+1,:)=l(i);
end

%% begin napthali sandlom
latent_heat_F=H_V(F(:,7)./100,330);

% form a vector for newton raphson from X
X_Vector=X(:,1)';
for i=2:19
    X_Vector= horzcat(X_Vector, X(:,i)');
end
iter=1;
tau = 1;
epsilon = inf;
while(true)
    % MEH equations
    M = zeros(c,n);
    E = zeros(c,n);
    H = zeros(c,n);
    X
    for j=1:n
        M(:,j)=M_j(X,F,j);
        E(:,j)=E_j(X,j);
        H(:,j)=H_j(X,j);
    end
    
    %function vector and matrix
    FunV=zeros(1,1); % dummy var
    FunMat=zeros((2*c+1),n);
    for j=1:n
        FunMat(1,j)=H(j);
        for i=1:c
            FunMat(1+i,j)=M(i,j);
            FunMat(1+c+i,j)=E(i,j);
        end
        FunV=horzcat(FunV,FunMat(:,j)');
    end
    FunV(1)=[]; % remove 1st dummy var
    
    % differential
    dfdx = Calc_dfdx(X,F);
    
    tau_old=tau;
    tau=sum(sum(H.^2)/(latent_heat_F^2)+sum(sum(M.^2))+sum(sum(E.^2)));
    epsilon=(sum(H.^2)+sum(sum(M.^2))+sum(sum(E.^2)))*n*(2*c+1)*(10^-10);
    
    Xnew=X_Vector'-dfdx\FunV';%% updating X
    X_Vector = Xnew';
    
    % converting X vector to matrix form
    for i=1:n
        X(:,i)=Xnew(1+(2*c+1)*(i-1):(2*c+1)*i)';
    end
    
    if abs(tau_old-tau)<(10^-5)
       break
    end
    iter=iter+1
end




