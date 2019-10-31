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
tau = 1;
epsilon = inf;
while(epsilon>tau)
    % MEH equations
    M = zeros(c,n);
    E = zeros(c,n);
    H = zeros(c,n);
    for j=1:n
        M(:,j)=M_j(X,F,j);
        E(:,j)=E_j(X,j);
    end
    
    % to close the loop for now
    epsilon = 0;
end

