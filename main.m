% given variables
n=17;
c=5;
feed_stage = 7;
feed = [10 30 5 20 20];
total_feed = sum(feed);
P = zeros(n);

% Pressure variable Pj
P(:) = 101325; %in pa

% Feed variable Fij
F = zeros(c,n);
F(:,feed_stage) = feed;

% Composition variable Zij
Z=F./total_feed;

% calc k from function
k=[0.2 2.54 1.2 0.3 1];
[L,V,x,y]=rachford_rice(total_feed, k, Z(:,feed_stage));

% molar flows
l=x.*L;
v=y.*V;
