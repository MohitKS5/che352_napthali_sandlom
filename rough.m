
%% graphs
%% graphs
figure(1);
grid on;
hold on
X(6,:)=X(6,:)+(n+1);
plot([1:19],X(6,:),'k-o','MarkerFaceColor','k');
title('Temperature profile');
xlabel('tray number')
ylabel('Temperature(kelvin)');
 
Ls=zeros(1,19);
Vs=zeros(1,19);
for i=1:19
    Ls(i)=sum(X(c+2:2*c+1,i));
end
for i=1:19
    Vs(i)=sum(X(1:c,i));
end

figure(2)
hold on
grid on
h1=plot([1:19],X(1,:)./Vs,'k-o');
h2=plot([1:19],X(2,:)./Vs,'r-o');
h3=plot([1:19],X(3,:)./Vs,'b-o');
h4=plot([1:19],X(4,:)./Vs,'g-o');
h5=plot([1:19],X(5,:)./Vs,'y-o');
legend([h1 h2 h3 h4 h5],{'methanol','acetone','methyl acetate','benzene','choloroform'});
title(' Composition profile  for vapor phase');
xlabel('tray number')
ylabel('Vapor molar flow (kg/h');
figure(3)
hold on
grid on
h1=plot([1:19],X(1+c+1,:)./Ls,'k-o');
h2=plot([1:19],X(2+c+1,:)./Ls,'r-o');
h3=plot([1:19],X(3+c+1,:)./Ls,'b-o');
h4=plot([1:19],X(4+c+1,:)./Ls,'g-o');
h5=plot([1:19],X(5+c+1,:)./Ls,'y-o');
legend([h1 h2 h3 h4 h5],{'methanol','acetone','methyl acetate','benzene','choloroform'});
title(' Composition profile  for liquid phase');
xlabel('tray number')
ylabel('Liquid molar flow (kg/h');
hold off


