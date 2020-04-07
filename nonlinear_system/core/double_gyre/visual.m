figure;
load result.dat
loglog(result(:,1),result(:,2),'.'),hold on
loglog(result(:,1),result(:,1).^3.39/200,'--')
loglog(result(:,1),result(:,1)*10,'--')
legend("\delta^2","t^3^.^3^9","t^1")
grid on
xlabel("t")
ylabel("\delta^2")
title("double gyre")
set(gca, "color","w")
set(gcf, "color","w")