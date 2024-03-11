figure;
%subplot(1,2,1)
p = load("123_0.25_0.0005/ccd.mean.dat");
loglog(p(:,1),gradient(p(:,2), 0.05),'-'),hold on

p = load("123_0.5_0.0005/ccd.mean.dat");
loglog(p(:,1),gradient(p(:,2), 0.05),'-'),hold on

p = load("123_1_0.0005/ccd.mean.dat");
loglog(p(:,1),gradient(p(:,2), 0.05),'-'),hold on

p = load("123_2_0.0005/ccd.mean.dat");
loglog(p(:,1),gradient(p(:,2), 0.05),'-'),hold on

p = load("123_4_0.0005/ccd.mean.dat");
loglog(p(:,1),gradient(p(:,2), 0.05),'-'),hold on

loglog(p(:,1),p(:,1).^2/10,'--')
loglog(p(:,1),p(:,1).^0*10,'--')

%%
legend("d=0.25","d=0.5","d = 1","d = 2", "d = 4","t^2","t^0","Location","SouthEast")
grid on
xlabel("t")
ylabel("d\delta^2/dt")
title("KS 3D 2048x16x1 random seed=123")
set(gca, "color","w")
set(gcf, "color","w")

% gradp

%%
subplot(1,2,2)

p = load("123_1_0.0005/ccd.mean.dat")
loglog(p(:,1)/p(1,2)^(1/3),p(:,2)/p(1,2),'.'),hold on

p = load("123_2_0.0005/ccd.mean.dat")
loglog(p(:,1)/p(1,2)^(1/3),p(:,2)/p(1,2),'.'),hold on

p = load("123_4_0.0005/ccd.mean.dat")
loglog(p(:,1)/p(1,2)^(1/3),p(:,2)/p(1,2),'.'),hold on

legend("d = 1","d = 2", "d = 4","t^3","t^1","Location","SouthEast")
grid on
xlabel("t//\delta_0^1^.^5")
ylabel("\delta^2/\delta_0^2")
title("Normalized \delta^2")
set(gca, "color","w")
set(gcf, "color","w")