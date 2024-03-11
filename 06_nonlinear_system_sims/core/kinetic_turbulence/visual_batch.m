batch_load;

%%
flist = batch_select(batch_select(...
paramlist,'INIT_DIST',paramlist(6).INIT_DIST),'TAU_INV', 10);

%%
%figure;
%subplot(1,3,1)

for i = 1:numel(flist)
    p = load([flist(i).OUTPUT_FOLDER+"/ccd.mean.dat"])
    loglog(p(:,1),p(:,2),'-'),hold on
end
loglog(p(:,1),p(:,1).^3/100,'--')
loglog(p(:,1),p(:,1)*10,'--')

%%
legend("has gravity","no gravity","Location","SouthEast")
grid on
xlabel("t")
ylabel("\delta^2")
title("1024x16x1 tau_inv = 10, initial distance = 1x")
set(gca, "color","w")
set(gcf, "color","w")

% gradp

%%
subplot(1,3,2)
flist = batch_select(batch_select(...
paramlist,'TOL',1e-3),'INIT_DIST',paramlist(4).INIT_DIST);

p = load([flist(1).OUTPUT_FOLDER+"/ccd.mean.dat"])
loglog(p(:,1),p(:,2),'-'),hold on

p = load([flist(2).OUTPUT_FOLDER+"/ccd.mean.dat"])
loglog(p(:,1),p(:,2),'.'),hold on

p = load([flist(3).OUTPUT_FOLDER+"/ccd.mean.dat"])
loglog(p(:,1),p(:,2),'+'),hold on
loglog(p(:,1),p(:,1).^3/100,'--')
loglog(p(:,1),p(:,1)*10,'--')

%%
legend("dt=0.005","dt=0.002", "dt=0.0005","t^3","t^1","Location","SouthEast")
grid on
xlabel("t")
ylabel("\delta^2")
title("2048x16x1 tol=1e-3, Init distance = 1x")
set(gca, "color","w")
set(gcf, "color","w")

% gradp

%%
flist = batch_select(batch_select(...
paramlist,'TOL',1e-3),'INIT_DIST',paramlist(9).INIT_DIST);

subplot(1,3,3)

p = load([flist(1).OUTPUT_FOLDER+"/ccd.mean.dat"])
loglog(p(:,1),p(:,2),'-'),hold on

p = load([flist(2).OUTPUT_FOLDER+"/ccd.mean.dat"])
loglog(p(:,1),p(:,2),'.'),hold on

p = load([flist(3).OUTPUT_FOLDER+"/ccd.mean.dat"])
loglog(p(:,1),p(:,2),'+'),hold on
loglog(p(:,1),p(:,1).^3/100,'--')
loglog(p(:,1),p(:,1)*10,'--')

%%
legend("dt=0.005","dt=0.002", "dt=0.0005","t^3","t^1","Location","SouthEast")
grid on
xlabel("t")
ylabel("\delta^2")
title("2048x16x1 tol=1e-2, Init distance = 2x")
set(gca, "color","w")
set(gcf, "color","w")

% gradp