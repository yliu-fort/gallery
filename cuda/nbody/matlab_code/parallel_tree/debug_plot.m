%% Plot
tree.origin.x = tree.bbox(:,1) - tree.bbox(:,4)/2;
tree.origin.y = tree.bbox(:,2) - tree.bbox(:,4)/2;
tree.origin.z = tree.bbox(:,3) - tree.bbox(:,4)/2;
tree.origin.w = tree.bbox(:,4);
figure(1), axis equal,axis([0 1 0 1 0 1]),hold on
%scatter3(p.x,p.y,p.z,[],'.')
scatter3(tree.mass(:,1),tree.mass(:,2),tree.mass(:,3),[],'.')
for i = 1:n-1
    if(tree.origin.w(i) == 0),continue;end
    plotcube([tree.origin.w(i) tree.origin.w(i) tree.origin.w(i)],[tree.origin.x(i) tree.origin.y(i) tree.origin.z(i)],0,[1 1 1])
    drawnow
end

axis equal, hold off