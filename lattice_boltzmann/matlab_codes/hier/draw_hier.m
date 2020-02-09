%% Draw grid
% define NX, NY
% define visibility array

% draw internal mesh
figure;
for i = 1:refine_level
    %subplot(1,refine_level,i);
    %draw_hier_grid(hier{i}.gridding,hier{i}.NX,hier{i}.NY,hier{i}.level)
    
    color = 0.*double(hier{i}.gridding == 1) + 0.5.*double(hier{i}.gridding == 100) + 1.0.*double(hier{i}.gridding == 10);
    color = repmat(color',4,1);
    height = i/refine_level*double(hier{i}.gridding~=0&hier{i}.gridding~=1000) + 0.2/refine_level*double(hier{i}.gridding==100) - 0.2/refine_level*double(hier{i}.gridding==1);
    height = repmat(height',4,1);
    patch(hier{i}.pX,hier{i}.pY,height,color,'FaceAlpha',0.25), hold on,axis tight,axis equal,view(3),grid on,zlim([0.2/refine_level 1.0])
end

figure;
for i = 1:refine_level
    %draw_hier_wireframe(hier{i}.gridding,hier{i}.NX,hier{i}.NY,hier{i}.level)
    patch('XData',hier{i}.pX,'YData',hier{i}.pY,'EdgeColor','k','FaceColor','none','LineWidth',5/hier{i}.level), hold on,axis tight,axis equal
end
