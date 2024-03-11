function draw_hier_grid(gridding,NX,NY,level)
for i = 2:NX-1
    for j = 2:NY-1
        if(gridding(sub2ind([NX NY],i,j)) == 10)
            rectangle('Position',[i-2,j-2,1,1]/level*2,'FaceColor',[0 0 1])
        else
            if(gridding(sub2ind([NX NY],i,j)) == 100) % to fine
                rectangle('Position',[i-2,j-2,1,1]/level*2,'FaceColor',[0.7 0.3 0])
            else
                if(gridding(sub2ind([NX NY],i,j)) == 1) % to coarse
                    rectangle('Position',[i-2,j-2,1,1]/level*2,'FaceColor',[0 0.7 0.3])
                %else
                    %if(gridding(sub2ind([NX NY],i,j)) == 1000) % occupied
                        %rectangle('Position',[i-2,j-2,1,1],'FaceColor',[0 1 0])
                    %end
                end
            end
        end
        hold on
    end
end
view(2), axis tight,axis equal
end