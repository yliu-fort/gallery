%% DEBUG 2: TREE TRAVERSAL

visited = false(2*n-1,1);
node = 1;
to_visit = [-1];
su = 0;
while(node ~= -1)
    
    visited(node) = true;
    if(checkSum(visited,su)),su = su +1;else, disp("Error!");end
    
    %fprintf("Enter node %d\n",node);
    if(node >= n), continue;end
    
    % 1-based in matlab
    childL = tree.connectivity(node,1)+1;
    childR = tree.connectivity(node,2)+1;
    
    traverseL = true;
    traverseR = true;
    
    if(childL >= n)
        
        visited(childL) = true;
        if(checkSum(visited,su)),su = su +1;else, disp("Error!");end
        
        traverseL = false;
    else
        if(tree.level(childL) > 29)
            for i = tree.segment(childL,1):tree.segment(childL,2)
            
                visited(n + i) = true;
                if(checkSum(visited,su)),su = su +1;else, disp("Error!");end
            
            end
            traverseL = false;
        end
    end
    
    if(childR >= n)
        
        visited(childR) = true;
        if(checkSum(visited,su)),su = su +1;else, disp("Error!");end
        
        traverseR = false;
    else
        if(tree.level(childR) > 29)
            for i = tree.segment(childR,1):tree.segment(childR,2)
                
                visited(n + i) = true;
                if(checkSum(visited,su)),su = su +1;else, disp("Error!");end
            
            end
            traverseR = false;
        end
    end
    
    if(~traverseL && ~traverseR)
        node = to_visit(end);
        to_visit = to_visit(1:end-1);
    else
        if(traverseL), node = childL;
        else, node = childR; end
        
        if(traverseL && traverseR)
            to_visit = [to_visit;childR];
        end
    end
end

% DISPLAYd
fprintf("Visited internal nodes = %d/%d\n",sum(visited(1:n-1)),n-1);
fprintf("Visited leaf nodes = %d/%d\n",sum(visited(n:2*n-1)),n);

% function
function out = checkSum(list,s)
    if(sum(list) == s + 1)
        out = true;
    else
        out = false;
    end
end