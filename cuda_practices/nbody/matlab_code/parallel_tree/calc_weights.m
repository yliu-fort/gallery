function weights = calc_weights(sid,x,y,z,w,isoct,seg)
    n = numel(sid);
    weights = zeros(n-1,4);
    
    for i = 1:n-1
        if(isoct(i))
            for j = (seg(i,1)+1):(seg(i,2)+1)
                idx = sid(j);
                weights(i,1) = weights(i,1) + x(idx)*w(idx);
                weights(i,2) = weights(i,2) + y(idx)*w(idx);
                weights(i,3) = weights(i,3) + z(idx)*w(idx);
                weights(i,4) = weights(i,4) +        w(idx);
            end
            if(weights(i,4) ~= 0)
                weights(i,1) = weights(i,1)/weights(i,4);
                weights(i,2) = weights(i,2)/weights(i,4);
                weights(i,3) = weights(i,3)/weights(i,4);
            end
        end
    end
    
end