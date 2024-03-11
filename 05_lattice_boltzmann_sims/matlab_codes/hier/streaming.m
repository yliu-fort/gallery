function out = streaming(f, dirx, diry,nx,ny,ndir,mask)

buf = reshape(f,nx,ny, ndir);
for i = 1:ndir
    buf(:,:,i) = circshift2d(buf(:,:,i),[dirx(i),diry(i)]);
end
buf = reshape(buf, nx*ny, ndir);

% Output
out = buf;

% If masked, only propagation happened inside masked region is allowed
if(nargin > 6)
    mask = reshape(mask,nx,ny);
    for i = 1:ndir
        mask_offset = mask&circshift2d(mask,[-dirx(i),-diry(i)]);
        mask_offset = mask_offset(:);
        out(:,i) = mask_offset.*buf(:,i) + ~mask_offset.*f(:,i);
    end
end

    function matrix = circshift2d(matrix, dir)
        if (dir(1) < 0)
            matrix = [matrix(2:end,:);matrix(1,:)];
        elseif (dir(1) > 0)
            matrix = [matrix(end,:);matrix(1:end-1,:)];
        end
        
        if (dir(2) < 0)
            matrix = [matrix(:,2:end) matrix(:,1)];
        elseif (dir(2) > 0)
            matrix = [matrix(:,end) matrix(:,1:end-1)];
        end
        
    end

end

