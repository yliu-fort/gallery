clear; clc; close all;

%% Build hier grid
build_hier;

%% Draw grid
draw_hier;

%% Post-processing
figure;
for i = 1:refine_level
    %subplot(1,refine_level,i);
    %draw_hier_grid(hier{i}.gridding,hier{i}.NX,hier{i}.NY,hier{i}.level)

    color = 0.*double(hier{i}.gridding == 1) + 0.5.*double(hier{i}.gridding == 100) + 1.0.*double(hier{i}.gridding == 10);
    color = repmat(color',4,1);
    height = i/refine_level*double(hier{i}.gridding~=0&hier{i}.gridding~=1000) ...
        + 0.55/refine_level*double(hier{i}.gridding==100) ...
        - 0.5/refine_level*double(hier{i}.gridding==1);
    height = repmat(height',4,1);
    fig{i} = patch(hier{i}.pX,hier{i}.pY,height,color,'FaceAlpha',1,'EdgeColor','none'); hold on,axis tight,axis equal,view(2),grid on,zlim([0.2/refine_level 1.0])
end
colormap(desaturated_rainbow);

% mass conservation
mass_conservation_ref = 0;
for i = 1:refine_level
    mass_conservation_ref = mass_conservation_ref + ...
        sum(hier{i}.rho(hier{i}.gridding==10|hier{i}.gridding==100));
end

%%  main simulation loop; take NSTEPS time steps
convergence = Inf;
n = 0;
while((~isfinite(convergence) || convergence > 1e-15))
    
    % record statistics
    uhist = hier{refine_level}.ux;
    vhist = hier{refine_level}.uy;
    
    hier = time_integrator(1,refine_level,hier,wi,dirx,diry,ndir);
    
    % Advance
    n = n+1;
    %if(n == 1) break;end
    
    % visualization
    if(mod(n,hier{1}.PLOTGAP) == 0 || n == hier{1}.NSTEPS)
        % Convergence
        conv = sum((hier{refine_level}.ux - uhist).*(hier{refine_level}.ux - uhist) ...
            + (hier{refine_level}.uy - vhist).*(hier{refine_level}.uy - vhist),'omitnan');
        convergence = sqrt(conv/sum(uhist.^2 + vhist.^2,'omitnan'));
        %fprintf('res: %e>\n', convergence);
        
        % conservation
        mass_conservation = 0;
        for i = 1:refine_level
            mass_conservation = mass_conservation + ...
                sum(hier{i}.rho(hier{i}.gridding==10|hier{i}.gridding==100));
        end
        fprintf('res: %e | mass conservation: %.15f>\n', convergence, mass_conservation/mass_conservation_ref);
        
        % Visualization
        for i = 1:refine_level
            %subplot(1,refine_level,i);
            %draw_hier_grid(hier{i}.gridding,hier{i}.NX,hier{i}.NY,hier{i}.level)
            
            color =sqrt(hier{i}.ux.^2 + hier{i}.uy.^2);
            color = repmat(color',4,1);
            fig{i}.CData = color;
        end
        
        %rh1 = reshape((rho1),NX1,NY1);
        %uh1 = reshape((Cu1*ux1),NX1,NY1);
        %vh1 = reshape((Cu1*uy1),NX1,NY1);
        %h21.CData = sqrt(uh1.^2 + vh1.^2);
        
        %[dvdx,~] = gradient(vh);[dudy,~] = gradient(uh);
        %h3.CData = dvdx - dudy;
        %h3.CData = max(abs(dudy),abs(dvdx));
        %l2.YData = uh(floor(NX/2),:);
        %h2.CData = vh;
        %h2.CData = rh-rho0;
        %subplot(1,2,1)
        %axis equal;axis tight;axis([0 ratio(1) 0 ratio(2)]);
        %subplot(1,2,2)
        %axis equal;axis tight;axis([0 ratio(1) 0 ratio(2)]);
        %title(sprintf('%d',n/PLOTGAP))%,caxis([0 2])
        drawnow;
        
    end
end

function hier = time_integrator(level,num_level,hier,wi,dirx,diry,ndir)
% Base case
if(level < 1) ,return;end
if(level > num_level) ,return;end

% Find current level by counting leading zeros
%fprintf("Enter level %d\n",level);
%fprintf("Apply collision on level %d\n",level);
%0: unused 1: to_coarse 10: fluid 100: to_fine 1000:occupied
isFluidOrToFineCell = hier{level}.gridding == 10|hier{level}.gridding == 100;

f = hier{level}.f;
rho = sum(f,2);
ux  = sum(bsxfun(@times,f,dirx),2)./rho;
uy  = sum(bsxfun(@times,f,diry),2)./rho;

%% apply collision only on fluid and to_fine cell
hier{level}.rho(isFluidOrToFineCell) = rho(isFluidOrToFineCell);
hier{level}.ux(isFluidOrToFineCell)  = ux(isFluidOrToFineCell);
hier{level}.uy(isFluidOrToFineCell)  = uy(isFluidOrToFineCell);

feq=compute_equilibrium(hier{level}.rho, hier{level}.ux, hier{level}.uy, dirx, diry, wi);
f(isFluidOrToFineCell,:) = (1.0-1.0/hier{level}.tau)*f(isFluidOrToFineCell,:) + (1.0/hier{level}.tau)*feq(isFluidOrToFineCell,:);

%% F-redistribution on level
if(level < num_level)
    %fprintf("Apply F-redistribution on level %d\n",level);
    isToCoarseCellNextLevel = hier{level+1}.gridding==1;
    hier{level+1}.f(isToCoarseCellNextLevel,:) = ...
        f(hier{level}.interface_to_fine_connectivity(isToCoarseCellNextLevel),:)/4;
    
end

%% Launch new kernels recursively
hier = time_integrator(level+1,num_level,hier,wi,dirx,diry,ndir);
hier = time_integrator(level+1,num_level,hier,wi,dirx,diry,ndir);

%% F - streaming on fine grid (swap buffer required)
%fprintf("Apply streaming on level %d\n",level);
f = streaming(f, dirx, diry, hier{level}.NX, hier{level}.NY, ndir);

%% C-redistribution on level
if(level < num_level)
    %fprintf("Apply C-redistribution on level %d\n",level);
    
    % a cell marker is needed in fine interface cells
    % to mark in which direction the particle is updated in current level
    % rest parts will be updated from next coarser level after streaming
    isFluidOrToFineCellNextLevel = hier{level+1}.gridding==10|hier{level+1}.gridding==100;
    isUpdated = repmat(isFluidOrToFineCellNextLevel,1,9);
    isUpdated1 = streaming(isUpdated, dirx, diry, hier{level+1}.NX, hier{level+1}.NY, ndir);
    isUpdated2 = streaming(isUpdated1, dirx, diry, hier{level+1}.NX, hier{level+1}.NY, ndir);
    isUpdated = (isUpdated1|isUpdated2)&(~isUpdated);
    
    % scan all interface_to_coarse cells
    for i = 2:9
        isUpdatedInterfaceCell = (hier{level+1}.gridding==1)&(~isUpdated(:,i));
        hier{level+1}.f(isUpdatedInterfaceCell,i) =...
            f(hier{level}.interface_to_fine_connectivity(isUpdatedInterfaceCell),i)/4;
    end
    
    % redistribution, operate->interface_cell on coarse grid
    isToFineCell=hier{level}.gridding == 100;
    f(isToFineCell,:) = ...
          hier{level+1}.f(hier{level+1}.interface_to_coarse_connectivity(isToFineCell,1),:)  ...
        + hier{level+1}.f(hier{level+1}.interface_to_coarse_connectivity(isToFineCell,2),:)  ...
        + hier{level+1}.f(hier{level+1}.interface_to_coarse_connectivity(isToFineCell,3),:)  ...
        + hier{level+1}.f(hier{level+1}.interface_to_coarse_connectivity(isToFineCell,4),:) ;
    
end

%% Output
hier{level}.f = f;

%fprintf("Leave level %d\n",level);

%out=hier;

end

