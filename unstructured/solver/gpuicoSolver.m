function gpuicoSolver(IterSettings,PhysicalProperties,mesh,SaveData)
%% Author: Yuxuan Liu

%% Read Data
dt = IterSettings.dt;

%% Read Data
ncells = mesh.ncells;
nfaces = mesh.nfaces;
verts = mesh.verts; % w/o trans
vol = mesh.vol;
areaf = mesh.areaf;
link_bface_to_face = mesh.link_bface_to_face;
owner = mesh.owner;
neighbour = mesh.neighbour;
is_bface = mesh.is_bface;
Sf = mesh.Sf;
gamma = mesh.lf_weighted;
df = mesh.df;

%% Read Data

[...
    nu,...
    U,...
    p,...
    phi,...
    INLET,...
    OUTLET,...
    WALL, ...
    MOVINGWALL,...
    setRef...
    ] = PhysicalProperties{:};

% gpu support
setRef.refCell = gpuArray(setRef.refCell);

%% new feature
fprintf('Start computation ...\n')
tic
[diffOperator,crossDiffOperator] = diffusionOperator(nu, mesh,true);
nodeInterpOperator = cell2nodeinterpOperator(mesh,true);
leastSquareOperator = reconstructOperator(mesh,true);
toc

%
[nodelist,~,~] = find(mesh.link_face_to_node);

%% gpu support
ncells = gpuArray(ncells);
nfaces = gpuArray(nfaces);

vol = gpuArray(vol);
areaf = gpuArray(areaf);
link_bface_to_face = gpuArray(link_bface_to_face);
owner = gpuArray(owner);
neighbour = gpuArray(neighbour);
is_bface = gpuArray(is_bface);
Sf = gpuArray(Sf);
gamma = gpuArray(gamma);
df = gpuArray(df);

nodelist = gpuArray(nodelist);

%% Iteration start point
loop_start_flag = 1;

t = IterSettings.dt:IterSettings.dt:IterSettings.tf;
stepnum = size(t,2);

%% Read saved data

if nargin > 3
    autosave = load (SaveData);
    U = autosave.('U');
    p = autosave.('p');
    phi = autosave.('phi');
    clear savedata
end

ts = tic;
%% Inner Loop start point
for i = 1:stepnum
    %% Assign value for history result
    
    U.old = U.field;
    p_old = p;
    
    %% PISO
    % Compute for cross-diffusion term

    Jf = gpuCrossDiffusion(U.field, nodeInterpOperator, crossDiffOperator, nodelist);

    %% Update matrix
    % Ax = b
    % construct A
    %CrankNicolson = 0.45; % 0 ~ 0.5
    % To improve performance, consider splitting sparse matrix to L,U,D or A, H
    
    if loop_start_flag
        o_diag = gpuArray([owner(~is_bface);owner(~is_bface)]);
        n_diag = gpuArray([neighbour(~is_bface);neighbour(~is_bface)]);
        o_offdiag = gpuArray([owner(~is_bface) ;neighbour(~is_bface)]);
        n_offdiag = gpuArray([neighbour(~is_bface);owner(~is_bface)]);
        
        % construct visc term
        VISC_d = sparse(o_offdiag,o_offdiag,[diffOperator(~is_bface);diffOperator(~is_bface)],ncells,ncells);
        VISC_h = sparse(o_offdiag,n_offdiag,[-diffOperator(~is_bface);-diffOperator(~is_bface)],ncells,ncells);
        VISC = VISC_d + VISC_h; clear VISC_d VISC_h
        
        % pre-allocate space for b term
        b_uv = zeros(ncells,2,'gpuArray');
        b_corr = zeros(ncells,2,'gpuArray');
        
        % pre-allocate space for ddt term
        I = spdiags(vol/dt,0,ncells,ncells);
    end
    
    % construct convection term
    COV = sparse(o_offdiag,n_offdiag,0.5*[phi(~is_bface);-phi(~is_bface)],ncells,ncells);
    
    Tri = COV + VISC;
    
    % construct b
    b_uv(:,1) = accumarray(o_offdiag,[-Jf(~is_bface,1);Jf(~is_bface,1)]);
    b_uv(:,2) = accumarray(o_offdiag,[-Jf(~is_bface,2);Jf(~is_bface,2)]);
    
    % correct boundary conditions
    ifc = link_bface_to_face;
    is_dirchlet = U.boundary.dirchlet;
    is_neumann = U.boundary.neumann;
    
    diag_corr = (diffOperator(ifc) - phi(ifc)/2).*is_dirchlet + (phi(ifc)/2).*is_neumann;
    
    Tri = Tri + sparse(owner(ifc),owner(ifc),diag_corr,ncells,ncells);
    
    b_corr(:,1) = accumarray(owner(ifc), - ...
        (U.boundary.field(:,1).*(-diffOperator(ifc) + phi(ifc)) - 2.*Jf(ifc,1)).*is_dirchlet...
        + Jf(ifc,1).*is_neumann,[ncells 1],[],0);
    b_corr(:,2) = accumarray(owner(ifc), - ...
        (U.boundary.field(:,2).*(-diffOperator(ifc) + phi(ifc)) - 2.*Jf(ifc,2)).*is_dirchlet...
        + Jf(ifc,2).*is_neumann,[ncells 1],[],0);
    
    b_uv = b_uv + b_corr;
    
    % (Implicit)dt & dpidxi
    DDT = I - 0.5*Tri;
    b_uv = b_uv + DDT*U.field;
    Tri = DDT + Tri;
    
    % (Linear multistep)store old flux
    %if loop_start_flag,f = b_uv - Tri*U.field;end
    %f_0 = f;f = b_uv - Tri*U.field;
    
    %% Solve U equation

    %Linear multistep
    %U.field = U.field + dt*(1.5*f - 0.5*f_0)./vol;
    
    % pbicgstab matrix solver
    [U.field(:,1),~,RELRES,ITER] = bicgstab(Tri,b_uv(:,1),1e-5,1000,[],[],U.field(:,1));
    fprintf('DILUPBiCGStab:Solving Ux ,residual: %e , iter: %i \n',RELRES,round(ITER));
    [U.field(:,2),~,RELRES,ITER] = bicgstab(Tri,b_uv(:,2),1e-5,1000,[],[],U.field(:,2));
    fprintf('DILUPBiCGStab:Solving Uy ,residual: %e , iter: %i \n',RELRES,round(ITER));
    
    %% Inner loop
    final_inner_iter_flag = 1;
    for inner_iter_flag = 1:final_inner_iter_flag
        
        %% Prepare to solve pressure correction equation
        
        if loop_start_flag
            rAU = gpuArray(dt);
            rAUf = zeros(nfaces,1,'gpuArray');
            rAUf(~is_bface) = rAU.*areaf(~is_bface)./df(~is_bface); % interpolate rAU to faces, df account for gradp
        end
        
        %% Update guessed surface normal flux
        
        phi_guess = gpuInterpPhi(U, owner, neighbour, is_bface, Sf, link_bface_to_face, nfaces); % 05/06 convergence becomes bad
        
        %% Adjust phi
        
        if ~isempty(INLET)
            if loop_start_flag
                inflow = link_bface_to_face(INLET);
                outflow = link_bface_to_face(OUTLET);
                inflow = gpuArray(inflow);
                outflow = gpuArray(outflow);
            end
            if abs(sum(phi(outflow))) > 1e-06 % edited @ 12/4
                phi(outflow) = -phi(outflow)*(sum(phi(inflow))/sum(phi(outflow)));
            end
            if any(phi(outflow) < 0)
                phi(outflow) = min(0,phi(outflow)); % block reversed flow
                fprintf('Reversed flow appears in outlet zone! \n');
            end
        end
        %fprintf('Adjust outflow, mass outflow = %e.
        %\n',sum(phi(outflow))); %debug
        
        %% Build pressure correction equation
        % Ax = b
        % construct A
        if loop_start_flag
            tmp_gradp_o = sparse(o_diag,o_offdiag,[rAUf(~is_bface);-rAUf(~is_bface)],ncells,ncells);
            tmp_gradp_n = sparse(n_diag,n_offdiag,[rAUf(~is_bface);-rAUf(~is_bface)],ncells,ncells);
            
            Tri_p = tmp_gradp_o + tmp_gradp_n;
            
            local_tri_p = gather(Tri_p);% caution: gather
            alpha = max(sum(abs(local_tri_p),2)./diag(local_tri_p))-2; 
            L_p = ichol(local_tri_p, struct('type','ict','droptol',1e-04,'diagcomp',alpha));

            clear local_* tmp_*
        end
        
        % construct b
        b_p = accumarray([o_offdiag;owner(~~is_bface)],[-phi_guess(~is_bface);phi_guess(~is_bface);-phi_guess(~~is_bface)]);
        
        %% correct p boundary condition
        
        if ~isempty(INLET)
            if loop_start_flag
                cell_outflow = owner(link_bface_to_face(OUTLET))';
                
                Tri_p = gather(Tri_p); % caution: gather
                for icell_o = cell_outflow
                    Tri_p(icell_o,:) = 0;
                    Tri_p(icell_o,icell_o) = 1;
                end
                Tri_p = gpuArray(Tri_p);
            end
            
            b_p(cell_outflow) = 0;
            
        else
            
            if loop_start_flag
                Tri_p = gather(Tri_p); % caution: gather
                Tri_p(setRef.refCell,:) = 0;
                Tri_p(:,setRef.refCell) = 0; % or not
                Tri_p(setRef.refCell,setRef.refCell) = 1;
                Tri_p = gpuArray(Tri_p);
            end
            
            b_p(setRef.refCell) = 0;
            
        end

        %% Solve pressure correction equation
        
        tol = [0.01 1e-06];
        if loop_start_flag, local_tri_p = gather(Tri_p);end
        
        % direct solve
        local_p = local_tri_p\gather(b_p);
        
        % PCG solver
        %[local_p,~,RELRES,ITER] = pcg(local_tri_p,gather(b_p),tol(end),1000,L_p,L_p',gather(p));
        %fprintf('DICPCG:Solving p ,residual: %e , iter: %4d \n',RELRES,sum(ITER));
        
        %GMRES solver
        %[local_p,~,RELRES,ITER] = gmres(local_tri_p,gather(b_p),10,tol(end),20,L_p,L_p',gather(p));
        %fprintf('DICGMRES:Solving p ,residual: %e , iter: %4d \n',RELRES,sum(ITER));
           
        %send back to gpu
        p = gpuArray(local_p);
        
        %% Correct phi
        
        phi(~is_bface) = phi_guess(~is_bface) ...
            - rAUf(~is_bface).*(p(neighbour(~is_bface)) - p(owner(~is_bface))); % sn.*lf -> df 05/28/2018
        
        phi(~~is_bface) = phi_guess(~~is_bface); % for stability
        
        %% Reconstruct pressure gradient

        gradp = gpuReconstruct(gamma, p, leastSquareOperator, o_offdiag, owner, neighbour, is_bface);
        
        %% Update p,U
        
        U.field = U.field - rAU.*gradp;
        %p = p - p(setRef.refCell) + setRef.refValue; % remove alpha_p 05/28/2018
        
        
    end
    
    %% Residual
    p_change = gather(norm(p - p_old, 2));
    
    % Residual Monitor
    if loop_start_flag
        clf,h = figure(1);
        set(h,'Units','pixels','Position',[1 1 1080 960]);
        TRI = delaunay(verts(:,1),verts(:,2));
        obj_1 = patch('Faces',TRI,'Vertices',verts,'FaceColor','interp','EdgeColor','None');
        %obj_1 = scatter(verts(:,1),verts(:,2),[],verts(:,1),'.');
        %hold on, [drawx,drawy] = drawPolygon; fill(drawx,drawy,[.7 .7 .7]), hold off
        obj_title =  title(['Umag, current time: ',num2str(t(i))]);
        colorbar,xlabel('X'),ylabel('Y'),axis equal,view(2)%,colormap jet
        if ~exist('image_output','dir'), mkdir image_output, end
    end
    
    if mod(i,IterSettings.plot_period) == 0||i == stepnum
        obj_1.CData = gather(mag(U.field,nodeInterpOperator));
        obj_title.String = ['Umag, current time: ',num2str(t(i))];
        drawnow
        saveas(h,[ '.\image_output\result' '.' sprintf('%06d', i) '.png'])
        fprintf('Current figure has been saved.\n');
    end
    
    fprintf('Current time: %f, pressure relTol: %e.\n',t(i),p_change);
    
    % autosave
    if mod(i,IterSettings.autosave_period) == 0||i == stepnum
        save('autosave.mat','U','p','phi')
        fprintf('Current computation has been saved in autosave.m.\n');
    end
    
    %% Something else
    if loop_start_flag, loop_start_flag = 0; end
end
telapsed = toc(ts);
fprintf([num2str(telapsed),' sec elasped for total computation.\n']);


end

%%
function varout = mag(inputArg1,inputArg2)

tmp_var = gpuCell2Verts(inputArg1,inputArg2); % caution: gather
varout = sqrt(tmp_var(:,1).^2 + tmp_var(:,2).^2);

end

function [x,y] = drawPolygon

origin = [0 0];
D = 1;
t = linspace(0,2*pi,100);
x = D/2*cos(t) + origin(1);
y = D/2*sin(t) + origin(2);

x = x(:); y = y(:);
end