function icoSolver(IterSettings,PhysicalProperties,mesh,SaveData)
%% Author: Yuxuan Liu

%% Read Data
dt = IterSettings.dt;

%% Read Data
ncells = mesh.ncells;
nfaces = mesh.nfaces;
nbfaces = mesh.nbfaces;
verts = mesh.verts;
vol = mesh.vol;
areaf = mesh.areaf;
link_bface_to_face = mesh.link_bface_to_face;
owner = mesh.owner;
neighbour = mesh.neighbour;
is_bface = mesh.is_bface;

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

%% new feature
fprintf('Start computation ...\n')
tic
[diffOperator,crossDiffOperator] = diffusionOperator(nu, mesh);
nodeInterpOperator = cell2nodeinterpOperator(mesh);
leastSquareOperator = reconstructOperator(mesh);
toc

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
    %U.old.y = U.y;
    p_old = p;
    
    %% PISO
    % Compute for cross-diffusion term

    Jf = ParallelCrossDiffusionTerm(U.field, nodeInterpOperator, crossDiffOperator, mesh);

    %% Update matrix
    % Ax = b
    % construct A
    %CrankNicolson = 0.45; % 0 ~ 0.5
    % To improve performance, consider splitting sparse matrix to L,U,D or A, H
    
    if loop_start_flag
        o_diag = [owner(~is_bface);owner(~is_bface)];
        n_diag = [neighbour(~is_bface);neighbour(~is_bface)];
        o_offdiag = [owner(~is_bface) ;neighbour(~is_bface)];
        n_offdiag = [neighbour(~is_bface);owner(~is_bface)];
        
        % construct visc term
        %VISC_o = sparse(o_diag,o_offdiag,[diffOperator(~is_bface);-diffOperator(~is_bface)],ncells,ncells);
        VISC_d = sparse(o_offdiag,o_offdiag,[diffOperator(~is_bface);diffOperator(~is_bface)],ncells,ncells);
        VISC_h = sparse(o_offdiag,n_offdiag,[-diffOperator(~is_bface);-diffOperator(~is_bface)],ncells,ncells);
        VISC = VISC_d + VISC_h; clear VISC_d VISC_h
        
        % pre-allocate space for b term
        b_uv = zeros(ncells,2);
        
        % pre-allocate space for ddt
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
    
    Tri = Tri + sparse(owner(ifc),owner(ifc),(diffOperator(ifc) - phi(ifc)/2).*is_dirchlet + (phi(ifc)/2).*is_neumann,ncells,ncells);% remove '+ diffOperator(ifc)' from neumann
    
    b_corr(:,1) = accumarray(owner(ifc), - ...
        (U.boundary.field(:,1).*(-diffOperator(ifc) + phi(ifc)) - 2.*Jf(ifc,1)).*is_dirchlet ...
        + Jf(ifc,1).*is_neumann,[ncells 1],[],0);
    b_corr(:,2) = accumarray(owner(ifc), - ...
        (U.boundary.field(:,2).*(-diffOperator(ifc) + phi(ifc)) - 2.*Jf(ifc,2)).*is_dirchlet ...
        + Jf(ifc,2).*is_neumann,[ncells 1],[],0);
    
    %b_uv(owner(ifc),:) = b_uv(owner(ifc),:) - (U.boundary.field.*(-diffOperator(ifc) + phi(ifc)) - 2.*Jf(ifc,:)).*is_dirchlet + Jf(ifc,:).*is_neumann;% edited @ 12/5
    b_uv = b_uv + b_corr;
    
    % (Implicit)dt & dpidxi
    DDT = I - 0.5*Tri;
    b_uv = b_uv + DDT*U.field;
    %b_v = b_v + DDT*U.field(:,2);
    Tri = DDT + Tri;
    
    % (Linear multistep)store old flux
    %if loop_start_flag,f.x = b_u - Tri*U.x;f.y = b_v - Tri*U.y;end
    %f_0 = f;f.x = b_u - Tri*U.x;f.y = b_v - Tri*U.y;
    
    %% Solve U equation
    
    % direct solve
    %U.x = Tri\b_u;U.y = Tri\b_v;
    
    %Linear multistep
    %U.x = U.x + dt*(1.5*f.x - 0.5*f_0.x)./vol;U.y = U.y + dt*(1.5*f.y - 0.5*f_0.y)./vol;
    
    % pbicgstab matrix solver
    [il,iu] = ilu(Tri);
    [U.field(:,1),~,RELRES,ITER] = bicgstab(Tri,b_uv(:,1),1e-6,1000,il,iu);
    fprintf('DILUPBiCGStab:Solving Ux ,residual: %e , iter: %i \n',RELRES,round(ITER));
    [U.field(:,2),~,RELRES,ITER] = bicgstab(Tri,b_uv(:,2),1e-6,1000,il,iu);
    fprintf('DILUPBiCGStab:Solving Uy ,residual: %e , iter: %i \n',RELRES,round(ITER));
    
    %% Inner loop
    final_inner_iter_flag = 1;
    for inner_iter_flag = 1:final_inner_iter_flag
        
        %% Prepare to solve pressure correction equation
        
        if loop_start_flag
            
            for icell = 1:ncells
                %rAU(icell,1) = vol(icell,1)/Tri(icell,icell);
                rAU(icell,1) = dt;
            end
            
            rAUf = zeros(nfaces,1);
            rAUf(~is_bface) = 0.5*(rAU(owner(~is_bface)) + rAU(neighbour(~is_bface))).*areaf(~is_bface)./mesh.df(~is_bface); % interpolate rAU to faces, df account for gradp
        end
        
        %% Update guessed surface normal flux
        
        phi_guess = ParallelSurfNormFlux(U, mesh); % 05/06 convergence becomes bad
        
        %% Adjust phi
        
        if ~isempty(INLET)
            inflow = link_bface_to_face(INLET);
            outflow = link_bface_to_face(OUTLET);
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
            gradp_o = sparse(o_diag,o_offdiag,[rAUf(~is_bface);-rAUf(~is_bface)],ncells,ncells);
            gradp_n = sparse(n_diag,n_offdiag,[rAUf(~is_bface);-rAUf(~is_bface)],ncells,ncells);
            
            Tri_p = gradp_o + gradp_n; clear gradp_o gradp_n
            
            [C1,C2] = precondition(Tri_p,1e-03,'iCHOL');

        end
        
        % construct b
        b_p = accumarray([o_offdiag;owner(~~is_bface)],[-phi_guess(~is_bface);phi_guess(~is_bface);-phi_guess(~~is_bface)]);
        
        %% correct p boundary condition

        if ~isempty(INLET)
            if loop_start_flag
                cell_outflow = owner(link_bface_to_face(OUTLET))';
            end
            
            for icell_o = cell_outflow
                if loop_start_flag,Tri_p(icell_o,:) = 0; Tri_p(icell_o,icell_o) = 1;end
                b_p(icell_o,1) = 0;
            end
            
        else
            
            if loop_start_flag,Tri_p(setRef.refCell,:) = 0;
                Tri_p(:,setRef.refCell) = 0; % or not
                Tri_p(setRef.refCell,setRef.refCell) = 1;
            end
            
            b_p(setRef.refCell) = 0;
            
        end

        %% Solve pressure correction equation
        tol = [0.01 1e-06];
        % direct solve
        %p = Tri_p\b_p;
        
        % PCG solver
        [p,~,RELRES,ITER] = pcg_solve(Tri_p,b_p,tol(end),1000,C1,C2,p);
        fprintf('DICPCG:Solving p ,residual: %e , iter: %4d \n',RELRES,ITER);
        
        %% Correct phi
        
        phi(~is_bface) = phi_guess(~is_bface) ...
            - rAUf(~is_bface).*(p(neighbour(~is_bface)) - p(owner(~is_bface))); % sn.*lf -> df 05/28/2018
        
        phi(~~is_bface) = phi_guess(~~is_bface); % for stability
        
        %% Reconstruct pressure gradient

        gradp = reconstruct([], p, leastSquareOperator, o_offdiag, mesh);
        
        %% Update p,U
        
        U.field = U.field - rAU.*gradp;
        %U.field(:,2) = U.field(:,2) - rAU.*gradp(:,2);
        %p = p - p(setRef.refCell) + setRef.refValue; % remove alpha_p 05/28/2018
        
        
    end
    
    %% Residual
    p_change = norm(p - p_old, 2);
    
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
        obj_1.CData = feval(@mag,U.field,nodeInterpOperator);
        obj_title.String = ['Umag, current time: ',num2str(t(i))];
        drawnow
        %saveas(h,[ '.\image_output\result' '.' sprintf('%06d', i) '.png'])
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

tmp_var = ParallelCell2Verts(inputArg1,inputArg2);
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