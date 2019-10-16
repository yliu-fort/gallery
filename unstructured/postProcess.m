%% Postprocessing
% will plot U magnitude, pressure and quiver for surface flux.

%%
sol = load('autosave');

U = sol.U;
p = sol.p;
phi = sol.phi;
                        
%% Read Data
IC = GeometryProperties.IC;
FC = GeometryProperties.FC;
verts = GeometryProperties.verts;
sn = GeometryProperties.sn;
link_face_to_node = GeometryProperties.link_face_to_node;
link_cell_to_node = GeometryProperties.link_cell_to_node;
link_bface_to_face = GeometryProperties.link_bface_to_face;

%% Read Data
[...
    nu,...
    ~,...
    ~,...
    ~,...
    INLET,...
    OUTLET,...
    WALL, ...
    MOVINGWALL,...
    setRef...
    ] = PhysicalProperties{:};
                        
%%
TRI = delaunay(verts(:,1),verts(:,2));
F =scatteredInterpolant(IC(:,1),IC(:,2),U.x);
F.Method = 'natural';
u_verts = F(verts(:,1),verts(:,2));
F.Values = U.y;
v_verts = F(verts(:,1),verts(:,2));
F.Values = p;
p_verts = F(verts(:,1),verts(:,2));
%% Plotting u
figure(2)

trisurf(TRI,verts(:,1),verts(:,2),u_verts)%,'EdgeColor','None'
title('u profile')
colorbar
xlabel('X')
ylabel('Y')
axis equal

%% Plotting v
figure(3)

trisurf(TRI,verts(:,1),verts(:,2),v_verts,'EdgeColor','None')
title('v profile')
colorbar
xlabel('X')
ylabel('Y')
axis equal

%% Plotting p
figure(4)

trisurf(TRI,verts(:,1),verts(:,2),p_verts,'EdgeColor','None','FaceColor','Interp')
title('pressure profile')
colorbar
xlabel('X')
ylabel('Y')
axis equal
view(2)
colormap jet

%% airfoil pressure profile

figure(5)
[node_list,~,~] = find(link_face_to_node(:,link_bface_to_face(WALL)));
node_surf = node_list(:);
node_surf_xcoord = verts(node_surf,1);
p_surf = p_verts(node_surf);
plot(node_surf_xcoord,p_surf,'.');

grid on
set(gca,'ydir','reverse')

%% plotting cell u vector and surface normal flux

figure(6)
patch('Faces',TRI,'Vertices',verts,'FaceColor','None');
hold on

quiver(IC(:,1),IC(:,2),U.x,U.y,1,'b')
quiver(FC(:,1),FC(:,2),phi.*sn(:,1),phi.*sn(:,2),1,'color','r');
hold off
title('Cartisian velocity and surface flux at each cell')
colorbar
xlabel('X')
ylabel('Y')
axis equal

%% Plotting u magnitude

figure(7)
trisurf(TRI,verts(:,1),verts(:,2),u_verts,'EdgeColor','None','FaceColor','Interp')
title('U magnitude profile')
colorbar
xlabel('X'),ylabel('Y'),axis equal
view(2)

%%
if strfind(filename,'cavity')
x = IC(:,1); y = IC(:,2);
[xq,yq] = meshgrid(0:1/128:1, 0:1/128:1);
uq = griddata(x,y,U.x,xq,yq,'natural');
vq = griddata(x,y,U.y,xq,yq,'natural');
figure(7)
plot(uq(:,65),yq(:,65),'LineWidth',1.4)
hold on
index_y = yq([129;126;125;124;123;110;95;80;65;59;31;23;14; 10;9;8;1],65);
ghia_u(:,1) = [1.00000;0.84123;0.78871;0.73722;0.68717;0.23151;0.00332;-0.13641;...
    -0.20581;-0.21090;-0.15662;-0.10150;-0.06434;-0.04775;-0.04192;-0.03717;0.00000];
ghia_u(:,2) = [1.00000;0.75837;0.68439;0.61756;0.55892;0.29093;0.16256;...
0.02135;-0.11477;-0.17119;-0.32726;-0.24299;-0.14612;-0.10338;-0.09266;...
-0.08186;0.00000];
ghia_u(:,3) = [1.00000;0.65928;0.57492;0.51117;0.46604;0.33304;0.18719;0.05702;-0.06080;...
-0.10648;-0.27805;-0.38289;-0.29730;-0.22220;-0.20196;-0.18109;0.00000];
% Ghia 1 = Re100,2 = Re400, 3 = Re1000
plot(ghia_u(:,3),index_y,'--o','LineWidth',1.4)
xlabel('u')
ylabel('y')
title('u velocity @x = 0.5, Re = 1000')
legend('My result','Ghia')
hold off

figure(8)
ghia_v(:,1) = [0.00000;-0.21388;-0.27669;-0.33714;-0.39188;-0.51550;-0.42665;-0.31966;
0.02526;0.32235;0.33075;0.37095;0.32627;0.30353;0.29012;0.27485;0.00000];
index_x = xq(65,[129;125;124;123;122;117;111;104;65;31;30;21;13;11;10;9;1]);
plot(xq(65,:),vq(65,:),'LineWidth',1.4)
hold on
plot(index_x,ghia_v(:,1),'--o','LineWidth',1.4)
xlabel('x')
ylabel('v')
title('v velocity @y = 0.5, Re = 1000')
legend('My result','Ghia')
hold off
end


%%
[xq,yq] = meshgrid(-3:.01:12, -3:.01:3);
uq = griddata(IC(:,1),IC(:,2),U.x,xq,yq,'natural');
vq = griddata(IC(:,1),IC(:,2),U.y,xq,yq,'natural');
streamslice(xq,yq,uq,vq,5,'noarrows');axis equal