function [diffOperator,crossDiffOperator,df,lf] = DiffusionTerm( nu,geo )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

%%
for ifc = 1:geo.nfaces
    if geo.is_bface(ifc)
        lf(ifc,:) = geo.FC(ifc,:) - geo.IC(geo.owner(ifc),:);
        %df(ifc,1) = geo.sn(ifc,:)*lf(:,ifc);
    else
        lf(ifc,:) = geo.IC(geo.neighbour(ifc),:) - geo.IC(geo.owner(ifc),:);
        %df(ifc,1) = geo.sn(ifc,:)*lf(:,ifc);
    end
end

df = dot(geo.sn,lf,2);

assert(min(df) > 0, 'normal distance between owner and neighbour can not be negative or 0.\n')

for iface = 1:geo.nfaces
    node_list = find(geo.link_face_to_node(:,ifc));
    %c1 = verts(link_face_to_node(iface,1),:);
    %c2 = verts(link_face_to_node(iface,2),:);
    c1 = geo.verts(node_list(1),:);
    c2 = geo.verts(node_list(2),:);
    tf(iface,:) = c2 - c1;
end

%tf = geo.verts(geo.link_face_to_node(:,2),:) - geo.verts(geo.link_face_to_node(:,1),:);

%%

diffOperator = nu*geo.areaf./df; %Df(gfc,:) 05/28/18 removed '-'
crossDiffOperator = -nu*dot(tf,lf,2)./df; % remove ./geo.areaf.^2.*diffOperator after "dot(tf,lf,2)" 06/03/18


end

