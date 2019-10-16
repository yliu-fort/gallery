function [ WBC,EBC,SBC,NBC,INT ] = CatchBoundary( bdrange,CatchFactor,GeometryProperties)
%   bdrange = [xmin ymin;xmax ymax]
%   Strongly suggest to have square outside boundaries attach to the box
%   which bdrange defined. any boundary face does not fit to the edge of
%   the box will be recognized as INT boundary.
%   CatchFactor: distance from a face to boundary, use to
%   tell if a face belongs to a specific boundary.
%   Further implementation may include functions.

%% Read data
%{
                            [...
                            ~, ~, nbfaces, ~, ...
                            ~, FC, ~, ...
                            ~, ~, ~, ~, ...
                            ~, ~, ~, ...
                            ~, ...
                            ~, ...
                            ~, ...
                            ~, ...
                            ~, ...
                            ~, ...
                            link_bface_to_face, ...
                            ] =   GeometryProperties{:};
                        %}

%%
nbfaces = GeometryProperties.nbfaces;
FC = GeometryProperties.FC;
link_bface_to_face = GeometryProperties.link_bface_to_face;

%% Some default settings
% further implementation required to avoid misread the settings

if ischar(bdrange) ~= 0
if strfind(bdrange,'cavity')
    bdrange = [0 0;1 1];
elseif strfind(bdrange, 'cylinder')
    bdrange = [-3 -3;12 3];
elseif strfind(bdrange, 'pitzDaily')
    bdrange = [-20.6 -25.4;290 25.4]*0.001;
end
end

%% If not default

if ischar(bdrange) == 0
xmin = bdrange(1,1); ymin = bdrange(1,2);
xmax = bdrange(2,1); ymax = bdrange(2,2);
end

%% Capture Boundaries

WBC = [];EBC = [];SBC = [];NBC = []; INT = [];
for ibfc = 1:nbfaces
    not_any_bc = 1;
        if  (abs(FC(link_bface_to_face(ibfc),1) - xmin) < CatchFactor)&&not_any_bc
            WBC = [WBC;ibfc];
            not_any_bc = 0;
        end
        if  (abs(-FC(link_bface_to_face(ibfc),1) + xmax) < CatchFactor)&&not_any_bc
            EBC = [EBC;ibfc];
            not_any_bc = 0;
        end
        if  (abs(FC(link_bface_to_face(ibfc),2) - ymin) < CatchFactor)&&not_any_bc
            SBC = [SBC;ibfc];
            not_any_bc = 0;
        end
        if  (abs(-FC(link_bface_to_face(ibfc),2) + ymax) < CatchFactor)&&not_any_bc
            NBC = [NBC;ibfc];
            not_any_bc = 0;
        end
        if  not_any_bc
            INT = [INT;ibfc];
            %fprintf(['Boundary face ',num2str(ibfc), ' has been classified as int boundary.\n'])
        else
            %fprintf(['Boundary face ',num2str(ibfc), ' has been classified.\n'])
        end
        
    
end

%% Check if all boundaries are assigned
Assertion3 = size(WBC,1) + size(EBC,1) + size(SBC,1) + size(NBC,1) + size(INT,1) == nbfaces;
assert(Assertion3,'BoundaryCatcherError: Some boundaries are failed to be classified.' )
fprintf('All boundaries has been settled.\n')

%% Sort the boundary
% orginal order of boundary face is random, so it is necessary to sort it in
% order to apply boundary condition.
% W/EBC starts from bottom to top.
% N/SBC starts from left to right.
% INT starts from random node and loops clockwise.

    [~,index] = sort(FC(link_bface_to_face(WBC),2));
    WBC = WBC(index);
    [~,index] = sort(FC(link_bface_to_face(EBC),2));
    EBC = EBC(index);
    [~,index] = sort(FC(link_bface_to_face(NBC),1));
    NBC = NBC(index);
    [~,index] = sort(FC(link_bface_to_face(SBC),1));
    SBC = SBC(index);
    CC_vert = mean(FC(link_bface_to_face(INT),:));
    INT_vecx = FC(link_bface_to_face(INT),1) - CC_vert(1,1);
    INT_vecy = FC(link_bface_to_face(INT),2) - CC_vert(1,2);
    ANG_INT = atan2d(INT_vecy, INT_vecx);
    [~,index] = sort(ANG_INT);
    INT = INT(index);

end

