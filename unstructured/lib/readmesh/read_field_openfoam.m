function output = read_field_openfoam(input,debug)
%% input management

if nargin < 2
    debug = 0;
else
    assert(any(find([0 1 2 3] == debug)),'Fatal error: debug level can only be set to 1, 2 or 3.\n');
    if ~debug == 0, fprintf('readOFdata_sup %d level debug mode on.\n',debug);end
end

%% read Information from headers
if debug == 1,disp('read_field_openfoam > filename:');disp(input);end % level 1 debug
fileID = fopen(input,'r');
Intro = textscan(fileID,'%s',25,'Delimiter','\n'); % Skip the header lines
% Get object information

str = strcat(Intro{:}{:});
ptn = 'class\s*(?<class>\w+).*object\s*(?<object>\w+).*?(?<length>\d+)[(]';
info = regexpi(str,ptn,'names');
info.length = str2num(info.length);

if debug == 1,disp('read_field_openfoam > info:');disp(info);end % level 1 debug

% Get headerlines before data
for i = 1:25, if str2num(Intro{1}{i}) == info.length, Headerline = i + 1 ;break, end, end
fseek(fileID,0,-1); % reset cursor

%% check data type
switch info.class
    case 'labelList'
        ptn = '%d';
        fieldType = {'value'};
    case 'faceList'
        ptn = '%s';
        fieldType = {'id','index'};
    case 'polyBoundaryMesh'
        ptn = '%s';
        fieldType = {'patch'};
    case 'vectorField'
        ptn = '(%f %f %f)';
        fieldType = {'x','y','z'};
    case 'volVectorField'
        ptn = '(%f %f %f)';
        fieldType = {'x','y','z'};
    case 'volScalarField'
        ptn = '%f';
        fieldType = {'value'};
    case 'surfaceScalarField'
        ptn = '%f';
        fieldType = {'value'};
    case 'volTensorField' % not implimented
        ptn = '(%f %f %f %f %f %f %f %f %f)';
        %Headers = {[info.object 'x'],[info.object 'y'],[info.object 'z']};
        fprintf('Tensor field processing not implemented.\n')
    otherwise
        fprintf('Unrecognized data type.Please add corresponding info in _sup.m.\n');
end

%% read data
indata = textscan(fileID,ptn,'Delimiter','\n', 'headerLines', Headerline);

% if unique data structure, apply modifications
if strcmp(info.class,'faceList'),indata = feval(@facelist,indata,info);end
if strcmp(info.class,'polyBoundaryMesh'),indata = feval(@polyBoundaryMesh,indata,info);end
if strcmp(info.class,'labelList'),indata = feval(@cpp2mat,indata,'labelList');end

%% Output struct
output.class = info.class;
output.object = info.object;
output.length = info.length;
output.field = cell2struct(indata,fieldType,2);

if debug == 1,disp('read_field_openfoam > output:');disp(output);end % level 1 debug

fclose(fileID);
end
%%
function output = facelist(str,info)
str = str{1}(1:info.length);
rawlist = split(str,["(",")"]);
rawlist = rawlist + " ";

id = rawlist(:,1);
id = str2num([id{:}]);
nodelist = rawlist(:,2);
nodelist = str2num([nodelist{:}]);

i = feval(@cpp2mat,cat(2,nodelist)); % c++ start from 0, do i = i + 1 for matlab indexing
j = repelem(1:size(id,2),id);

index = sparse(i,j,true,max(i),max(j));
output = {id' index};

end

function output = polyBoundaryMesh(str,info)
str = strcat(str{:}{:});
rawlist = split(str,["{","}"]);
rawlist = reshape(rawlist(1:2*info.length),2,[]);

for i = rawlist
    A = regexpi(i{2},'(?<attribute>\w+)\s+(?<value>(\d|\w)+)','names');
    output.(i{1}) = cell2struct({A.value}',{A.attribute}');
end

output = {output};

end

function output = cpp2mat(index,type)
if nargin < 2, type = 'double';end
switch type
    case 'double'
        output = index + 1;
    case 'labelList'
        output = {index{:} + 1};
    otherwise
        fprintf('cpp2mat: not a supported data type. function will return null.\n')
        return
end
end