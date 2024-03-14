function [dc,IA,IC,p1,p2,p3,p4,a1,a2,a3,a4,d] = ImmersedBoundary2D(X,Y,meshsize)
rn = rand(300,3);
[d,marker, gp(:,1),gp(:,2)] = arrayfun(@(x,y,m)(calcDist2D(x,y,m,300, rn)),X(:),Y(:),meshsize+0*X(:));
Xmin = min(X(:));Xmax = max(X(:));
Ymin = min(Y(:));Ymax = max(Y(:));
IA = find(marker==1 & all(isfinite(gp),2) & all(gp>[Xmin Ymin],2)  & all(gp<[Xmax Ymax],2));
IC = find(marker==2);
dc = d(IA);
gp = gp(IA,:);

for i = 1:numel(IA)
    pmin(i,1) = find(X<gp(i,1)&Y<gp(i,2),1,'last');
    pmax(i,1) = find(X>gp(i,1)&Y>gp(i,2),1,'first');
end

p1 = pmin;
p2 = pmin+1;
p3 = pmax-1;
p4 = pmax;

[a1,a2,a3,a4] = arrayfun(@bilinearCoeffs,gp(:,1),gp(:,2),X(p1),Y(p1),X(p4),Y(p4));

% Flux correction coefficients
dc = (abs(dc) + 1.75*meshsize)/(1.75*meshsize);
end

function [a1,a2,a3,a4] = bilinearCoeffs(x,y,x1,y1,x4,y4)
s1 = (x4-x)/(x4-x1);
s2 = (y4-y)/(y4-y1);

a1 = s1*s2;
a2 = (1-s1)*s2;
a3 = (1-s2)*s1;
a4 = (1-s1)*(1-s2);
end

function [d, type, gpx, gpy] = calcDist2D(x, y, meshsize, n, rn)
%METHOD1 Summary of this method goes here
%   Detailed explanation goes here
%f = @(p)(sdNoisy(p, n, rn));
f = @sdCircle;
d = f([x;y]);
n = calcNormal2D(f,[x;y]);

if(d > 0)
    type = 0; % fluid cell
else
    if(d < -2*sqrt(2)*meshsize)
        type = 2; % object cell
    else
        type = 1; % interface cell
    end
end

gpx = x + n(1)*(abs(d)+1.75*meshsize);
gpy = y + n(2)*(abs(d)+1.75*meshsize);

end

function n = calcNormal2D( f,p )
eps = 0.0001; % or some other value
n =  normalize([f(p+[eps;0]) - f(p-[eps;0]);...
    f(p+[0;eps]) - f(p-[0;eps])],'norm');
end

function n = calcNormal3D( f, p )
h = 0.0001; % replace by an appropriate value
%k = [1;-1];
%return normalize( k.xyy*f( p + k.xyy*h ) +
%              k.yyx*f( p + k.yyx*h ) +
%              k.yxy*f( p + k.yxy*h ) +
%              k.xxx*f( p + k.xxx*h ) );
n = normalize( [1;-1;-1]*f( p + [1;-1;-1]*h ) + ...
    [-1;-1;1]*f( p + [-1;-1;1]*h ) + ...
    [-1;1;-1]*f( p + [-1;1;-1]*h ) + ...
    [1;1;1]*f( p + [1;1;1]*h ) );
end

% if crashed due to f inlining, try this.
%{
function n = calcNormalAlt( pos )
            % inspired by tdhooper and klems - a way to prevent the compiler from inlining map() 4 times
            n = zeros(3,1);
            for i=1:3
                e = 0.5773*(2.0*[(((i+3)>>1)&1),((i>>1)&1),(i&1)]-1.0);
                n = n + e*map(pos+0.0005*e);
                %if( n.x+n.y+n.z>100.0 ) break;end
            end
            n = normalize(n);
        end
%}

function d = sdParabola( p )

k =  0.2; % width
p = p - 1;
p(2) = abs(p(2));
p(1) = p(1)-0.5;
p(2) = p(2)-0.1;

p(1) = abs(p(1));

ik = 1.0/k;
ip = ik*(p(2) - 0.5*ik)/3.0;
q = 0.25*ik*ik*p(1);

h = q*q - ip*ip*ip;
r = sqrt(abs(h));

if(h > 0)
    x = (q+r)^(1.0/3.0) - (abs(q-r))^(1.0/3.0)*sign(r-q);% 1 root
else
    x = 2.0*cos(atan(r/q)/3.0)*sqrt(ip);% 3 roots
end

d = norm(p-[x;k*x*x]) * sign(p(1)-x);
end

function d = sdCircle(p)
p = p - [1.5;1.5];
c = [0.0;0.0];
r = 0.2;
d = norm(p-c) - r;
end

function d = sdManyCircle(p)
p = p - [1.0;1.0];
c = [0.0;0.0];
r = 0.3;
d = norm(p-c) - r;
d = d + 0.5*r*randn();
end

function d = sdNoisy(p,n, rn)
% random position and size
%d = randn()+3.0;
d = Inf;
for i = 1:n
c = [0.0;0.0]+6*[rn(i,1);rn(i,2)];
r = 0.02+0.01*rn(i,3);
%d = min(d, norm(p-c) - r);
dir = deg2rad(30);
pi = [cos(dir) -sin(dir);...
    sin(dir) cos(dir)]*p;

di = abs(pi-c)-r;
di = norm(max(di,0.0),2) + min(max(di(1),di(2)),0.0);
d = min(d, di);
end

p = p - [3.0;1.5];
c = [-0.65;-0.4];
r = 0.4;
di = -min(0.0,norm(p-c) - r);
di(di==0) = -Inf;
d = max(d, di);

end

function d = sdEquilateralTriangle( p )

loc = [1.0;1.0];
scale = [0.2;0.1];
p  = (p-loc)./scale; 
dir = deg2rad(30);
p = [cos(dir) -sin(dir);...
    sin(dir) cos(dir)]*p;

    k = sqrt(3.0);
    p(1) = abs(p(1)) - 1.0;
    p(2) = p(2) + 1.0/k;
    if( (p(1)+k*p(2))>0.0 ), p = [p(1)-k*p(2),-k*p(1)-p(2)]/2.0;end
    p(1) = p(1) - min(max( p(1), -2.0), 0.0 );
    %p = p * 0.1;
    d = -norm(p,2)*sign(p(2));

end

function d = sdBox( p )

loc = [1.0;1.0];
scale = [0.2;0.1]/4;
p  = (p-loc)./scale; 

b = 1.0;
d = abs(p)-b;
d = norm(max(d,0.0),2) + min(max(d(1),d(2)),0.0);

end

function d = sdUnevenCapsule( p )
loc = [1.0;1.0];
scale = 0.2;
p  = (p-loc)/scale; 
dir = deg2rad(90);
p = [cos(dir) -sin(dir);...
    sin(dir) cos(dir)]*p;

    r1 = 0.4;
    r2 = 0.15;
    h = 1.0;
    p(1) = abs(p(1));
    b = (r1-r2)/h;
    a = sqrt(1.0-b*b);
    k = dot(p,[-b;a]);
    if( k < 0.0 ), d = norm(p) - r1;return;end
    if( k > a*h ), d = norm(p-[0.0;h]) - r2;return;end
    d = dot(p, [a;b])  - r1;
end

function d= sdArc( p)

% animation
p = p - 1;
iTime = 10.0;
ta = 3.14*(0.5+0.5*cos(iTime*0.52+2.0));
tb = 3.14*(0.5+0.5*cos(iTime*0.31+2.0));
rb = 0.15*(0.5+0.5*cos(iTime*0.41+3.0));

% distance
sca = [sin(ta) cos(ta)];
scb = [sin(tb) cos(tb)];
ra = 0.5;

p = [sca(1) sca(2);-sca(2) sca(1)]*p;
p(1) = abs(p(1));
T = (scb(2)*p(1)>scb(1)*p(2));
k = T* dot(p,scb) + (1-T)* norm(p);
d= sqrt( dot(p,p) + ra*ra - 2.0*ra*k ) - rb;
end

function d= sdHorseshoe( p )

% animation
%p = 2*p-1;
p = p - 1;
dir = deg2rad(90);
p = [cos(dir) -sin(dir);...
    sin(dir) cos(dir)]*p;
%p = p - [0;-0.1];

iTime = 2.57;
t = 3.14* (0.3+0.3*cos(iTime*0.5));
c = [cos(t);sin(t)];
r = 0.5;
w = [0.750;0.25].*(0.5+0.5*cos(iTime*[0.7;1.1]+[0.0;2.0]));

% distance
p(1) = abs(p(1));
l = norm(p);
p = [-c(1)  c(2);  c(2)  c(1)]*p;
T = double(p > 0.0);
p = [T(2)*p(1)+(1-T(2))*l*sign(-c(1));T(1)*p(2)+(1-T(1))*l ];
p = [p(1);abs(p(2)-r)]-w;
d = norm(max(p,0.0)) + min(0.0,max(p(1),p(2)));
%d = d/2;
end