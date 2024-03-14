function d = tsdEquilateralTriangle( px,py )
loc = [0.5;0.5];
scale = 0.2;
p  = ([px;py]-loc)/scale; 
    k = sqrt(3.0);
    p(1) = abs(p(1)) - 1.0;
    p(2) = p(2) + 1.0/k;
    if( (p(1)+k*p(2))>0.0 ), p = [p(1)-k*p(2),-k*p(1)-p(2)]/2.0;end
    p(1) = p(1) - min(max( p(1), -2.0), 0.0 );
    %p = p * 0.1;
    d = -norm(p,2)*sign(p(2));

end