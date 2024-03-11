function out = smoothstep (edge0, edge1, x) 
   % Scale, and clamp x to 0..1 range
   x = clamp((x - edge0) ./ (edge1 - edge0), 0, 1);

   out = x .* x .* (3.0 - 2.0 * x);

end

function x = clamp(x, lowerlimit, upperlimit) 
    x = min(max(x, lowerlimit), upperlimit);
end