function result = ClampRadius(r, Rg, Rt)
    result = min(max(r, Rg), Rt);
end