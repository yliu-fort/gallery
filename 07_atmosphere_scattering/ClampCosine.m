function result = ClampCosine(mu)
    result = min(max(mu, -1.0), 1.0);
end
