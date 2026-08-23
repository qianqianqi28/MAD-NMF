function ang = spectral_angle(x, y)
    ang = acos( (x' * y) / (norm(x) * norm(y)) );
end