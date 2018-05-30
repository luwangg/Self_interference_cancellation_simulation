function A = DF_SGD_function(A,r_r,r_t)
    mu = 0.0005; % Step
    A = A - mu*r_r*r_t'; % [1 x Order]
end