function a_GSD = DF_SGD_function(a_GSD,y_SGD,x)
    mu = 0.01; % Step
    a_GSD = a_GSD + mu*conj(y_SGD)*x; % [Order x 1]
end