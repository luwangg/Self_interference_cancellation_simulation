function [a_RLS,P] = DF_RLS_function(a_RLS,y_RLS,x,P)
    lambda = 0.9999; % Forgetting factor
    g = (lambda^-1*P*x) / (1+x'*lambda^-1*P*x);
    P = (1/lambda)*(P-g*x'*P);
    a_RLS = a_RLS + g*conj(y_RLS);
end