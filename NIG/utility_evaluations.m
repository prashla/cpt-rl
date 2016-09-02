function [pos_utility] = utility_evaluations(alpha, beta, mu, delta, b , N)
%
% Calculating cpt value using gaussian quadrature
%
% alpha:shape parameter;
% beta: skewness parameter;
% mu:location parameter;
% delta:scale parameter;
% a : lower bound of the interval
% b : upper bound of the interval
% N : number of nodes in Gaussian quadrature
    


    nigpdfu = @(u)(nigpdf(u, alpha, beta, mu, delta));
    
    function z = pos_integrand(u)
        z = nigpdfu(u) *  u;
    end     
    
    
    % set up the nodes and weights for Gaussian quadrature:
    [pos_n, pos_w] = lgwt(N,0,b);
    pos_qua_evaluations = zeros(size(pos_n));
    for i = 1:N;
       pos_qua_evaluations(i) = pos_integrand(pos_n(i));
    end
    pos_utility = sum(pos_w.* pos_qua_evaluations);
    
    % negative part
    

    
end









