% Produce the sequence of cpt-utilities and normal utilities
alphas = 3:0.005:6;
mu = 2; 
delta = 3;
a = -8;
b = 8;
N = 30;
beta = 2;
cpt_values = zeros(size(alphas));
for i = 1:numel(alphas)
    cpt_values(i) = cpt_nig(alphas(i), beta, mu, delta, a, b , N);
end

utility_values = zeros(size(alphas));
for i = 1:numel(alphas)
    utility_values(i) = utility_nig(alphas(i), beta, mu, delta, a, b , N);
end




