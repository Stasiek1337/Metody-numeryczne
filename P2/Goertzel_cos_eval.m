function [p_eval] = Goertzel_cos_eval(p, t_eval)
% funkcja przyjmuje wielomian trygonometryczny p ("cosinusowy"!) oraz punkt
% t_eval i zwraca wartość tego wielomianu w punkcie t_eval - p_eval 
% wyznaczonego za pomocą algorytmu Goretzela

% przeprowadzamy prodecurę Goertzela

a_0 = p(1);
N = max(length(p)) - 1;

b = zeros(1, N + 2);
b(N + 2) = 0;
b(N + 1) = p(N + 1);

p_daszek = 2 * cos(t_eval);
q_daszek = -1;

for i = (N):(-1):2
    b(i) = p(i) + p_daszek * b(i + 1) + q_daszek * b(i + 2);
end

p_eval = a_0 + cos(t_eval) * b(2) + q_daszek * b(3);
% p_eval = sin(t_eval) * b(2);

end


