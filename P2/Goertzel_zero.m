function [p_zero] = Goertzel_zero(p, t_seed, eps, N_stop)
% funkcja przyjmuje wielomian trygonometryczny p i zwraca miejsce zerowe p 
% - p_zero wyznaczone za pomocą metody Newtona. Wartości p i p` w punktach 
% wyznaczane są za pomocą funkcji Goretzel_eval i Goertzel_cos_eval 
% odpowiednio. Pochodną wyznaczamy tak długo, aż spełniony zostanie jeden 
% z warunków:
% 1) różnica w kolejnych krokach wynosi mniej niż eps
% 2) wykonaliśmy więcej niż N_stop kroków
% Wówczas przerywamy algorytm i zwracamy otrzymaną wartość.
% Aby wystartować metodę Newtona funkcja przyjmuje argument t_seed.

% przed wystartowaniem metody Newtona odpowiednio przygotujmy p

N = max(length(p)) - 1;
p_falka = p;
p_falka(1) = 0;

for i = 1:1:N
    p_falka(i) = p(i) * i;
end

% przeprowadzamy procedurę Newtona:

for i = 1:1:N_stop
    p_zero = t_seed - Goertzel_eval(p, t_seed) / Goertzel_cos_eval(p_falka, t_seed);
    if abs(p_zero - t_seed) < eps
        break
    end
    t_seed = p_zero;
end

end