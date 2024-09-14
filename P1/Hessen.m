function [xh, detH] = Hessen(n)
% zadaniem tej pięknej funkcji jest rozwiązanie układu równań liniowych
% postaci Hx=b za pomocą eliminacji Gaussa, gdzie H - dolna macierz
% Hessenberga
% przy okazji policzymy wyznacznik macierzy H

% generujemy losową macierz Hessenberga i wektor wyrazów wolnych
% (poprawność "n" sprawdzą wewnętrzne funkcje")

H = Hessenberg_gen(n);

b = www_gen(n);

H = [H b];

detH = 1;

xh = zeros(n, 1);

% kluczową obserwacją, która pozwoli nam oszczędzić wielu operacji będzie
% to iż dolna macierz Hessenberga jest bardzo zbliżona do macierzy
% dolno-trójkątnej

% klasyczny Gaussik

for i = 1:n
    if H(i, i) == 0
        continue
    end
    for j = i:(n - 1)
        if H(j + 1, i) == 0
            continue
        end
        h = 1/ H(i, i) * H(j + 1, i);
        H(j + 1,1:(i + 1)) = H(j + 1,1:(i + 1)) - h * H(i,1:(i + 1));
        H(j + 1, n + 1) = H(j + 1, n + 1) - h * H(i, n + 1);
    end
end

% redukujemy superprzekątną i wektor wyrazów wolnych; część operacji jest
% wykonywana w celu sprowadzenia macierzy do wizualnie ładnej postaci,
% tutaj już operujemy na pojedyńczych miejscach a nie na wierszach więc
% możnaby to pominąć

for i = 0:(n - 2)
    h = 1 / H(n - i, n - i) * H(n - i - 1, n - i);
    H(n - i - 1, n - i) = 0;
    H(n - i - 1, n + 1) = H(n - i - 1, n + 1) - h * H(n - i, n + 1);
end

% spisujemy wyniki i liczymy wyznacznik

for i = 1:n
    xh(i) = H(i, n + 1) / H(i, i);
    detH = detH * H(i, i);
end
