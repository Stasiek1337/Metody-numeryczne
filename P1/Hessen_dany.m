function [xh, detH] = Hessen_solve(H, b)
% zadaniem tej pięknej funkcji jest rozwiązanie układu równań liniowych
% postaci Hx=b za pomocą eliminacji Gaussa, gdzie H - dolna macierz
% Hessenberga przy ustalonej macierzy H i wektorze b; 
% przy okazji policzymy wyznacznik macierzy H

n = max(size(b));

H = [H b];

detH = 1;

xh = zeros(n, 1);

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

for i = 0:(n - 2)
    H(n - i - 1, n + 1) = H(n - i - 1, n + 1) - 1 / H(n - i, n - i) * H(n - i - 1, n - i) * H(n - i, n + 1);
end

for i = 1:n
    xh(i) = H(i, n + 1) / H(i, i);
    detH = detH * H(i, i);
end