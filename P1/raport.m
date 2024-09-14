% RAPORT PROJEKT_1 STANISŁAW FRELIK

% OPIS MATEMATYCZNY:
% Macierze Hessenberga w naturalny sposób pojawiają się w problemach
% znajdowania wektorów własnych oraz w obszarze matematyki zwanym teorią
% sterowania (1). Dowolna macierz jest podobna do macierzy Hessenberga a
% algorytm sprowadzający jest wykonywalny w skończenie wielu operacjach (w
% przeciwieństwie do macierzy trójkątnych). Niemniej macierze Hessenberga
% mają postać zbliżoną do macierzy trójkątnych, to znaczy dolna macierz
% Hessenberga to macierz trójkątna dolna, której superdiagonala jest
% niekoniecznie zerowa. Analogicznie definiujemy górną macierz Hessenberga.
% Nietrudno wywnioskować zatem, iż podobnie jak dla macierzy trójkątnych,
% algorytm eliminacji Gaussa będzie bardzo efektywnym sposobem
% rozwiązywania układów równań liniowych postaci Hx=b, gdzie H jest
% macierzą Hessenberga (w naszym przypadku będzie to macierz dolna). W celu
% zapewnienia wykonywalności eliminacji Gaussa będziemy rozważać wyłącznie
% dolne macierze Hessenberga, które będą diagonalnie dominujące kolumnowo.
% Nie stanowi to znaczącego ograniczenia, gdyż "ulepszona" metoda
% eliminacji Gaussa, jak na przykład GECP, poradziłaby sobie z tym
% problemem, jednakże nie jest to obiektem naszych zainteresowań w tych
% rozważaniach. Ponadto po sprowadzeniu de facto macierzy Hessenberga do
% postaci diagonalnej, co posłuży do rozwiązania układu równań, będziemy
% obliczać wyzacznik zadanej macierzy w oczywisty sposób mnożąc elementy na
% diagonali.

% OPIS IMPLEMENTACJI: 
% W celu zachowania porządku, prostoty i jednocześnie funkcjonalności
% postanowiłem napisać 3 funkcje:
% - www_gen(n) -> generuje wektor wyrazów wolnych o długości n z rozkładu
% jednostajnego z przedziału (0, 1).
% - Hessenberg_gen(n) -> generuje dolną, diagonalnie dominującą kolumnowo
% macierz Hessenberga o wymiarach nxn, o wyrazach z rozkładu jednostajnego
% z przedziału (0, 1).
% - Hessen_solve(H, b) -> rozwiązuje układ równań liniowych postaci Hx=b,
% gdzie H jest dolną macierzą Hessenberga, oraz wyznacza det(H).
% Funkcja Hessen_solve korzysta w nietrywialny sposób z postaci dolnej
% macierzy Hessenberga, to znaczy wykonywanie operacje odejmnowania od
% siebie wierszy są zwektoryzowane w taki sposób, aby nie odejmować od
% elementów w niższym wierszu elementów z wyższego wiersza, które wiemy, że
% są na pewno zerami, stąd w każdej iteracji algorytmu wykonujemy bliską
% optymalnej liczbę działań arytmetycznych. Pozostałe elementy są w 
% standardowy sposób zaimplementowane zgodnie z procedurą GE.

% SKRYPT TESTUJĄCY:
% Poniższy kod stanowi poligon testowy dla powyższych funkcji, w
% szczególności dla Hessen_solve. "Przeciwnikiem" zaimplementowanej funkcji
% w rozwiązywaniu układów równań liniowych będzie zbudowana w Matlaba
% metoda znana jako znaczek: "A\b" (choć u nas będzie się pojawiać H\b).
% Przetestujemy zarówno poprawność Hessen_solve w rozwiązywaniu zadanego
% problemu na losowo wygenerowanych macierzach (za pomocą Hessenberg_gen 
% i www_gen) jak i czas wykonywania obu algorytmów w unormowanych testach.

% Zajmijmy się najpierw poprawnością. Będą nas interesowały maksymalny i
% średni (w sensie arytmetycznym) błąd bezwzględny (oczywiście błąd średni
% otrzymamy dla każdej współrzędnej (dla n współrzędnych)) zatem policzymy
% również jego średnią).

%%
k = 1;
l = 150;
steps = 100;
blad_max_v = zeros([1 (l - k + 1)]);
blad_sredni_v = zeros([1 (l - k + 1)]);
blad_sredni_det_v = zeros([1 (l - k + 1)]);
for n = k:1:l
    blad_sredni = 0;
    blad_max = zeros([1 steps]);
    blad_sredni_det = 0;
    for i = 1:steps
        H = Hessenberg_gen(n);
        b = www_gen(n);
        [x, y] = Hessen_solve(H, b);
        r = abs(x - H\b);
        blad_sredni = blad_sredni + r;
        blad_max(i) = max(r);
        blad_sredni_det = blad_sredni_det + abs(det(H) - y);
    end
    blad_max_v(n - k + 1) = max(blad_max);
    blad_sredni_v(n - k + 1) = sum(blad_sredni / steps);
    blad_sredni_det_v(n - k + 1) = sum(blad_sredni_det / steps);
end

blad_sredni_f = sum(blad_sredni_v) / (l - k)
blad_max_f = sum(blad_max_v) / (l - k)
blad_det_f = sum(blad_sredni_det_v) / (l - k)
%%
x = k:1:l;
figure
subplot(3, 1, 1)
plot(x, nthroot(blad_sredni_v, 8) / (l - k))
title(['Błąd średni na przestrzeni testów dla n=',num2str(k),':',num2str(l)])
xlabel('Numer próby')
ylabel('Wartość błędu średniego')
subplot(3, 1, 2)
plot(x, nthroot(blad_max_v, 8) / (l - k))
title(['Błąd maksymalny na przestrzeni testów dla n=',num2str(k),':',num2str(l)])
xlabel('Numer próby')
ylabel('Wartość błędu maksymalnego')
subplot(3, 1, 3)
plot(x, nthroot(blad_sredni_det_v,3) / (l - k))
title(['Błąd średni wyznacznika na przestrzeni testów dla n=',num2str(k),':',num2str(l)])
xlabel('Numer próby')
ylabel('Wartość błędu średniego')

% WNIOSKI:
% Widzimy, iż dla rozmiaru macierzy poniżej około 60-70 bardzo dobrze (z
% dokładnością do błędu maszynowego) rozwiązuje zadany problem. Z przyczyn
% szczerze mówiąc bliżej mi nieznancych od około rozmiaru 70 różnica w
% algorytmach zaczyna drastycznie rosnąć (o rząd wielkości wraz ze
% zwiększeniem rozmiaru macierzy o 1). Ponadto Matlab zaczyna niekiedy
% zwracać błąd mówiący o tym, iż macierz jest "źle wyskalowana" lub "bliska
% macierzy singularnej". Możliwym wyjaśnieniem takiego zachowania byłoby
% to, iż kumulujące się błędy arytmetyki generowane przez algorytm
% powodują, że nie radzi on sobie z powstałymi obiektami. Co ciekawe, błąd
% wyznacznika wydaje się maleć wraz ze wzrostem rozmiaru macierzy...

%%

% W następnej części skryptu testującego postaramy się zbadać "szybkość"
% funkcji Hessen_solve w porównaniu do Matlabowego "A\b". Test szybkości 
% zostanie przeprowadzony w podobny sposób jak test błędu.

k = 1;
l = 125;
steps = 100;
czas_sredni_base_v = zeros([1 (l - k + 1)]);
czas_sredni_hess_v = zeros([1 (l - k + 1)]);
for n = k:1:l
    czas_sredni_base = 0;
    czas_sredni_hess = 0;
    for i = 1:steps
        H = Hessenberg_gen(n);
        b = www_gen(n);
        tic
        H\b;
        det(H);
        dt_base = toc;
        tic
        Hessen_solve(H, b);
        dt_hess = toc;
        czas_sredni_base = czas_sredni_base + dt_base;
        czas_sredni_hess = czas_sredni_hess + dt_hess;
    end
    czas_sredni_base_v(n - k + 1) = czas_sredni_base / steps;
    czas_sredni_hess_v(n- k + 1) = czas_sredni_hess / steps;
end

sum(czas_sredni_base_v) / (l - k)
sum(czas_sredni_hess_v) / (l - k)

x = k:1:l;
figure
subplot(3, 1, 1)
plot(x, czas_sredni_base_v)
title(['Średni czas dla Matlaba na przestrzeni testów dla n=',num2str(k),':',num2str(l)])
xlabel('Numer próby')
ylabel('Czas')
subplot(3, 1, 2)
plot(x, czas_sredni_hess_v)
title(['Średni czas dla Hess solve na przestrzeni testów dla n=',num2str(k),':',num2str(l)])
xlabel('Numer próby')
ylabel('Czas')
subplot(3, 1, 3)
plot(x, czas_sredni_hess_v - czas_sredni_base_v)
title(['Średnia różnica w czasie na przestrzeni testów dla n=',num2str(k),':',num2str(l)])
xlabel('Numer próby')
ylabel('Średnia różnica w czasie')

% WNIOSKI:
% Wygenerowane wykresy pokazują niezaskakujące prawidłowości: 
% - wbudowana metoda jest około 100 razy szybsza od zaimplementowanej (i
% tak całkiem nieźle!)
% - czas wykonywania (w obu przypadkach) rośnie wraz ze zwiększaniem
% rozmiarów macierzy 
% - różnica w czasie wykonywania również rośnie jednak dzieje się to w
% sposób względnie liniowy, nie "wybucha" on jak można było się spodziewać
% (co również uznajemy za plus)

%%
% PODSUMOWANIE:
% Po wykonaniu (jak mniemam) dość wyczerpujących (również pod względem
% zakodowania) testów możemy śmiało wnioskować, iż zaimplementowana przez
% nas metoda staniowi całkiem "niezły" sposób na rozwiązywanie układów
% równań liniowych postaci Hx=b, gdzie H jest dolną macierzą Hessenberga,
% przy pomocy metody GE oraz wyznaczania wyznacznika macierzy H.
% Nierozstrzygniętą kwestią pozostaje to co dzieje się dla jeszcze
% większych (niż n ~= 70) macierzy, których "rozwiązywanie" bez wątpliwości
% znalazłoby zastosowanie w rzeczywistych problemach. Jednak każdą pracę 
% należy w jakimś miescu zakończyć a pozostawione czytelnikowi narzędzia
% mam nadzieję, że zachęcą do dalszej zabawy z macierzami Hessenberga.
% (1) - https://www.cs.cornell.edu/~bindel/class/cs6210-f16/lec/2016-10-21.pdf

% Dziękuję za uwagę!
% Stanisław Frelik
