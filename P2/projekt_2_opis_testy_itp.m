% Raport do projektu nr 2 
% Laboratorium z Metod Numerycznych
% Stanisław Frelik, 21.05.2023
%
%    Paragraf 1: Wstęp.
%
%                           --- Problem ---
% Należy wykorzystać metodę Newtona do wyznaczania zer wielomianu
% trygonometrycznego, przy czym do obliczania wartości tego wielomianu oraz
% jego pochodnych należy zastosować metodę Goertzela.
% Przez wielomian trygonometryczny będziemy rozumieli wielomian postaci:
% p(t) = SUMA(od n=0 do N) a_n sin(nt) dla N naturalnych i będziemy go 
% oznaczać przez p. 
%                      --- Jak go rozwiązać? ---
% Do problemu podeszliśmy modułowo, to znaczy najpierw zajęliśmy się
% napisaniem funkcji do obliczania wartości p, następnie do obliczania p' a
% na koniec połączyliśmy dwie poprzednie wraz z metodą Newtona. Użyliśmy
% implementacji dobrze nam znanych z wykładu.
%                     --- O algorytmie Goertzela ---
% Pierwsze wpisanie frazy algorytm Goertzela w "Googla" nie przyniosło
% znaczących rezultatów. Trochę niepokojącym wydały się przejawiające się w
% wielu nagłówkach frazy takie jak: dyskretna transformata Fouriera,
% analiza sygnałów czy analiza pojedyńczych częstotliwości. Lektura
% pokazała, iż faktycznie jest to pewna modyfikacja do znanej każdemu
% inżynierowi metody FFT (Fast Fourier Transform) [3] poprawiająca jej
% "osiągi" w pewnych sytuacjach. Tak również została zaimplementowana w 
% Matlabie metoda 'goertzel()'. De facto ta informacja nie powinna dziwić,
% w końcu mamy go zastosować do wielomianu trygonometrycznego p obciętego
% do samych sinusów - bliskiego krewnego szeregów Fouriera. Wówczas jednak
% poszukiwania informacji o algorytmie Goertzela przeszły na materiały z
% wykładu z Metod Numerycznych. 
% 
% Algorytm Goertzela pozwala nam w swojej ogólności na ewaluowanie
% skończonych szeregów potęgowych z ciała liczb zespolonych. Algorytm
% opiera się na kroku 'podzielenia' wejściowej sumy przez odpowiedni
% trójmian, wyznaczenie 'od końca' współczynników nowego szeregu i spisaniu
% szukanej odpowiedzi. Jest on bardzo zbliżony do algorytmu Hornera, jednak
% dla pewnych szczególnych wielomianów (np. o rzeczywistych
% współczynnikach i ewaluacji w punkcie takim, że |z|=1) złożoność
% algorytmu Goertzela wynosić może nawet N (w porównaniu do stałego 4N dla
% algorytmu Hornera). Co więcej algorytm Goertzela potrafi ewaluować
% interesujące nas wielomiany trygonometryczne składające się zarówno z
% sinusów jak i cosinusów. Te drugie wykorzystamy później w dość "sprytny"
% sposób do ewaluacji p'.
%
%   Paragraf 2: Technikalia.
%
%                       --- Pewne wyjaśnienia ---
% Bedą nas interesować wielomiany trygonometryczne postaci:  
% p(t) = SUMA(od n=0 do N) a_n sin(nt) dla N naturalnych. Statystyką
% dostateczną takiej sumy jest zdecydowanie wektor: p = [a0, a1, ..., aN],
% gdzie an to pewne stałe zespolone. 
%                       --- O liczeniu pochodnej ---
% Algorytm Goertzela pozwala nam na ewaluację zarówno wielomianów
% trygonometrycznych z sinusami jak i cosinusami. Nietrudno policzyć na
% kartce, że jeżeli p(t) = SUMA(od n=0 do N) a_n sin(nt) to 
% p'(t) = SUMA(n=1 do N) n a_n cos(nt). Ta jakże niepodniosła wiadomość
% będzie stanowić podstawę naszych rozważać.
%                           --- O zerach ---
% Łatwo pokazać, że wielomian p jest funkcją okresową o okresie 2pi. Trochę
% trudniej pokazać, że wielomian p posiada dokładnie 2N zer na przedziale
% [0, 2pi). Ta informacja oznacza, że zera p będą dość mocno zagęszczone
% dla dużych N. Ponadto biorąc pod uwagę fakt sumowania sinusów, w 
% okolicach bliskich początku lub końcu okresu będziemy doświadczać bardzo 
% duże skoki w wartościach p. Co za tym idzie metoda Newtona może nie
% stanowić najlepszego narzędzia do znajdowania miejsc zerowych p co
% przetestujemy w dalszej części.
%                       --- O współczynnikach ---
% W implementacji algorytmu z wykładu dostajemy wolną rękę co do a_n - są
% to pewne zespolone skalary. Postaramy się sprawdzić dla jak bardzo
% skomplikowanych a_n jesteśmy w stanie skutecznie zaimplementować
% rozwiązania problemu oraz jak ich wybór wpływa na 'osiągi'.
%
% Na początku skupimy się na zdegenerowanej wersji problemu, gdzie dla
% każdego n = 1:N a_n = 1. Przejdzmy zatem do przypadku "podstawowego".
%
%   Paragraf 3: Miejsca zerowe w przypadku podstawowym.
%
%                           --- Redukcje ---
% Nietrudno zauważyć, że t = 0 jest miejscem zerowym wielomianu p oraz, że 
% dla n = 0 dodajemy wyłącznie funkcję stale równą 0 co nie wpływa na sumę.
% Można dowieść, iż wielomian p posiada 2N miejsc zerowych na przedziale 
% [0, 2pi). Ponadto dla dużych N zera p będą coraz bardziej gęste co 
% ilustruje poniższy obrazek:
%%
x = 0:1e-3:2*pi;
figure
subplot(3, 1, 1)
N = 3;
y = zeros(1, max(length(x)));
for i = 1:max(length(x))
    for j = 0:N
        y(i) = y(i) + sin(x(i) * j);
    end
end
plot(x, y)
grid on
title(['Wykres przykładu podstawowego dla N=',num2str(N)])
subplot(3, 1, 2)
N = 10;
y = zeros(1, max(length(x)));
for i = 1:max(length(x))
    for j = 0:N
        y(i) = y(i) + sin(x(i) * j);
    end
end
plot(x, y)
grid on
title(['Wykres przykładu podstawowego dla N=',num2str(N)])
subplot(3, 1, 3)
N = 50;
y = zeros(1, max(length(x)));
for i = 1:max(length(x))
    for j = 0:N
        y(i) = y(i) + sin(x(i) * j);
    end
end
plot(x, y)
grid on
title(['Wykres przykładu podstawowego dla N=',num2str(N)])

%%
% Powoduje to, iż metoda Newtona może mieć znaczący problem w znajdywaniu
% niektórych zer. Poniżej sprawdzimy czy faktycznie tak jest i do jakiego
% stopnia jest to problem.
%                       --- Znaleść je wszystkie ---
% Spróbujmy "rzucić" Matlabowi wektor różnych punktów początkowych w celu
% sprawdzenia ile zer nam "wypluje". Zacznijmy od przypadku gdy N=3
% (na przedziale [0, 2pi) mamy 6 zer w punktach 0, pi/2, 2pi/3, pi, 4pi/3 
% i 3pi/2).
%%
% Czy wogóle znajdziemy jakiekolwiek zero?
p = [1 1 1 1]; % przykład podstawowy dla N=3
t_seed = 1.34; eps = 1e-20; N_stop = 1e7;
disp(Goertzel_zero(p, t_seed, eps, N_stop))
% Istotnie działa! :)
%%
% Wybierzmy losowe 6 seedów z przedziału (0, 2pi) z rozkładu jednostajnego:
rng(1337)
p_zero = zeros(1, 6);
tv = rand(1, 6) * 2 * pi;
p = [1 1 1 1]; eps = 1e-10; N_stop = 1e7;
for i = 1:6
    p_zero(i) = Goertzel_zero(p, tv(i), eps, N_stop);
end
disp(p_zero)

% Jak widzimy dostaliśmy coś czego się (niestety) spodziewaliśmy, przy czym
% pierwiastek pi pojawił się aż 4 razy a raz dolecieliśmy aż do 9pi/2.
% Co z tym faktem zrobimy? Spróbujemy najpierw "na siłę"...
%%
rng(420)
p = [1 1 1 1]; eps = 1e-20; N_stop = 1e8;
p_zero = [];
for i = 1:1000000
    t_seed = rand * 2 * pi;
    p_cand = round(Goertzel_zero(p, t_seed, eps, N_stop), 3);
    if ~ismember(p_cand, p_zero)
        if p_cand >= 0 && p_cand < 2 * pi
            p_zero = horzcat(p_zero, p_cand);
        end
    end
end
disp(sort(p_zero))
% no niestety jest to prakycznie niewykonalne,  nietrafiamy zupełnie w dwa
% pierwiastki pomimo znaczącej liczby prób (co ciekawe pomimo potężnej
% pętli cała sekcja wykonuje się zaskakująco szybko!).
% Spróbujmy teraz dla większych N, np N = 15 (czyli 30 pierwiastków na
% przedziale [0, 2pi)) i sprawdzmy powyższy kod:
%%
rng(420)
p = [1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1]; eps = 1e-10; N_stop = 1e8;
p_zero = [];
for i = 1:100000
    t_seed = rand * 2 * pi;
    p_cand = round(Goertzel_zero(p, t_seed, eps, N_stop), 3);
    if ~ismember(p_cand, p_zero)
        if p_cand >= 0 && p_cand < 2 * pi
            p_zero = horzcat(p_zero, p_cand);
        end
    end
end
disp(sort(p_zero))
disp(max(length(p_zero)))

% Jest niewiele lepiej. Następnie sprawdzimy pewną hipotezę badawczą: czy
% pewne pierwiastki p mają większą szansę na bycie odnalezionym? Pytamy się
% czy niektóre wartości "przyciągają" metodę Newtona co powodowałoby
% trudności z odnajdywaniem ich znaczącej liczby? 
%%
p = [1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1]; eps = 1e-30; N_stop = 1e6;
p_zero = [];
for i = 1:1000
    t_seed = rand * 2 * pi;
    p_cand = round(Goertzel_zero(p, t_seed, eps, N_stop), 3);
    if p_cand >= 0 && p_cand < 2 * pi
        p_zero = horzcat(p_zero, p_cand);
        disp(i)
    end

end
[p_unique, ~, index] = unique(p_zero);
p_count = histcounts(index, length(p_unique));
figure
subplot(2,1,1);
format('shortG')
disp([p_unique; p_count]);
bar(p_unique, p_count)
title('Ile razy trafiliśmy zera')
xlabel('Zera p')
ylabel('Liczba trafień')
subplot(2,1,2);
x = 0:1e-3:2*pi;
N = 15;
y = zeros(1, max(length(x)));
for i = 1:max(length(x))
    for j = 0:N
        y(i) = y(i) + sin(x(i) * j);
    end
end
plot(x, y)
grid on
title(['Wykres przykładu podstawowego dla N=',num2str(N)])

% Na wygenerowany obrazku obserwujemy bardzo ciekawą sytuację, okazuje się
% że dla N = 15 pierwiastki w okolicach pi (np 3pi/4, 14pi/15, pi, 16pi/15,
% 9pi/8) są praktycznie nietrafiane w porównaniu do ich bardziej skrajnych
% kolegów. Ponadto im bardziej skrajny pierwiastek tym mniej razy został
% trafiony. Wyjątkowo niespodziewanie wygląda spadek liczby trafień w
% okolicach pi/2 i 3pi/2 a sam wykres wydaje się być symetryczny! 
%
% Powodów takiego zachowania może być wiele. Patrząc na zamieszczony 
% poniżej wykres rozważanej funkcji możemy się domyślać, że dla starterów
% metody Newtona bliskich pi pierwsza styczna będzie stosunkowo pozioma co
% sprawi, że zostaniemy "wystrzeleni" na pojawiające się przy krańcach
% przedziału [0, 2pi) górki, które następnie sprowadzą nas na okoliczne
% miejsca zerowe. Możemy zatem twierdzić, iż pierwiastki "środkowe" mają
% własność odpychającą, natomiast skrajne - "przyciągającą".
%%
%   Paragraf 4: Zabawa z a_n.
%
%                   --- O pewnych limitach ---
% W ogólności chcielibyśmy uniknąć przypadku, że dla pewnego n, a_n będzie
% niezdefiniowane, np. a_n = 1/n. Niemniej obrazek właśnie dla a_n = -1/n
% (przy czym będziemy sumować od n = 1) wygląda interesująco...
%                   --- O szeregach Fouriera ---
figure
x = 0:1e-3:8*pi;
subplot(2,1,1);
N=2;
y = zeros(1, max(length(x)));
for i = 1:max(length(x))
    for j = 1:N
        y(i) = y(i) - sin(x(i) * j) / j;
    end
end
plot(x, y)
grid on
title(['Wykres dla a_n = -1/n dla N=',num2str(N)])
subplot(2,1,2);
N=25;
y = zeros(1, max(length(x)));
for i = 1:max(length(x))
    for j = 1:N
        y(i) = y(i) - sin(x(i) * j) / j;
    end
end
plot(x, y)
grid on
title(['Wykres dla a_n = -1/n dla N=',num2str(N)])

% Przyglądając się w szczególności wykresowi dla N = 25 dostrzegamy na nim
% zapewne znajomy niektórym obrazek - częściowy szereg Fouriera dla funkcji
% mantysy. Istotnie:
%%
figure 
x = 0:1e-3:8*pi;
subplot(3, 1, 1)
N = 1;
y = zeros(1, max(length(x)));
for i = 1:max(length(x))
    for j = 1:N
        y(i) = y(i) - sin(x(i) * j) / j;
    end
end
plot(x, y)
hold on
x_mod = x / (2 * pi);
z = pi*(x_mod-floor(x_mod)-1/2);
plot(x, z)
grid on
title(['Częściowy szereg Fouriera dla mantysy dla N=',num2str(N)])
hold off
subplot(3, 1, 2)
N = 7;
y = zeros(1, max(length(x)));
for i = 1:max(length(x))
    for j = 1:N
        y(i) = y(i) - sin(x(i) * j) / j;
    end
end
plot(x, y)
hold on
x_mod = x / (2 * pi);
z = pi*(x_mod-floor(x_mod)-1/2);
plot(x, z)
grid on
title(['Częściowy szereg Fouriera dla mantysy dla N=',num2str(N)])
hold off
subplot(3, 1, 3)
N = 50;
y = zeros(1, max(length(x)));
for i = 1:max(length(x))
    for j = 1:N
        y(i) = y(i) - sin(x(i) * j) / j;
    end
end
plot(x, y)
hold on
x_mod = x / (2 * pi);
z = pi*(x_mod-floor(x_mod)-1/2);
plot(x, z)
grid on
title(['Częściowy szereg Fouriera dla mantysy dla N=',num2str(N)])
hold off

%%
%                   --- Zera dla różnych a_n ---
% Powróćmy do naszego zadania szukania miejsc zerowych p. Rozważymy kila
% przykładów dla różnych a_n i różnych wartości N. Tam gdzie to będzie
% konieczne będziemy sumować od n = 1.
%
% Przykład 1) Niech a_n = sqrt(abs(sin(n))):
% Przełamując delikatnie czwartą ścianę chciałbym powiedzieć, że bardzo
% lubie ładne i ciekawe wykresy, stąd będą one wciąż obecne, nie chciałbym
% będąc na miejscu szanownego czytelnika patrzeć tylko na "numerki" czyli
% wyliczone miejsca zerowe. 
figure 
x = 0:1e-3:4*pi;
subplot(2,1,1);
N=3;
y = zeros(1, max(length(x)));
for i = 1:max(length(x))
    for j = 0:N
        y(i) = y(i) + sin(x(i) * j) * sqrt(abs(sin(j)));
    end
end
plot(x, y)
grid on
title(['Wykres dla a_n = a_n = sqrt(abs(sin(n))) dla N=',num2str(N)])
subplot(2,1,2);
N=31;
y = zeros(1, max(length(x)));
for i = 1:max(length(x))
    for j = 0:N
        y(i) = y(i) + sin(x(i) * j) * sqrt(abs(sin(j)));
    end
end
plot(x, y)
ylim([-10 10])
grid on
title(['Wykres dla a_n = a_n = sqrt(abs(sin(n))) dla N=',num2str(N)])

% Dostaliśmy zupełny chaos ale o to nam chodziło! Co ciekawe powyższe
% przykłady różnią się w zasadniczy sposób, zasadniczy dla metody Newtona -
% to znaczy dla N = 3 obserwujemy bardzo małe wahania wartości pochodnej
% funkcji - jest ona tam dość płaska co może skutkować wystrzeleniem ciągu
% przybliżeń pomimo bliskści przybliżenia początkowego do pewnego zera.
% Natomiast odwrotnej sytuacji spodziewamy się dla N = 50.
%%
rng(420)
p = [0 sqrt(abs(sin(1))) sqrt(abs(sin(2))) sqrt(abs(sin(3)))]; eps = 1e-10;
N_stop = 1e8; p_zero = [];
for i = 1:1000
    t_seed = rand * 2 * pi;
    p_cand = round(Goertzel_zero(p, t_seed, eps, N_stop), 3);
    if p_cand >= -100 && p_cand < 100
        p_zero = horzcat(p_zero, p_cand);
    end
end
[p_unique, ~, index] = unique(p_zero);
p_count = histcounts(index, length(p_unique));
figure
subplot(2,1,1);
format('shortG')
%disp([p_unique; p_count]);
p_count = log(p_count + 1);
bar(p_unique, p_count)
N = 3;
title(['Ile razy trafiliśmy zera dla N = ',num2str(N)])
xlabel('Zera p')
ylabel('Liczba trafień (zlogarytmowana)')
N = 1:50;
p = sqrt(abs(sin(N))); eps = 1e-7;
N_stop = 1e3; p_zero = [];
for i = 1:1e5
    t_seed = rand * 2 * pi;
    p_cand = round(Goertzel_zero(p, t_seed, eps, N_stop), 3);
    if p_cand >= -100 && p_cand < 100
        p_zero = horzcat(p_zero, p_cand);
    end
end
[p_unique, ~, index] = unique(p_zero);
p_count = histcounts(index, length(p_unique));
subplot(2,1,2);
format('shortG')
%disp([p_unique; p_count]);
disp(min(p_unique))
disp(max(p_unique))
p_count = log(p_count + 1);
bar(p_unique, p_count, 'BarWidth', 150)
N = 50;
title(['Ile razy trafiliśmy zera dla N = ',num2str(N)])
xlabel('Zera p')
ylabel('Liczba trafień (zlogarytmowana)')

% Okazuje się, że metoda Newtona nie sprawuje się najlepiej w tym
% wymagającym przypadku. Ciekawy kształt otrzymanych barplotów zapewne
% wprawiłby niejednego probabilistyka w euforię, obserwując coś
% przypominającego gęstość rozkładu normalnego (co prawda z dość ciężkimi
% ogonami). Jednakże losowaliśmy początkowe przybliżenia z przedziału (0,
% 2pi) więc normalnym jest to, że w tych okolicach trafimy najwięcej zer.
% Niemniej jest co coś nad czym warto się na chwilę zatrzymać...
% Spełniając kronikarski obowiązek rozważymy jeszcze jeden ciąg
% współczynników a_n = sin(n^n)cos(1/n).
%%
p = zeros(1, 75);
for i = 1:75
    p(i) = sin(i^i)*cos(1/i);
end
figure 
x = 0:1e-3:4*pi;
N=75;
y = zeros(1, max(length(x)));
for i = 1:max(length(x))
    for j = 1:N
        y(i) = y(i) + sin(x(i) * j) * p(j);
    end
end
plot(x, y)
grid on
title(['Wykres dla a_n = sin(n^n)cos(1/n) dla N=',num2str(N)])

% Wyznaczenie analitycznych zer powyższej funkcji poza trywialnymi
% przypadkami byłoby zdecydowanie wymagające. Niemniej przychodzą nam z
% pomocą metody numeryczne:
%%
rng(420)
p = zeros(1, 75);
for i = 1:75
    p(i) = sin(i^i)*cos(1/i);
end
eps = 1e-10; N_stop = 1e5; p_zero = [];
for i = 1:10000
    t_seed = rand * 2 * pi;
    p_cand = round(Goertzel_zero(p, t_seed, eps, N_stop), 3);
    if p_cand >= -100 && p_cand < 100
        p_zero = horzcat(p_zero, p_cand);
    end
end
[p_unique, ~, index] = unique(p_zero);
p_count = histcounts(index, length(p_unique));
figure
format('shortG')
%disp([p_unique; p_count]);
p_count = (p_count + 1);
bar(p_unique, p_count, 'BarWidth', 150)
N = 75;
title(['Ile razy trafiliśmy zera dla N = ',num2str(N)])
xlabel('Zera p')
ylabel('Liczba trafień (zlogarytmowana)')

% Otrzymaliśmy dość ciekawy rozkład zer, powiedziałbym że przypomina wieżę
% Barad-dûr z powieści Tolkiena.
%%
%   Paragraf 5: Porównania.
%
%                    --- Matlab vs Newton ---
% Policzyliśmy już dużo zer ale czy te wszystkie obliczenia miały jakiś
% sens? Odpowiedź brzmi tak i nie. Testowanie "na oko" poprawności
% wyznaczania zer tak skomplikowanych funkcji jak p jest daremne. Ponadto
% wyznaczenie wszystkich zer w analityczny sposób jest możliwie wykonalne 
% wyłącznie dla przypadku podstawowego. Dlatego w tym paragrafie postaramy
% się porównać wyznaczone przez nas zera z Matlabową funkcją fzero - opartą
% na metodzie bisekcji wspartej odwrotną interpolacją kwadratową (idea za
% funkcją fzero sprawia, iż porównanie jej do metody Newtona powinno mieć
% sens). Pewnym problemem jest to, że fzero jest w stanie znaleść tylko
% jedno zero funkcji w danym wywołaniu. Stąd przeprowadzimy w pętli pewną
% liczbę testów polegających na porównaniu wyników znajdywania zera dla
% tych samych argumentów początkowych i sprawdzimy czy i w jak istotnym
% stopniu otrzymane wyniki różnią się .
% 1) Niech p(x) = sin(x) + sin(2x) + sin(3x) + sin(4x) + sin(5x + sin(6x)

Q = 5000; eps = 1e-50; N_stop = 1e6; h = (2 * pi) / Q; p = [1 1 1 1 1 1 1];
newton_zeros = zeros(1, Q); darboux_zeros = zeros(1, Q); t_seed = h;
for i = 1:Q
    newton_zeros(i) = Goertzel_zero(p, t_seed, eps, N_stop);
    darboux_zeros(i) = fzero(@(x) sin(x)+sin(2*x)+sin(3*x)+sin(4*x)+sin(5*x)+sin(6*x), t_seed);
    t_seed = t_seed + h;
    disp(i)
end
D = abs(newton_zeros - darboux_zeros);
Dd = ['max: ', num2str(max(D)), ', mean: ', num2str(mean(D)), ', std: ', num2str(std(D))];
disp(Dd)
x = h:h:(2 * pi);
figure
subplot(2, 1, 1);
plot(x, newton_zeros)
title("Znalezione zera a punkt startowy - metoda Newtona");
subplot(2,1,2);
plot(x, darboux_zeros)
title("Znalezione zera a punkt startowy - metoda fzero");

%%
% Jak widzimy... dużo nie widzimy... Niestety okazuje się, że
% zaimplementowana przez nas metoda nie radzi sobie zbyt dobrze w tym
% przypadku. Spróbujmy na podstawie otrzymanych wyników z fzero lekko
% zmodyfikować wektor newton_zeros w celu otrzymania bardziej czytelnego
% obrazu sytuacji:

newton_zeros_mod = newton_zeros;
for i = 1:Q
    if abs(newton_zeros_mod(i) - darboux_zeros(i)) > max(darboux_zeros)
        newton_zeros_mod(i) = 0;
    end
end
for i = 1:Q
    if abs(newton_zeros(i) - darboux_zeros(i)) > 1e-10
        newton_zeros(i) = 0;
    end
end
figure
subplot(2, 1, 1);
plot(x, newton_zeros_mod, x, darboux_zeros)
title("Wykes otrzymanych zer zmodyfikowany dla dużych błędów")
subplot(2, 1, 2);
plot(x, newton_zeros, x, darboux_zeros)
title("Wykes otrzymanych zer zmodyfikowany dla małych błędów")

% Teraz dobrze widać, że faktycznie zaimplementowana przez nas metoda (poza
% pewnymi obszarami, na których "daje radę" - czyli na tych gdzie "schodki"
% się pokrywają) nie sprawdza się najlepiej. Powodów takiego zachowania
% może być wiele. Najprościej próbować temu zaradzić można przez
% zwiększenie dokładności metody Newtona, to znaczy zwiększenie parametru
% N_stop i zmiejszenie parametru eps. Niestety biorąc pod uwagę liczne
% próby zabawy z ich wartościami okazuje się, że implementacja metody
% Newtona jest na tyle niestablina, iż nawet takie zabiegi nie przynosiły
% znaczącego efektu. Ponadto czas wykonywania wzrastał wówczas znacząco, co
% kwestionowało całość podjętych wysiłków. Serdecznie zachęcamy czytelnika
% do podjęcia własnych prób i testów. Niewykluczonym jest, że dla wyjątkowo
% skomplikowanych przypadków zaimplementowana przez nas metoda okazałaby
% się mniej gorsza niż w pokazanym wyżej przykładzie. 
%%
%   Paragraf 6: Podsumowanie, myśli końcowe.
%
%                        --- Podsumowanie --- 
% Jak się można było spodziewać, szukanie zer wielomianów
% trygonometrycznych stanowi duże wyzwanie, do którego zastosowanie
% metod numerycznych nie zawsze okazuje się proste, przyjemne i przede
% wszystkim skuteczne. Oczywiście powinniśmy docenić to co udało się nam
% osiągnąć a zapewne drobne poprawki, ulepszenia czy nowe pomysły
% sprawiłyby, że nasza implementacja rozwiązania problemu stałaby się
% jeszcze lepsza.
%                       --- Myśli końcowe ---
% Chciałbym podzielić się w dosłownie jednym zdaniu swoim wrażeniem na
% temat wykonanej powyżej pracy. Pomimo tego, że momentami nie stanowi ona
% perfekcyjnie merytorycznej i "na temat" roboty to stanowiła prawdziwą
% przyjemność w jej spełnianiu, dostarczyła wielu zagwozdek, w pozytywny
% sposób napędziła szare komórki pod kątem matematyczno/informatycznym jak 
% i pozwoliła na wolny potok myśli (czasami luźno) wziązanych z problemem.
% 
% Dziękuję bardzo za uwagę, mam nadzieję, że macie Państwo jakieś pytania
% ;)
% 
%   Paragraf 7: Bibliografia.
%
% [1] - wykład 
% [2] - algorytm Goertzela, też wykład
% [3] - https://www.mstarlabs.com/dsp/goertzel/goertzel.html