# sin²ψ i eksperymentalny model Wilkensa

Wersja 20.3.0, 2026-10-10. Implementacja: `xrd_advanced_analysis.py`.

## Naprężenia sin²ψ

Program dopasowuje `d = intercept + slope sin²ψ` i oblicza
`sigma_phi = slope / (d0 * (S2/2))`. S2/2 podaje się w 1/MPa,
d oraz d0 w Å, kąty w stopniach. Wynik dotyczy jednego azymutu φ
i płaskiego stanu naprężenia; nie jest pełnym tensorem naprężeń.

Wymagane są co najmniej trzy różne wartości sin²ψ, jedna faza/hkl,
zgodny φ, niezależne d0 i rentgenowskie S2/2. S1 można podać przez API
jako metadane; nie wchodzi do nachylenia w tym modelu. Nie zastępujemy
d0 wyrazem wolnym. Nie ma automatycznej zamiany XEC na makroskopowe E/ν.
Wejściem mogą być d i u(d), albo 2θ i u(2θ) z λ; Bragg propaguje
niepewności. Wszystkie dane w GUI serii muszą używać wspólnej λ1.

Niepewności wszystkich punktów dają regresję ważoną. Brak wszystkich
u oznacza jawnie warunkową regresję nieważoną; częściowy brak blokuje
wynik. Kowariancję regresji oraz niezależne niepewności XEC/d0
przenosi się na σ i 95% CI. Wspólne systematyki λ, XEC, geometrii
i korelacje między pomiarami wymagają osobnej analizy.

Diagnostyka obejmuje residua, leverage, pomijanie całego nachylenia ψ,
liczbę niezależnych ψ, porównanie +ψ/−ψ, test zgodności liniowego modelu
z u(d) i zakres odkształceń. Znaczące splitting, krzywizna lub
niezgodne fazy/hkl/φ blokują wynik MPa. Surowy wynik próbny zachowano
oddzielnie do audytu. Status VALIDATED oznacza spełnienie jawnych
warunków tego modelu i dostarczonych deklaracji; nie certyfikuje pomiaru.

Tło metody i wymogi pomiarowe opisuje
[NPL Good Practice Guide 52](https://eprintspublications.npl.co.uk/2391/1/mgpg52.pdf).
Ograniczenia tekstury, stanu naprężenia i jakości XEC pozostają
odpowiedzialnością analizy eksperymentalnej, a nie samego testu regresji.

## Wilkens: co zaimplementowano

To ograniczony, eksperymentalny model whole-profile, **nie pełny program CMWP**.
Wykorzystuje co najmniej trzy rozdzielone refleksy jednej fazy, pełne
profile intensywności, znane centra Kα1, wektor Burgersa b, C_hkl i IPF.
Wspólnie dopasowuje ρ [m⁻²], Re [nm], medianę liczbowego lognormalnego
rozkładu średnic sferycznych domen [nm] i sigma_ln. Dla każdego refleksu
dopasowuje dodatnie pole oraz liniowe tło. Obsługuje jawny dublet Kα.

Iloczyn współczynników Fouriera rozmiaru, dyslokacji i aparatury
jest transformowany do przestrzeni odwrotnej. Stosujemy
`A_dislocation = exp(-pi/2 * rho_nm^-2 * b_nm² * g_nm^-2 * C * L_nm² * f)`.
Gęstość w m⁻² jest mnożona przez 1e-18. `g=2 sinθ/λ` w nm⁻¹,
`M=Re[m]*sqrt(rho[m^-2])`. Konwencja promienia odcięcia:
`eta=exp(7/4)*L/(2*Re)`, przez co dla małego L, `f ~ ln(Re/L)`.

Funkcja f* używa równań 3.6.35–37 z
[International Tables, rozdział 3.6](https://onlinelibrary.wiley.com/iucr/itc/Ha/ch3o6v0001/).
Całkę oblicza dylorarytm; rozwinięcie małych argumentów stosowane jest
tylko poniżej 0.001, dla stabilności numerycznej. Test porównuje wynik
z niezależną kwadraturą. Rozmiar uwzględnia uśrednioną objętościowo
autokorelację sfer po liczbowym rozkładzie lognormalnym.

Model dyslokacji wywodzi się z
[Wilkens (1970)](https://doi.org/10.1002/pssa.19700020224).
Szerszy zakres kompletnego CMWP, w tym defekty planarne, opisano w
[publikacji autorów metody (2020)](https://www.mdpi.com/2073-4352/10/7/623).
Samych FWHM lub 1/D² nie używa się do wyznaczania ρ.

## Kontrole i ograniczenia Wilkensa

Implementacja stosuje kilka startów, ograniczenia parametrów,
macierz kowariancji i korelacji parametrów logarytmicznych, CI,
kontrolę identyfikowalności, reszty każdego profilu, test liczby punktów
w szerokości oraz niezależne zagęszczenie i wydłużenie siatki Fouriera.
Niestabilny wynik, aktywne granice lub podobne dopasowania z różnymi ρ
powodują REJECTED i `rho_m2=None`. Udany wynik zawsze pozostaje CONDITIONAL.

Nie zaimplementowano pełnego modelu nakładających się faz, defektów
planarnych, mieszanych układów poślizgu, rafinacji C_hkl/centrów/IPF,
bootstrapu/profile likelihood ani automatycznej analizy czułości
na tło, IPF i zakres. CI są lokalne i warunkowe przy ustalonych tych
wejściach. Przed użyciem bezwzględnej ρ trzeba porównać wynik z niezależnym
CMWP/wzorcem i powtórzyć analizę z uzasadnionymi wariantami wejść.
Eksport whole-profile JSON jest własnym, udokumentowanym formatem wejścia,
nie deklarowanym formatem importu oficjalnego CMWP.

## Walidacja

Testy `test_advanced_analysis.py` obejmują naprężenia rozciągające i
ściskające, d0 różne od interceptu, brak danych/XEC, błędne fazy/hkl/φ,
splitting, krzywiznę, propagację u, funkcję Wilkensa względem kwadratury,
SI, granicę sfery monodyspersyjnej, FFT względem niezależnej całki
cosinusowej oraz odzyskanie parametrów whole-profile na innej siatce,
z szumem i dubletem. To walidacja numeryczna/syntetyczna.
Brak danych laboratoryjnych i zewnętrznego silnika CMWP uniemożliwia
deklarację walidacji eksperymentalnej. Bieżący wynik pełnego zestawu
testów jest zapisany w AUDIT_REPORT.md.
