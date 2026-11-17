# Analiza profili linii XRD

Data: 2026-10-09. Moduł: `xrd_profile_analysis.py`.

Moduł wykonuje lokalne dopasowanie jednego refleksu z dubletem, korekcję
instrumentalną, warunkowego Scherrera, Williamson–Halla i dopasowanie parametru
sieci kubicznej. Nie zastępuje walidacji na danych wzorca i próbek rzeczywistych.
Obliczenia nie korzystają z wizualnego wygładzania ani z szerokości
half-prominence z detektora pików. Surowe tablice wejściowe nie są modyfikowane.

## Jednostki i model profilu

- `two_theta*`, `beta_deg`, `*fwhm_deg`: **stopnie osi 2θ**.
- `wavelength_A`, `d_A`, `a_A`: Å.
- `D_nm`: nm; wewnątrz Scherrera i WH następuje jawna konwersja λ Å → nm.
- `epsilon`: bezwymiarowe pozorne mikroodkształcenie; procent to `100*epsilon`.
- `amplitude`: pole składowej Kα1, w jednostkach intensywność × stopień.
- `eta`: udział pola Lorentza w pseudo-Voigcie, nie stosunek wysokości.

Pseudo-Voigt jest sumą unormowanych Gaussa i Lorentza o tej samej FWHM `H`.
`G` i `L` oznaczają szerokości składowych splotu przybliżanego przez TCH:

```text
H^5 = G^5 + 2.69269 G^4 L + 2.42843 G^3 L^2
            + 4.47163 G^2 L^3 + 0.07842 G L^4 + L^5
q = L/H
eta = 1.36603 q - 0.47719 q^2 + 0.11116 q^3
```

To przybliżenie kształtu Voigta opisane przez
[Thompson, Cox i Hastings (1987)](https://doi.org/10.1107/S0021889887087090).
Współczynniki podano również w
[dokumentacji technicznej TOPAS](https://xraysrv.wustl.edu/web/xrd/bruker/man/TOPAS%205%20Technical%20Reference.pdf).
Test niezależnie porównuje otrzymaną szerokość z dokładnym splotem
`scipy.special.voigt_profile` dla kilku mieszanin.

Składowe dubletu mają wspólne `d`, `G`, `L` i ustalony stosunek pól.
Pozycja Kα2 wynika z Bragga, nigdy z niezależnego przesunięcia:

```text
d = lambda1 / (2 sin(theta1))
2theta2 = 2 asin((lambda2/lambda1) sin(theta1))
y = B0 + B1 (x-xref) + A PV1(x) + ratio A PV2(x)
```

Wartości początkowe Cu z raportów: λ1 = 1.540598 Å, λ2 = 1.544426 Å,
Kα2/Kα1 = 0.5. Są parametrami wejściowymi, a nie kalibracją aparatu.
`kalpha2_ratio=0` włącza model pojedynczej linii. Model wspólnej szerokości
jest lokalnym przybliżeniem dla bliskich składowych dubletu.

`fit_peak_profile` wymaga wybranego zakresu obejmującego jeden refleks/dublet
i otoczenie tła, minimum 12 punktów, rosnącej osi i skończonych wartości.
Dopasowuje sześć parametrów: pozycję, pole, G, L, stałe i liniowe tło.
Ograniczenia wymuszają dodatnie szerokości/pole, geometrię Bragga i pozycję
w zakresie. Stosowane są dwa starty kształtu oraz ważone/nieważone
[ograniczone najmniejsze kwadraty SciPy](https://docs.scipy.org/doc/scipy/reference/generated/scipy.optimize.least_squares.html).

Wynik `fwhm_obs_deg` dotyczy jednej składowej Kα1 modelu. Nie jest szerokością
całego nierozdzielonego dubletu ani half-prominence width.

## Niepewności i ograniczenia

Z podanym `sigma_y` stosowane są bezwzględne niepewności standardowe
intensywności i macierz `(JᵀJ)^-1`. Bez `sigma_y` macierz skaluje wariancja
residuów, a χ²red jest niedostępne. CI profilu/sieci/WH to przybliżenie
normalne przy znanych bezwzględnych σ, a przy σ oszacowanym z residuów —
przedział t-Studenta z właściwą liczbą stopni swobody. To lokalna liniowa
aproksymacja; nie jest bootstrapem ani profile likelihood.

Sprawdzany jest rząd i uwarunkowanie Jacobiego po skalowaniu kolumn.
Słaba identyfikowalność, parametry na ograniczeniach, wąski zakres i zbyt
mała liczba punktów na szerokość generują ostrzeżenia. Wynik obejmuje
residua, RMSE, R², χ²red (jeśli ma sens), pełną macierz kowariancji i CI.
Samo wysokie R² nie nadaje statusu `VALIDATED`.

`bragg_d` propaguje niezależne `u(2θ)` oraz `u(λ)`, przy czym pochodna
względem pozycji używa **radianów 2θ**. `fit_peak_profile(...,
u_wavelength_A=...)` przekazuje tę niepewność do `d_A`. Niepewność λ nie
jest dodatkowym dopasowywanym parametrem dubletu.

Niepewności nie obejmują nieznanej asymetrii, nakładania faz, błędów
kalibracji osi, tekstury, błędu wybranego zakresu ani nieznanych korelacji
między pomiarami. Parametry i CI przy aktywnym ograniczeniu wymagają
osobnej kontroli; symetryczne CI mogą wtedy nie być poprawne.

## Korekcja instrumentalna

`InstrumentBroadeningModel` przechowuje punkty `(2θ, G, L)`, opcjonalne
`u(G),u(L)` oraz metadane wzorca/aparatury. Wymaga co najmniej dwóch
rosnących kątów. Pomiędzy punktami stosowana jest liniowa interpolacja;
**ekstrapolacja jest zabroniona**. Niepewności niezależnych węzłów są
łączone kwadratowo z wagami interpolacji.

```text
G_sample = sqrt(G_observed^2 - G_instrument^2)
L_sample = L_observed - L_instrument
H_sample = TCH(G_sample, L_sample)
```

Nie odejmuje się całkowitych mieszanych FWHM. Ujemne składowe po
dekonwolucji odrzucają wynik, zamiast dawać sztucznie zerowy wkład.
Brak pokrycia kąta przez wzorzec daje `INSUFFICIENT_DATA`.
Niezgodna λ daje `REJECTED`. Zgodność optyki i geometrii wymaga
potwierdzenia przez użytkownika na podstawie metadanych pomiaru.

[NIST SRM 660c](https://www.nist.gov/publications/certification-standard-reference-material-660c-powder-diffraction)
służy kalibracji pozycji i kształtu linii/IPF. Program nie zawiera
uniwersalnych szerokości dla SRM 660c: trzeba zmierzyć wzorzec we własnej
konfiguracji i uwzględnić jego ewentualny wkład w poszerzenie.

Korekcja propaguje kowariancję G/L profilu oraz dostępne niepewności
instrumentu. Gdy ich brak, `u_fwhm_corrected_deg` (jeśli dostępne) obejmuje
tylko część budżetu; `uncertainty_complete=False` i ostrzeżenie są jawne.
`uncertainty_complete=True` oznacza wyłącznie obecność tych losowych
składowych w zadanym modelu; nie obejmuje systematyk. Przy `G_sample`
bliskim zeru liniowa propagacja jest osobliwa, więc u nie jest raportowane.

## Scherrer i Williamson–Hall

Scherrer używa skorygowanej instrumentalnie FWHM w radianach:
`D_nm=K*lambda_A*0.1/(beta_rad*cos(theta))`. `K=0.9` jest domyślnym
założeniem. Można przekazać u(K), u(λ), u(β) i u(2θ). Wynik zawsze
opisuje warunkową długość domen koherentnego rozpraszania.

K zależy m.in. od definicji szerokości, kształtu i rozkładu rozmiarów;
zależności te omawia
[publikacja IUCr o interpretacji równania Scherrera](https://journals.iucr.org/j/issues/2024/05/00/oc5037/index.html).
Poszerzenie odkształceniowe nie jest automatycznie oddzielane przez
jednorefleksowe równanie Scherrera. Jeśli 95% CI β obejmuje zero,
`D_nm=None` zamiast nieograniczonej pozornej precyzji.

WH stosuje `x=4 sin(theta)`, `y=beta_rad cos(theta)`; nachylenie jest ε,
a wyraz wolny `K lambda / D`. To przybliżony izotropowy model FWHM-WH
zainspirowany klasyczną analizą
[Williamsona i Halla (1953)](https://www.xray.cz/xray/csca/kol2011/kurs/Dalsi-cteni/clanky/Williamson-ActaMet-1953-1-22-WH-Plot.pdf).
FWHM i szerokość całkowa nie są wymiennymi definicjami β; program używa
konsekwentnie pierwszej z nich i nie deklaruje pełnej analizy Fourierowskiej.

Wymagane są minimum trzy niezależne refleksy tej samej fazy. Niepewności
β są przenoszone do y i wag 1/u(y)². Niepewności kąta/X oraz korelacje
wspólnego wzorca nie są uwzględnione w tej regresji; wynik jest warunkowy.
Gdy σ nie podano, obliczana jest jawnie oznaczona regresja nieważona.
Raportowane są SE, 95% CI, R² i skorygowane R², RMSE, residua, leverage,
standaryzowane residua, odległość Cooka i pominięcie każdego refleksu
po kolei (dla 3 refleksów takie dopasowania kontrolne mają tylko 2 punkty).

`intercept<=0` blokuje D_WH. Blokuje go również CI wyrazu wolnego
obejmujące zero. Ujemne ε odrzuca interpretację fizyczną.
Progi R²<0.8, |standaryzowane residuum|>3 i Cook>4/n są ostrzeżeniami
diagnostycznymi, nie prawdopodobieństwami poprawności. Wynik może mieć
`REJECTED` i nadal zawierać parametry matematycznego dopasowania do audytu.
Pole `valid_for_quantitative_interpretation` pozostaje False bez osobnej
walidacji eksperymentalnej i fizycznej.

## Parametr sieci kubicznej

`fit_cubic_lattice` dopasowuje wszystkie pozycje do modelu Bragga
`2theta = 2 asin(lambda sqrt(h²+k²+l²)/(2a)) + zero_shift`.
Wymagane są co najmniej 3 rodziny refleksów przy wolnym zero shift
i 2 bez niego. Przyjmuje hkl jako trzy liczby całkowite albo tekst
`111`, `(1,1,1)`; dla wielocyfrowych indeksów należy używać separatorów.
Nie identyfikuje automatycznie hkl/fazy. Zwraca u(a), CI, zero shift,
kowariancję, korelację, predykcje i residua.
u(a) nie obejmuje u(λ) ani systematycznych błędów osi.

## Publiczne API i zapis

```python
fit = fit_peak_profile(raw_x_in_roi, raw_y_in_roi, center_hint=43.4,
                       hkl="111", sigma_y=None, kalpha2_ratio=0.5)
saved_fit = fit.to_dict()        # tylko JSON-safe list/dict/liczby/None
restored = PeakFitResult(**saved_fit)
line_y = evaluate_peak_profile(raw_x_in_roi, restored)

instrument = InstrumentBroadeningModel(**calibration_dictionary)
corrected = correct_instrument_broadening(fit, instrument)
if corrected.fwhm_corrected_deg is not None:
    size = scherrer_size(fit.two_theta_kalpha1_deg,
                        corrected.fwhm_corrected_deg,
                        wavelength_A=fit.wavelength_A,
                        u_beta_deg=corrected.u_fwhm_corrected_deg,
                        u_two_theta_deg=fit.u_two_theta_deg)

wh = fit_williamson_hall(positions, corrected_widths,
                         sigma_beta_deg=width_uncertainties, hkls=hkls)
lattice = fit_cubic_lattice(positions, hkls, fit_zero_shift=True,
                            sigma_2theta_deg=position_uncertainties)
```

W kalibracji wymagane są klucze `two_theta_deg`, `gaussian_fwhm_deg`,
`lorentzian_fwhm_deg`. Opcjonalne: `u_gaussian_fwhm_deg`,
`u_lorentzian_fwhm_deg`, `name`, `standard`, `measurement_date`,
`radiation`, `wavelength_A`, `optics_id`, `source_file`.
Tablice muszą mieć tę samą długość. Nie należy wpisywać umownych
szerokości z przykładu jako rzeczywistego wzorca.

## Wykonana walidacja

`python -m unittest discover -s tests -p test_profile_analysis.py -v`:
**19/19 PASS**, Python 3.13, NumPy 2.5.3, SciPy 1.18.1.

Testy obejmują granice Gaussa/Lorentza, FWHM TCH względem niezależnego
splotu Voigta, normalizację pól, geometrię Bragga dubletu, dopasowanie na
nierównym próbkowaniu z tłem i szumem, brak modyfikacji raw, round-trip
JSON, dekonwolucję składowych, zakaz ekstrapolacji, ujemne składowe,
niepewność Bragga, konwersje Å/nm i stopnie/radiany, Scherrera o zadanym
D, WH o zadanym ε/D, poprawę dzięki ważeniu, granice fizyczne WH,
leverage/LOO oraz odzyskanie a i zero shift.

Nie wykonano porównania z surowymi danymi referencyjnymi 316L/CuNi ani
pomiarami NIST, ponieważ nie były załączone do tej aktualizacji. Liczby
z opisowych raportów nie są zakodowane jako oczekiwane wyniki naukowe.
