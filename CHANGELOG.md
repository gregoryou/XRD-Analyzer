# Historia zmian

## 20.3.5 PL — 2026-11-17: COD Offline

- Dodano moduł `xrd_cod_offline.py` oraz panel w `Bazy faz → COD .db3`.
- Synchronizacja metadanych COD przez oficjalne API (`format=csv`) z filtrami pierwiastków i numerów COD.
- Zapisywanie i aktualizacja rekordów w kompatybilnej lokalnej bazie SQLite, bez duplikatów.
- Zabezpieczono transfer: limity 48 MB i 30 tys. rekordów na odpowiedź, jawne anulowanie i zachowanie danych przy błędach.
- Opcjonalny import bardzo dużego archiwum CIF `.zip` bez wypakowywania, z wyliczaniem metadanych w `gemmi` i zapisem partiami.
- Możliwość pobrania całego oficjalnego archiwum COD dopiero po świadomym zatwierdzeniu transferu.
- Dodatkowa konsola `COD_OFFLINE_CLI.py` i testy bez połączenia z internetem.
- Nie zmieniono algorytmów naukowych ani formatu projektów `.xrdproj`.
- Przegląd ograniczeń i linków: `COD_OFFLINE_INSTRUKCJA_PL.md`.

## 20.3.4 PL — 2026-10-17

### Kalibracja z wzorca krzemowego (nowa funkcja)

- **Analiza naukowa → Kalibracja z wzorca Si… (Ctrl+Shift+K).** Skan Si zmierzony razem z próbkami (wczytany do tego
  samego projektu) służy do poprawienia wyników pozostałych próbek. Program dopasowuje profile wszystkich refleksów Si
  w zakresie skanu i wyznacza (1) błąd osi 2θ aparatu: Δ(2θ) = p0 + p1·cosθ + p2·sin2θ (zero, przesunięcie próbki,
  transparencja; człony wybierane ręcznie lub kryterium AICc), (2) szerokość linii aparatury (profil IPF).
- Korekta położeń: 2θ_popr = 2θ − Δ(2θ); d jest przeliczane. Dotyczy próbek o tej samej długości fali, jest zapisana
  w projekcie, a wynik dla każdej próbki jest zapisywany w sposób odwracalny (powtórna analiza lub nowa kalibracja nie
  nakłada poprawki dwa razy; zmiana λ próbki wycofuje poprawkę i oznacza to w podsumowaniu).
- Okno analizy naukowej dopasowuje profile na osi skorygowanej, więc parametr sieci, WH i naprężenia dostają poprawione
  położenia, a przycisk „IPF z wzorca Si” (automatycznie, jeśli próbka nie ma jeszcze IPF) wprowadza profil aparatury
  zmierzony na wzorcu do korekcji Scherrera i Williamsona–Halla.
- Zastrzeżenia widoczne w oknie i w zapisanych uwagach: wzorzec niecertyfikowany (a = 5,43102 Å z literatury, błąd a
  przesuwa człon zerowy); człon cosθ to przesunięcie próbki WZORCA — przenosi się na próbki tylko przy tej samej wysokości
  mocowania, dlatego domyślnie stosowane jest samo zero; profil aparatury nie jest ekstrapolowany poza zakres refleksów Si.
- Dla wzorca Si dostarczonego z aparatem (PW3064, Cu Kα): Δ(2θ) = −0,0085 ± 0,0009° (tylko zero, RMS reszt 0,0024°,
  8 refleksów), czyli próbki przesuwają się o +0,0085° (d o 0,0001–0,001 Å). Linie wzorca są czysto lorentzowskie:
  FWHM 0,117° przy 28° do 0,177° przy 107°.
- Położenia na liście pików to nadal maksimum obwiedni Kα1+Kα2 (kilka setnych stopnia powyżej Kα1); kalibracja usuwa
  błąd aparatu, nie tę różnicę. Dopasowania w analizie naukowej są liczone dla samego Kα1.

### Import .rd sprawdzony na rzeczywistym pliku

- Plik wzorca Si (PW3064, PW3050/60 Theta/Theta) potwierdził układ: nagłówek 810 B, zliczenia 16-bitowe, długości fal
  jako float64 (przesunięcia 94/102/110), krok/start/koniec jako float64 (214/222/230). Nagłówek opisuje 2992 punkty
  (start 15,046°), plik zawiera 2711, czyli brakuje 281 punktów na początku skanu. Refleksy Si wskazały, że zapisane
  punkty kończą się na kącie „koniec” — czytnik kotwiczy wtedy oś na końcu i ostrzega (osiem refleksów w ±0,013° od
  teorii, a = 5,43132 ± 0,00004 Å). Wniosek opiera się na JEDNYM pliku; V5RD i inne aparaty pozostają niesprawdzone.
- Plik wzorca jest dołączony jako test regresji (tests/data/Si_polycrystal_PW3064.rd).

### Wykrywanie pików

- Szerokie piki i piki na zboczu silnego refleksu nie były wykrywane, bo próg SNR dotyczył prominencji pojedynczego
  punktu. Teraz pik jest przyjmowany także wtedy, gdy ma prominencję ≥ 5σ, istotność uwzględniającą szerokość
  (prominencja·√(0,752·szerokość w punktach)/σ) ≥ min SNR oraz istotność filtru dopasowanego w wielu skalach ≥ min SNR.
  Dotychczas wykrywane piki zawsze pozostają. Na czystym szumie wychodzi ok. 1 fałszywy pik na 200 skanów.
- Na Twoich próbkach: nowe piki to wyłącznie 44,52° w material_podstaw_316L i 96,11° (refleks 222 austenitu) w drucie 316L.

### Oznaczanie ferrytu z jednego refleksu

- Pik silny (SNR ≥ 10), którego nie wyjaśnia żadna faza, a który leży przy refleksie 110 ferrytu/martenzytu stali
  zgodnej z dozwolonymi pierwiastkami, dostaje etykietę prowizoryczną (np. „δ-ferrite?” / „… (niepotwierdzony, 1 refleks)”,
  źródło „internal-tentative”). Reguła potwierdzania faz BCC wieloma refleksami nie została zmieniona. Notatka podaje,
  czy refleks 200 jest widoczny: w material_podstaw_316L zmierzono 420 zliczeń przy ok. 567 oczekiwanych, co zgadza się
  z ferrytem. W pozostałych próbkach projektu nic nie jest oznaczane.

## 20.3.3 PL — 2026-10-15

### Korekta błędu z 20.3.2 (linia tła AsLS)

- Odniesieniem dla skalowania sztywności tła jest krok 0,05° (skany PANalytical, dla których λ=2e6 sprawdzono
  na rzeczywistych danych), a nie 0,02°. Dla skanów o kroku 0,05° wynik jest taki sam jak w 20.3.1: piki z
  projektów użytkownika zgadzają się co do położenia, a wysokości różnią się o najwyżej 0,16% (UDF o kroku
  0,0501° ma λ skalowaną o 1%). Skany o drobniejszym kroku nadal nie pochłaniają szerokich refleksów.
- Błąd wyszedł dopiero na rzeczywistych danych; testy syntetyczne z 20.3.2 go nie wykrywały, bo przyjmowały
  to samo błędne założenie. Test regresyjny jest teraz oparty na kroku 0,05°.

### Import plików Philips .rd (V3RD) — NIEPRZETESTOWANY NA PRAWDZIWYM PLIKU

- Nowy moduł xrd_rd_import.py i rozszerzenie *.rd w oknie wczytywania. Format jest binarny i nie miałem ani
  specyfikacji, ani pliku wzorcowego, dlatego czytnik niczego nie zgaduje: oś kątowa jest przyjmowana tylko
  wtedy, gdy trzy sąsiadujące liczby w nagłówku (start, krok, koniec) zgadzają się z liczbą punktów w pliku,
  a dane wyglądają jak zliczenia detektora. W przeciwnym razie program zgłasza błąd z instrukcją.
- Długość fali i anoda są szukane w nagłówku po znanych liniach (Cu, Co, Cr, Fe, Mo, Ni, Ag); brak = λ domyślna
  z ostrzeżeniem w uwagach próbki.
- V5RD jest rozpoznawany, ale nieobsługiwany (układ nieznany bez próbki). `python xrd_rd_import.py plik.rd`
  wypisuje diagnostykę do zgłoszenia.
- Testy (15) używają plików syntetycznych zbudowanych z założonego układu, więc potwierdzają logikę walidacji i
  obsługę błędów, a nie zgodność z plikami z aparatu.

### Wnioski z rzeczywistych danych (nie zmieniono kodu; opis w validation/dane-rzeczywiste-20.3.3.txt)

- Szeroki garb tła (próbka material_podstaw_316L, 20–35°) jest wykrywany jako „pik” 26,5° o szerokości 4,2° i
  dostaje identyfikację fazową; brak ustawienia maksymalnej szerokości piku.
- Pik na zboczu silnego refleksu (44,5° obok 43,7°) nie jest wykrywany, bo próg SNR dotyczy prominencji
  (742 zliczeń = 6,4σ < 8σ). Obniżenie „Min. sygnał/szum” do 5 go wykrywa.
- Położenie piku to maksimum punktu pomiarowego; dla szerokich, słabych refleksów różni się od dopasowania
  profilu nawet o 0,2°, co przy tolerancji 0,35° potrafi zmienić identyfikację (311 w próbce 316L przypisany
  do Cr z COD).
- Wygładzanie z zachowaniem pików: na skanach 0,05° praktycznie bez różnicy względem klasycznego (piki
  szerokie na 6–25 punktów); na wzorcach LaB6 i NIST 660 klasyczny filtr obniża piki do 67–84%, nowy zachowuje
  98–103%.

## 20.3.2 PL — 2026-10-14

### Logika obliczeń

- Wzory sprawdzone niezależnie (wartości liczone osobno lub obliczenia zwrotne na danych syntetycznych):
  pseudo-Voigt TCH, Scherrer wraz z niepewnością, Williamson–Hall, parametr sieci i przesunięcie zera,
  sin²ψ, niepewność d z prawa Bragga, funkcja Wilkensa, położenie Kα2. Nie wymagały zmian.
- Poprawione: sztywność tła AsLS była zadana w punktach, a nie w stopniach. Przy małym kroku pomiaru tło
  „zjadało” szerokie refleksy (krok 0,005°: szeroki pik odzyskany w 28% zamiast ok. 90%; krok 0,002°: 14%
  zamiast 85%). Parametr λ skaluje się teraz z krokiem jako (krok odniesienia / krok)^4.
  **KOREKTA (20.3.3): w wydaniu 20.3.2 krokiem odniesienia było błędnie 0,02°. Na rzeczywistych skanach o
  kroku 0,05° (316L, CuNi, Monel) tło unosiło się pod piki, zmieniało ich wysokości o 2–10% i w jednej
  próbce usuwało pik. Wydanie 20.3.2 nie powinno być używane do analizy; poprawne odniesienie to 0,05°.** Zapisane piki w projekcie pozostają, dopóki analiza nie zostanie powtórzona;
  linia tła jest przeliczana przy otwarciu projektu.
- Uwagi bez zmian w kodzie (decyzje projektowe do rozważenia): pozycja piku to maksimum punktu pomiarowego
  (rozdzielczość = krok), skala prominencji liczona od najwyższego punktu (pojedynczy skok sygnału zawyża
  próg), jedno odniesienie fazowe może zostać dopasowane do dwóch pików, a parametr sieci w identyfikacji
  faz jest średnią ważoną intensywnością (pełne dopasowanie jest w Analizie naukowej).

### Wygładzanie linii w Kreatorze wykresów

- Wada: klasyczny filtr miał stałe okno 0,25° i stopień 2, więc piki węższe niż okno były spłaszczane
  (przy 100%: pik o szerokości ok. 0,06° tracił 30–60% wysokości, ok. 0,12° — ok. 7%). Przy kroku 0,05°
  okno miało tylko 5 punktów.
- Nowa metoda „Z zachowaniem pików” (domyślna dla nowych wierszy). Tło i szum są wygładzane tak samo jak
  dotąd, ale w sąsiedztwie wąskich pików filtr jest dobierany tak, by zachować co najmniej 99% wysokości
  (filtr 4. stopnia o możliwie szerokim oknie albo wartości zmierzone dla pików o szerokości ok. 3 punktów).
  Szerokie piki są traktowane identycznie jak w filtrze klasycznym. Poziom szumu jest szacowany lokalnie,
  więc metoda działa także przy szumie Poissona. Dane bez szumu nie są zmieniane. Wynik nie wychodzi poza
  zakres danych; procent siły nadal oznacza udział krzywej wygładzonej.
- Pomiar na syntetycznych dyfraktogramach (100%): wysokość wąskich pików 0,98–1,00 zamiast 0,40–0,94;
  błąd względem krzywej bez szumu 11–21 zamiast 42–1430 zliczeń; szum tła taki sam.
- Nowa kolumna „Metoda wygładzania” w oknie kolejności materiałów. Projekty zapisane wcześniej zachowują
  metodę klasyczną, aby istniejące rysunki wyglądały tak samo — wybierz „Z zachowaniem pików”, aby
  zastosować poprawę.

### Testy

- Nowe: tests/test_smoothing_and_baseline.py (18 testów, uruchomione) oraz 2 testy Qt w
  tests/test_plot_smoothing_geometry.py (NIEURUCHOMIONE w środowisku przygotowania poprawek).
- Wynik: validation/tests-20.3.2-bez-qt.txt (112 testów bez Qt: 111 OK, 1 pominięty).

## 20.3.1 PL — 2026-10-14

### Błędy danych i obliczeń

- Import CSV/ASCII: pliki rozdzielane przecinkiem („2θ,intensywność”) były czytane z błędnymi
  intensywnościami (1000 → 0,1) lub z gubieniem wierszy. Nowy parser obsługuje przecinek, średnik,
  tabulator i spację oraz przecinek dziesiętny; wiersze nagłówka zawierające liczby nie są danymi;
  malejąca oś 2θ jest odwracana (inne uporządkowania nie są zmieniane).
- Pole λ w oknie głównym zatwierdza wartość po zakończeniu edycji. Wcześniej każda wpisana cyfra
  przeliczała d i kasowała wyniki identyfikacji faz oraz analizy naukowej.
- Pozostałe pola liczbowe detekcji i etykiet również zatwierdzają się po edycji (mniej zapisów projektu).

### Odporność

- Globalna obsługa wyjątków: błąd w funkcji obsługi zdarzeń nie kończy programu; szczegóły trafiają do
  ~/XRDAnalyzerData/XRD_Analyzer_error.log.
- Eksport tabeli i zapis wykresu: błąd zapisu (np. plik otwarty w Excelu) daje komunikat; xlsx zawiera
  liczby zamiast tekstu; rozszerzenie wynika z wybranego filtru.
- Poprawny komunikat po anulowaniu analizy wsadowej.
- Brak uprawnień do katalogu profilu nie blokuje startu programu.

### Interfejs

- Menu Plik / Analiza / Analiza naukowa / Wykres / Bazy faz / Pomoc ze skrótami klawiszowymi, pasek stanu,
  okno „Skróty klawiszowe i wskazówki”, „O programie”.
- Przyciski w trzech logicznych rzędach z podpowiedziami; zapamiętywanie ostatniego katalogu; domyślne
  nazwy plików z nazwy próbki/projektu.
- Tabela pików sortuje liczby numerycznie (wcześniej alfabetycznie).
- Lewy panel: zawijanie etykiet i krótsze listy rozwijane, aby nie wymuszać poziomego przewijania
  (zmiana układu niesprawdzona wizualnie).
- Okno analizy naukowej: przecinek dziesiętny w tabelach i CSV, czytelne komunikaty dla brakujących lub
  błędnych wartości (C_hkl, ψ, d0, zakres), krótsze liczby i zapis „wartość ± niepewność” zamiast notacji
  wykładniczej, podpowiedzi do kolumn i przycisków; Enter w polu nie uruchamia przycisków.
- Polskie komunikaty błędów parserów plików i bazy COD.

### Testy

- Nowe: tests/test_import_and_formatting.py (20 testów, uruchomione) i tests/test_ui_usability.py
  (7 testów Qt, NIEURUCHOMIONE w środowisku przygotowania poprawek).
- Wynik: validation/tests-20.3.1-bez-qt.txt. Pełny unittest z Qt należy uruchomić na docelowym komputerze.

## 20.3.0 PL — 2026-10-11

### P0 — obliczenia i semantyka

- Odstęp pików w rzeczywistych stopniach 2θ zamiast mediany dx.
- Half-prominence width w tabeli/projekcie; oddzielne modelowe FWHM.
- Sprzężony Kα1/Kα2 + pseudo-Voigt/TCH, IPF, Scherrer i WH z niepewnościami oraz diagnostyką.
- Wielorefleksowy fit a i zero shift struktury kubicznej.

### P1 — projekty

- Wersja schematu 13, wsteczny odczyt, zapis nowych wejść/wyników i smoothingu.
- Ochrona niezapisanej pracy przy otwieraniu projektu, rollback błędnego odczytu.
- Zmiana λ unieważnia nowe wyniki naukowe.

### P2 — nowe funkcje

- Okno Analiza naukowa z profilami/residuami, eksportem JSON/CSV i kalibracją IPF.
- Naprężenia sin²ψ, XEC/d0, kontrola ±ψ, import CSV i zbieranie hkl z próbek.
- Eksperymentalny ograniczony Wilkens whole-profile; nie jest pełnym CMWP.
- Indeks dowodów fazowych 0–100, jawnie nieprobabilistyczny.

### P3 — GUI

- Wygładzanie wizualne osobno dla próbki; procent jako mieszanka raw/SG.
- Odtwarzanie geometrii/maksymalizacji, paski przewijania, przyciski w 3 rzędach.
- Walidacja uszkodzonych opcjonalnych pól Kreatora.

### Testy

- 65 nowych testów + 36 wcześniejszych = 101/101 PASS; kompilacja, Gemmi/COD, Qt offscreen.
- Instrukcja NOWE_FUNKCJE_PL.md, dwa opisy naukowe, projekt DEMO_SYNTHETIC i raport audytu.
- Brak walidacji laboratoryjnej, pełnego CMWP i ręcznej sesji GUI Windows.


## 20.2.1 PL — 2026-10-03

### P1 — Plot Creator

- Zapis projektu pobiera aktualny stan również z otwartego Kreatora.
- Zamknięcie przyciskiem, X i Escape synchronizuje stan; odtwarzanie blokuje pośrednie sygnały i odrysowania.
- Ostatnia pozycja przeciągniętej legendy jest przechwytywana bez dodatkowego odświeżania. Osobne panele zachowują własne pozycje legend.
- Dodano trwałe dataset_id. Zmiana nazwy, zbieżne nazwy i przeniesienie projektu nie zmieniają powiązań ustawień.
- Zachowano odczyt wcześniejszych projektów ZIP/JSON, separatorów i jednoznacznych starych kluczy etykiet.
- Nieudany zapis zachowuje zmiany w pamięci i stary plik, sygnalizuje błąd i pozwala anulować zamknięcie.
- Samo otwarcie/zamknięcie bez edycji nie zapisuje projektu ponownie.

### P0 — obliczenia

- Scentralizowano zmianę λ w profilu aparatury i edytorach metadanych. Przeliczane jest d_hkl; nieaktualne fazy/hkl/Δ2θ i dopasowania zostają unieważnione z komunikatem.
- Zachowano zmierzone pozycje, intensywności, szerokości i własne teksty pików przy zmianie λ.
- Szerokość piku jest przeliczana przez interpolację rzeczywistej osi 2θ. Usunięto użycie mediany kroku całego pomiaru dla FWHM.
- Osie malejące, powtórzone, niefinitywne i o niewłaściwej długości są odrzucane przed zmianą wyników.

### Walidacja i stan

- 27 nowych testów + 9 zachowanych = 36/36 PASS.
- Restart testowany w nowym procesie Python/QApplication z rzeczywistym .xrdproj.
- Dodano cztery pliki kontynuacji, validation/ i zachowany rozszerzony plan.
- APP_VERSION=20.2.1; PROJECT_VERSION=12 z opcjonalnymi polami dataset_id i manual_legend_positions.
- Oryginalny ZIP pozostał nienaruszony. Nie wykonano pełnego audytu naukowego.

## Stan źródłowy 20.2

Baseline tej sesji: składnia/import/zależności PASS, 9/9 istniejących testów PASS.
Historyczne VALIDATION_* zachowano; nie są wynikami testów nowej wersji.
