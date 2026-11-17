# XRD Analyzer 20.3.0 PL — nowe funkcje

Rozpakuj cały ZIP do nowego katalogu i uruchom `START_XRD_ANALYZER.bat`.
Nie kopiuj samego głównego .py: potrzebuje czterech nowych modułów z paczki.
Zachowaj kopię dotychczasowych projektów. Format projektu 13 czyta starsze
projekty; ponowne zapisanie w starej aplikacji może usunąć nowe ustawienia.

## Profile, FWHM, Scherrer i Williamson–Hall

1. Wczytaj dyfraktogram i wykryj piki jak dotychczas.
2. Otwórz **Analiza naukowa → Profile / Scherrer / WH / sin²ψ / Wilkens…**.
3. Sprawdź λ1, λ2 i stosunek pól Kα2/Kα1. Dla pojedynczej linii ustaw stosunek 0.
4. W tabeli podaj fazy/hkl i zakres każdego refleksu. Usuń lub odznacz
   nakładające się refleksy. Zakres powinien objąć dublet oraz tło po obu stronach.
5. Naciśnij **Dopasuj zaznaczone profile**. Wybór wiersza pokazuje pomiar,
   dopasowanie i residua. Diagnostyka zawiera niepewności, 95% CI i ostrzeżenia.
6. Wczytaj zmierzony profil instrumentalny IPF. Wtedy pojawią się szerokości
   skorygowane i warunkowy Scherrer. Zaznacz ≥3 różne hkl jednej fazy do WH.
7. Przycisk parametru sieci dopasowuje a i zero shift dla struktury kubicznej.

Wzorzec IPF można wyeksportować z fitów pomiaru wzorca (minimum dwa kąty)
przyciskiem **Eksportuj fity jako wzorzec**, a następnie wczytać dla próbki.
Plik JSON zawiera metadane i tablice; schemat opisano w SCIENTIFIC_PROFILE.md.
CSV wymaga kolumn `two_theta_deg, gaussian_fwhm_deg, lorentzian_fwhm_deg`,
opcjonalnie `u_gaussian_fwhm_deg, u_lorentzian_fwhm_deg`. Akceptowane separatory
to przecinek, średnik i tabulator; liczby używają kropki dziesiętnej.
CSV przypisuje bieżącą λ1; sprawdź zgodność konfiguracji. Brak ekstrapolacji IPF.

Tabela detekcji pokazuje teraz **szerokość przy połowie prominencji**.
Ta wielkość nie jest modelowym FWHM. Scherrer określa domeny koherentnego
rozpraszania, nie ziarna metalograficzne. WH pokazuje apparent microstrain;
niefizyczny albo nierozdzielony intercept blokuje D_WH.

## Naprężenia sin²ψ

Użyj zakładki **Naprężenia sin²ψ**. Wpisz serię pomiarów lub wczytaj CSV:
`psi_deg, phi_deg, two_theta_deg, u_two_theta_deg, hkl, phase`.
Opcjonalna kolumna `sample` identyfikuje próbkę. Wszystkie pozycje używają
λ1 z zakładki profili. Dodaj bieżący fit albo zbierz wybrane hkl z dopasowań
wszystkich wczytanych próbek; uzupełnij ψ/φ i usuń powtórzenia.
Wymagane są rentgenowskie S2/2 w 1/MPa, d0 w Å i ich niepewności.
Źródło XEC oraz zgodność fazy/hkl są jawne. Zwykły pojedynczy skan nie daje MPa.

## Wilkens — funkcja eksperymentalna

Wybierz ≥3 całe, rozdzielone profile jednej fazy; podaj w tabeli C_hkl,
a w zakładce Wilkensa b [nm], strukturę, źródła IPF/kontrastów i układ poślizgu.
Program dopasuje ograniczony model Wilkensa z lognormalnymi domenami,
wyświetli residua i warunkowe ρ [m⁻²]. Nie jest to pełny CMWP.
Szczegóły i niewdrożone elementy podano w SCIENTIFIC_ADVANCED.md.

## Wygładzanie wykresu i geometria

W Kreatorze wykresów otwórz konfigurację kolejności materiałów.
Każda próbka ma **Wygładź** oraz **Siła [%]**; pole siły jest aktywne po
zaznaczeniu checkboxa. 0% zachowuje przebieg, 100% używa filtru SG;
wartość pośrednia to `(1-p/100)*y + (p/100)*SG(y)`.
Domyślny zakres okna to 0.25° 2θ, wielomian stopnia 2, okno nieparzyste
ograniczone długością danych. Dla nierównej osi filtr używa pomocniczej
siatki równomiernej. Dla zbyt krótkich lub błędnych danych zwraca kopię bez zmian.
Wygładzanie wpływa wyłącznie na wykres i jest zapisywane osobno dla każdej próbki.

Główne okno zapisuje normalną geometrię i maksymalizację lokalnie w QSettings.
Po zmianie monitora mieści się w dostępnym obszarze. Przy małym ekranie dostępne
są paski przewijania. Projekt nie przenosi geometrii między komputerami.

## Zapis i interpretacja

Ustawienia oraz wyniki nowych analiz zapisują się w metadanych próbki w .xrdproj.
JSON zachowuje pełne diagnostyki; CSV eksportuje tabelę profili i ostrzeżenia.
Zmiana λ unieważnia zależne wyniki. Edycja wejść analizy wymaga przeliczenia.
Indeks dowodów fazowych ma skalę 0–100, a nie prawdopodobieństwo procentowe.
Statusy wyników są jawne; żaden test syntetyczny nie zastępuje walidacji
na rzeczywistych pomiarach i wzorcu tej samej aparatury.

## Zmiany w 20.3.1

- Pliki CSV/TXT/DAT/ASC: obsługiwane są separatory przecinek, średnik, tabulator i spacja oraz przecinek
  dziesiętny (np. `35,02;1000`). Wiersze nagłówka są pomijane; malejąca oś 2θ jest odwracana.
- Menu **Plik / Analiza / Analiza naukowa / Wykres / Bazy faz / Pomoc**. Skróty: Ctrl+O, Ctrl+S, Ctrl+Shift+S,
  Ctrl+Shift+O (wczytaj pliki), F5 / Shift+F5 (analizuj bieżący / wszystkie), Ctrl+K (kreator wykresów),
  Ctrl+E (zapisz wykres), Ctrl+T (eksport tabeli), Ctrl+Shift+A (analiza naukowa), F1 (pomoc).
- Pola liczbowe zatwierdzaj klawiszem Enter lub przejściem do innego pola; dotyczy zwłaszcza λ,
  której zmiana unieważnia wyniki identyfikacji faz i analizy naukowej.
- W oknie analizy naukowej dziesiętny przecinek jest akceptowany; wyniki podawane są jako „wartość ± niepewność”.
- Nieoczekiwane błędy nie zamykają programu; zapisują się w `~/XRDAnalyzerData/XRD_Analyzer_error.log`.

## Zmiany w 20.3.2

- **Kreator wykresów → „Konfiguracja kolejności”**: nowa kolumna *Metoda wygładzania*. „Z zachowaniem pików
  (zalecane)” wygładza tło jak dotąd, ale nie spłaszcza wąskich refleksów; „Klasyczny Savitzky–Golay” odtwarza
  wygląd z wersji do 20.3.1. Projekty zapisane wcześniej zachowują metodę klasyczną.
- Wygładzanie jest nadal tylko wizualne — nie zmienia danych, pików ani wyników analizy.
- **Linia tła (AsLS)**: w 20.3.2 wprowadzono skalowanie sztywności z krokiem pomiaru z błędnym odniesieniem 0,02°
  — patrz 20.3.3.

## Zmiany w 20.3.3

- **Linia tła (AsLS)** ma taką samą sztywność w stopniach 2θ dla każdego kroku pomiaru, z odniesieniem 0,05°.
  Dla skanów 0,05° wyniki są takie same jak w 20.3.1; drobniejszy krok nie pochłania już szerokich refleksów.
- **Pliki Philips `.rd` (V3RD)** można wczytać jak pozostałe formaty. Odczyt jest sprawdzany spójnością nagłówka
  z liczbą punktów; przy niepowodzeniu program podaje powód zamiast zgadywać. Po wczytaniu sprawdź zakres 2θ i
  długość fali (uwagi próbki). Obsługa V5RD wymaga pliku wzorcowego.

## Zmiany w 20.3.4

### Jak skorzystać z kalibracji wzorcem Si

1. Wczytaj do jednego projektu skany próbek i skan wzorca Si (formaty CIF, UDF, XRDML, CSV oraz `.rd`).
2. Analiza naukowa → **Kalibracja z wzorca Si…** (Ctrl+Shift+K). Wybierz skan Si, sprawdź λ Kα2, parametr sieci Si
   (domyślnie 5,43102 Å — wartość literaturowa; wpisz wartość z certyfikatu, jeśli ją masz, i zaznacz to) i naciśnij
   „Dopasuj refleksy wzorca”.
3. Przejrzyj tabelę i wykres. Odznacz refleksy, które się nie dopasowały. Wybierz model korekty (domyślnie „Tylko zero”).
4. „Zastosuj do próbek”. Próbki z tą samą λ dostają poprawione 2θ i d; próbki bez pików zostaną skorygowane przy analizie
   (F5 / Shift+F5). Podsumowanie próbki pokazuje zastosowaną poprawkę.
5. W oknie analizy naukowej profile są dopasowywane na osi skorygowanej, a profil aparatury z wzorca można wczytać
   przyciskiem „IPF z wzorca Si” (ustawia się sam, jeśli próbka nie ma jeszcze IPF).

Jedna kalibracja przypada na projekt: nowa zastępuje poprzednią. „Usuń kalibrację ze wszystkich” przywraca zmierzone
położenia pików.

### Pozostałe zmiany

- **Pliki `.rd`** są wczytywane (sprawdzone na rzeczywistym pliku wzorca Si; patrz CHANGELOG).
- **Wykrywanie pików**: szerokie piki i piki na zboczu silnego refleksu są wykrywane, a dotychczasowe wyniki się nie
  zmniejszają.
- **Ferryt z jednego refleksu** jest oznaczany prowizorycznie (znak „?”), jeśli pasuje do refleksu 110 i skład pozwala.
