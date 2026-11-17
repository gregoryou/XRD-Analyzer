# XRD Analyzer 20.3.5 — COD Offline

## Cel

Wersja 20.3.5 dodaje pobieranie i aktualizowanie metadanych struktur krystalicznych z oficjalnej **Crystallography Open Database (COD)** do lokalnego pliku SQLite zgodnego z istniejącym czytnikiem programu. Nie zmienia sposobu wyznaczania pików, kalibracji instrumentu, geometrii sieci, faz, naprężeń ani szerokości profili.

## Szybka synchronizacja podzbioru COD — zalecana

1. Rozpakuj nową paczkę do osobnego folderu; **nie nadpisuj działającej 20.3.4**. Uruchom `START_XRD_ANALYZER.bat`.
2. Otwórz **Bazy faz → COD .db3**.
3. Wybierz „Pierwiastki (wymagane)” i wpisz np. `Fe Cr Ni`, `Cu Ni`, `Al Mg` albo `Ti Al V`. Zaznaczenie **Tylko podane pierwiastki** wyłącza rekordy o *większej liczbie różnych elementów*; to filtr API, nie pełne sprawdzenie podzbioru składu.
4. Kliknij **Pobierz / aktualizuj COD Online → SQLite**. To wymaga internetu. Jeżeli wyniki przekraczają 30 tys. rekordów albo 48 MB, zawęź filtr (np. więcej pierwiastków) lub wybierz import archiwum.
5. Pomyślnie pobrana baza znajduje się domyślnie w `%LOCALAPPDATA%\XRDAnalyzer\COD_Offline.db3`. Program otworzy ją automatycznie. Później można korzystać z niej **offline**, przez istniejące „Szukaj w COD”, „Tylko zewnętrzna baza COD” lub tryb łączony.
6. Powtórzenie synchronizacji tego samego filtra aktualizuje istniejące rekordy według **COD ID**. Inne rekordy nie są usuwane. To podzbiór, a nie kompletna migawka całej COD.

Wybierz „Numery COD” i podaj do 100 siedmiocyfrowych identyfikatorów, jeżeli interesują Cię konkretne struktury.

## Pełna baza — opcjonalnie

Oficjalna COD zaleca do utrzymywania pełnej, często aktualizowanej kopii **rsync lub Subversion**. W XRD Analyzer dodałem też dwa tryby awaryjne o dużych wymaganiach dyskowych:

- **Importuj lokalne archiwum CIF ZIP** — wskaż pobrane osobno archiwum zawierające pliki `XXXXXXX.cif`. Program strumieniowo czyta CIF, wydobywa metadane przez `gemmi` i partiami zapisuje je do SQLite. Nie wypakowuje całej bazy na dysk. Pojedyncze zbyt duże lub niepoprawne CIF mogą zostać pominięte.
- **Pobierz pełne archiwum COD (duże, opcjonalne)** — po świadomym potwierdzeniu program pobiera archiwum ze źródła `https://www.crystallography.net/archives/cod-cifs-mysql.zip` i próbuje je zaindeksować; rozmiar, dostępność i czas zależą od serwera. Nie uruchamia się automatycznie. W razie błędu sieciowego nie usuwa istniejącej bazy.

Archiwa serwisu mogą być mniej aktualne niż COD online. Pełny snapshot nie jest automatycznie odświeżany przez niewielkie zapytania API; do wiernej synchronizacji pełnej bazy zalecane pozostają `rsync`/`svn`.

## CLI (tryb zaawansowany)

```powershell
cd 'C:\sciezka\XRD_Analyzer_GUI_COMPLETE_V20_3_5_PL'
.\.venv\Scripts\python.exe COD_OFFLINE_CLI.py elements --elements 'Fe Cr Ni' --strict
.\.venv\Scripts\python.exe COD_OFFLINE_CLI.py ids --ids '1000000 2002926'
.\.venv\Scripts\python.exe COD_OFFLINE_CLI.py archive --file 'D:\COD\cod-cifs-mysql.zip'
```

## Dokładność i bezpieczeństwo

- COD podaje parametry sieci krystalicznej i metadane literaturowe. Program **nie otrzymuje automatycznie eksperymentalnych profili intensywności** i nie wykonuje tym samym pełnej analizy Rietvelda.
- Wprowadź wiarygodną chemię próbki i sprawdzaj przypisane fazy, pasujące refleksy oraz przywołaną publikację, zwłaszcza przy fazach o nakładających się pikach.
- SQLite i metadane COD nie zawierają w tym module oryginalnych surowych pomiarów XRD. Przesłana wcześniej rzeczywista seria Si w `tests/data` jest zachowana w pakiecie **lokalnym**; przed osobną publikacją open source wymaga ustalenia praw do dystrybucji.
- Oficjalne dane COD są CC0; przy użyciu naukowym należy wskazywać autorów oryginalnych danych krystalograficznych.
- Serwer może zwrócić odpowiedź inną niż CSV lub ograniczyć zapytania; wtedy program zgłasza błąd i zachowuje dotychczasową bazę. W tej sesji nie przetestowano połączenia z serwerem na żywo.

## Oficjalne źródła

- API COD: https://wiki.crystallography.net/RESTful_API/
- Pobieranie COD (SVN/rsync/archiwa): https://wiki.crystallography.net/howtoobtaincod/
- Licencja CC0 i cytowanie: https://www.crystallography.net/cod/
