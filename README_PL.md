# XRD Analyzer 20.3.5.1 — opis programu (PL)

XRD Analyzer to program badawczy do analizy **dyfrakcji rentgenowskiej (X-ray Diffraction, XRD)**. Metoda XRD mierzy intensywność promieniowania ugiętego przez materiał w funkcji kąta 2θ i pozwala badać strukturę krystaliczną oraz identyfikować kandydatów na fazy.

Program odczytuje dyfraktogramy, koryguje tło, wykrywa refleksy, pozwala analizować ich położenie, szerokość i profile, porównywać próbki, dopasowywać struktury z COD oraz tworzyć wykresy do dokumentacji i publikacji. Dostępne są także moduły analizy profili, Scherrera, Williamsona–Halla i eksperymentalne workflow sin²ψ. Metody mają ograniczenia i wymagają weryfikacji naukowej.

**Import:** `.udf`, `.rd` (zweryfikowany wariant V3RD), `.xrdml`, wybrane `.xml`, proszkowy `.cif` z intensywnościami, `.csv`, `.txt`, `.dat`, `.asc`.

**Eksport:** CSV, XLSX, PNG, TIFF, PDF, SVG, EPS. **Projekty:** `.xrdproj` i zgodne starsze `.json`.

**Baza struktur:** bezpłatna [COD](https://www.crystallography.net/cod/). Gotowe pliki SQLite `.db3` bazujące na COD można znaleźć w zasobach [Profex](https://www.profex-xrd.org/download/); wymagany odpowiedni schemat danych. **Baza COD nie jest częścią repozytorium.**

**Start Windows:** uruchom `START_XRD_ANALYZER.bat` po rozpakowaniu do nowego katalogu. Program tworzy środowisko `.venv` i doinstalowuje zależności.

**Ograniczenia:** Wynik rankingu COD nie jest prawdopodobieństwem identyfikacji ani ilościową analizą fazową Rietvelda. Wyniki wymagają kontroli danych, kalibracji i interpretacji przez badacza. Kod źródłowy: GPL-3.0-only; wykorzystanie AI podczas opracowywania udokumentowano osobno.
