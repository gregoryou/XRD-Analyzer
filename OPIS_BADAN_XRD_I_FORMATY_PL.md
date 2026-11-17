# XRD Analyzer 20.3.5 — badania XRD, formaty i baza COD

**Cel programu.** XRD (X-ray Diffraction) to dyfrakcja rentgenowska. Analiza
położeń, szerokości i intensywności refleksów umożliwia badanie faz
krystalicznych, zmian parametrów sieci, rozmiaru krystalitów i mikroodkształceń.
XRD Analyzer służy do obróbki dyfraktogramów, porównywania próbek oraz
przygotowywania wykresów (intensywność–2θ, profile pików, Williamson–Hall,
sin²ψ) i zestawień do raportów oraz artykułów naukowych.

## Formaty wejściowe

| Rozszerzenie | Rodzaj danych / uwagi |
|---|---|
| `.udf` | PANalytical UDF (RawScan) |
| `.rd` | Philips / PANalytical, obsługiwany wariant V3RD |
| `.xrdml`, `.xml` | Zgodne pliki pomiarowe XRDML |
| `.cif` | pdCIF z profilem intensywności; zwykły strukturalny CIF nie jest dyfraktogramem |
| `.csv`, `.txt`, `.dat`, `.asc` | Dwukolumnowe dane tekstowe 2θ [°], intensywność |
| `.xrdproj`, `.json` | Projekty XRD Analyzer (to nie są surowe dane pomiarowe) |
| `.db3`, `.sqlite`, `.sqlite3`, `.db` | Zewnętrzna baza faz COD w zgodnym schemacie SQLite |

**Nie należy utożsamiać rozszerzeń z uniwersalną kompatybilnością.** Wersja
20.3.5 nie importuje bezpośrednio `.raw`, `.brml` ani `.xy`. Plik o zgodnym
rozszerzeniu, lecz niezgodnej strukturze może zostać odrzucony.

## Eksport

- Wykresy: `.png`, `.tif`, `.tiff`, `.pdf`, `.svg`, `.eps`.
- Tabele danych i pików: `.csv`, `.xlsx`.
- Projekty: `.xrdproj` i `.json`.

## COD — źródła danych

- Oryginalna otwarta baza danych: https://www.crystallography.net/cod/
- Gotowe pliki `.db3` dla Profex: https://www.profex-xrd.org/download/
- Przykładowe wydanie bazy: `cod-260101.db3` (01.01.2026).

Baza przygotowana dla Profex może wymagać konwersji przed wczytaniem do
XRD Analyzer. Alternatywnie w `Bazy faz → COD .db3` można utworzyć własną
bazę SQLite z danych pobranych przez oficjalne API COD.
Przy publikacji wyników należy wskazać oryginalne źródła strukturalne.

## Ograniczenia naukowe

Jest to oprogramowanie badawcze, nie certyfikowany system metrologiczny.
Identyfikacja faz, wynik Scherrera / Williamsona–Halla oraz analiza sin²ψ
wymagają właściwych danych wejściowych, korekcji aparaturowych i weryfikacji.
Struktury COD nie dostarczają automatycznie zmierzonych intensywności
refleksów proszkowych.

**W aplikacji:** `Pomoc → Badania XRD — analiza, wykresy i formaty plików` oraz
`Bazy faz → COD .db3 → Skąd pobrać COD?`.
