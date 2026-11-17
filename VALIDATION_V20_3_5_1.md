# XRD Analyzer 20.3.5.1 — poprawki po walidacji COD na próbce Si

Data: 17.11.2026. Aktualizacja wersji 20.3.5 z zainstalowaną poprawką Pomoc/XRD.

## Zakres

1. **Dublet Kα₁ / Kα₂**: konserwatywne oznaczanie prawdopodobnych par
   tylko po wczytaniu dwóch długości fal z nagłówka pliku. Użytkownik
   może wyłączyć pomijanie kandydatów Kα₂ w rankingu COD. Wszystkie
   zmierzone piki pozostają zapisane w projekcie i tabeli.
2. **Eksport XLSX**: zachowane dotychczasowe arkusze oraz nowe
   „Pochodzenie pikow” (ID COD, grupa przestrzenna, źródło, hkl,
   Δ2θ, status dubletu) i „Fazy COD” (źródło rankingu).
3. **Ranking chemiczny**: poprzednie 108% przedstawiane jest jako
   mnożnik 1,080×. To NIE jest prawdopodobieństwo identyfikacji.
   Sama formuła rankingu COD i obliczania refleksów NIE zmieniła się.
4. **Oznaczenia wykresów**: wyświetlanie rodzin hkl dla kubicznych
   faz COD ujednolicono (np. 022→220). Zbliżone etykiety są
   rozdzielane pionowo; możliwy dublet Kα₂ oznaczany jest „Kα₂?”.
   Oryginalne hkl oraz wyniki w projekcie pozostają bez zmian.
5. **Importer .rd**: po wczytaniu danych V3RD z brakującymi punktami
   względem nagłówka wyświetla ostrzeżenie. Liczba brakujących
   punktów i sposób zakotwiczenia osi są zachowywane jako metadane.
   Oryginalny algorytm parsowania .rd pozostaje bez zmian.
6. **Starsze projekty**: jeśli oryginalny .rd istnieje i jego oś
   odpowiada projektowi, można odtworzyć źródłowe metadane Kα₂,
   bez zmiany zapisanych wartości pomiarowych.

## Ograniczenia

- Dublety są wykrywane heurystycznie i wymagają kontroli operatora;
  przy niepewnej geometrii/widmie nie należy usuwać ani dopasowywać
  automatycznie zidentyfikowanych składowych.
- Wyszukiwanie COD nadal jest dopasowaniem położeń, NIE analizą
  ilościową Rietvelda i nie wykorzystuje modelu intensywności.
- Wynik rankingu 100% nie jest naukowym stopniem pewności.
- Testy kodu i parsera wykonano bez GUI, ponieważ środowisko
  przygotowania nie ma PyQt6. GUI i XLSX należy przetestować
  na docelowym Windows (skrypt TEST_2035_1_WINDOWS.ps1).

## Bezpieczeństwo

Instalator nie dotyka .rd, .db3, .xrdproj ani plików wynikowych.
Sprawdza sumy SHA-256 źródłowych plików oraz payloadu, zatrzymuje się
przy nieoczekiwanych modyfikacjach i wykonuje kopie zapasowe.
