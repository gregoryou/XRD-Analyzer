# -*- coding: utf-8 -*-
"""Help-only XRD and COD documentation. No scientific analysis algorithms."""
from __future__ import annotations

# Exactly the suffixes accepted by load_dataset() in XRD Analyzer 20.3.5.
IMPORT_SUFFIXES = ('.cif', '.udf', '.rd', '.xrdml', '.xml', '.csv', '.txt', '.dat', '.asc')
PROJECT_SUFFIXES = ('.xrdproj', '.json')
COD_SUFFIXES = ('.db3', '.sqlite', '.sqlite3', '.db')
EXPORT_PLOT_SUFFIXES = ('.png', '.tif', '.tiff', '.pdf', '.svg', '.eps')
EXPORT_TABLE_SUFFIXES = ('.csv', '.xlsx')

XRD_OVERVIEW_HTML = """
<h2>Badania XRD — od dyfraktogramu do interpretacji materiału</h2>
<p><b>XRD (X-ray Diffraction)</b> to dyfrakcja rentgenowska: metoda badania
materiałów krystalicznych. Mierzy się natężenie promieniowania ugiętego
w funkcji kąta <b>2θ</b>, otrzymując <b>dyfraktogram</b>. Położenia i kształty
refleksów pomagają identyfikować fazy i analizować strukturę materiału.</p>
<h3>Do czego służy XRD Analyzer 20.3.5?</h3>
<ul>
<li>Wczytuje dyfraktogramy, pokazuje wykresy intensywność–2θ oraz pozwala
porównywać serie pomiarowe.</li>
<li>Wykrywa piki i szacuje ich położenia, intensywności i szerokości;
przetwarza tło i profile dyfrakcyjne.</li>
<li>Wspomaga identyfikację faz z baz wewnętrznych oraz zewnętrznej COD.</li>
<li>Umożliwia analizę Scherrera, Williamsona–Halla i sin²ψ, zgodnie z
założeniami metod i metadanymi aparatury.</li>
<li>Tworzy wykresy i zestawienia do dokumentacji technicznej i artykułów:
dyfraktogramy, oznaczenia refleksów, porównania próbek, profile pików,
wykresy Williamsona–Halla i sin²ψ.</li>
</ul>
<h3>Import danych XRD — rzeczywiście obsługiwane rozszerzenia</h3>
<table border='1' cellpadding='5' cellspacing='0'>
<tr><th>Rozszerzenie</th><th>Rodzaj danych i ograniczenia</th></tr>
<tr><td><b>.udf</b></td><td>Tekstowy format PANalytical UDF z sekcją RawScan.</td></tr>
<tr><td><b>.rd</b></td><td>Binarny Philips/PANalytical, rozpoznawany wariant V3RD; walidacja układu przed odczytem.</td></tr>
<tr><td><b>.xrdml, .xml</b></td><td>Pliki pomiarów XRDML o zgodnej strukturze XML.</td></tr>
<tr><td><b>.cif</b></td><td>Dyfraktogram proszkowy pdCIF zawierający profil intensywności; zwykły CIF strukturalny nie wystarcza.</td></tr>
<tr><td><b>.csv, .txt, .dat, .asc</b></td><td>Tekstowe dwie kolumny: 2θ [°] i intensywność; separator spacja, tabulator, średnik lub przecinek.</td></tr>
</table>
<p><b>Nieobsługiwane bez konwersji:</b> m.in. .raw, .brml, .xy oraz
inne formaty binarne aparatów XRD. Zapis samego rozszerzenia nie gwarantuje
poprawnego importu: liczy się rzeczywista struktura pliku.</p>
<h3>Pozostałe formaty</h3>
<ul>
<li><b>Projekty:</b> .xrdproj i .json (otwieranie/zapis projektu).</li>
<li><b>Baza faz COD:</b> .db3, .sqlite, .sqlite3 i .db (SQLite o zgodnym schemacie).</li>
<li><b>Eksport wykresów:</b> .png, .tif, .tiff, .pdf, .svg i .eps.</li>
<li><b>Eksport tabel:</b> .csv oraz .xlsx.</li>
</ul>
<p><b>Przebieg pracy:</b> import danych → sprawdzenie λ i geometrii pomiaru →
analiza pików i faz → tworzenie wykresów → weryfikacja interpretacji → eksport.</p>
<p><b>Ograniczenia:</b> program ma charakter badawczy; sam wynik automatycznego
dopasowania nie jest potwierdzeniem składu fazowego. Wyniki Scherrera,
Williamsona–Halla i sin²ψ wymagają odpowiedniej korekty instrumentalnej,
założeń modelu oraz oceny przez badacza.</p>
"""

COD_HELP_HTML = """
<h2>Bazy krystalograficzne COD — pobieranie i użycie</h2>
<p><b>Crystallography Open Database (COD)</b> jest otwartą bazą danych
struktur krystalicznych. Udostępnia m.in. parametry komórek elementarnych,
skład, grupy przestrzenne i źródła publikacji.</p>
<p><b>Oficjalna witryna i dane COD:</b><br>
<a href='https://www.crystallography.net/cod/'>https://www.crystallography.net/cod/</a></p>
<p><b>Gotowa baza SQLite (.db3) przygotowywana dla Profex:</b><br>
<a href='https://www.profex-xrd.org/download/'>https://www.profex-xrd.org/download/</a></p>
<p>Przykładowy plik: <b>cod-260101.db3</b> (wydanie datowane
1 stycznia 2026 r.). Wybierz <b>Bazy faz → COD .db3 → Wczytaj bazę COD .db3</b>.
Baza Profex może mieć inny schemat niż obsługiwany przez XRD Analyzer;
jeżeli pojawi się błąd wczytywania, konieczna będzie konwersja.
Nie nadpisuj pliku źródłowego.</p>
<h3>COD Offline w XRD Analyzer 20.3.5</h3>
<p>W zakładce COD można także pobierać dane przez oficjalne API według
pierwiastków lub identyfikatorów, a następnie zbudować własną bazę SQLite.
Wyszukiwanie w pobranej bazie działa później bez internetu.
Opcjonalny import pełnego archiwum CIF wymaga dużo miejsca i czasu.
</p>
<p><b>Różnica:</b> rekordy strukturalne CIF nie są gotowymi eksperymentalnymi
dyfraktogramami — położenia refleksów można oszacować na podstawie
struktury, ale nie oznacza to zgodności intensywności z rzeczywistym pomiarem.</p>
<p>COD jest udostępniana na zasadach <b>CC0</b>. Przy wykorzystaniu wyników
w publikacjach należy podawać oryginalne źródła i identyfikatory COD,
a także cytować właściwe publikacje krystalograficzne.</p>
"""


def show_help_dialog(parent, title: str, html: str):
    """Open scrollable rich-text, clickable external links; no app state changes."""
    from PyQt6.QtWidgets import QDialog, QVBoxLayout, QTextBrowser, QDialogButtonBox
    dlg = QDialog(parent)
    dlg.setWindowTitle(title)
    dlg.resize(860, 680)
    layout = QVBoxLayout(dlg)
    browser = QTextBrowser(dlg)
    browser.setOpenExternalLinks(True)
    browser.setHtml(html)
    layout.addWidget(browser)
    buttons = QDialogButtonBox(QDialogButtonBox.StandardButton.Close, parent=dlg)
    buttons.rejected.connect(dlg.reject)
    layout.addWidget(buttons)
    return dlg.exec()
