# Third-party components and data

This repository's own code uses external libraries which have independent licenses and are **not** sublicensed under the source-code GPL.

- **PyQt6 / Qt** — GPL/commercial dual licensing; the public edition of this application uses GNU GPL-3.0-only for its original source. Review your distribution obligations, including Qt components.
- **NumPy**, **SciPy**, **Matplotlib**, **openpyxl**, **Gemmi** — each distributed under its own upstream license; check installed-version license metadata before redistributing binary packages.
- **Crystallography Open Database (COD)** — crystallographic data made publicly available under COD's own terms (CC0 data dedication). Original publications should be cited when using structures in research.
- **Profex-provided COD SQLite** — third-party conversion/distribution. Not bundled here; refer to the provider for current terms and schema.
- **Instrument file formats** — may be vendor-specific; interoperability with one version does not establish support for all variants.

This repository contains no PyQt/Qt, COD mirror, instrument-proprietary libraries, or third-party binary distributions.
