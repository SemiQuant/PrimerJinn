# Fonts and logo

Copied from the DrDx.Me offline tools (`drdx_me_combined/assets/fonts` and `assets/img`, 2026-09-26) so the GUI looks
like the other DrDx.Me dashboard tools and works offline: `gui.py` inlines them into the page as `data:` URIs.

| File | Family | Weight | Source | Licence |
|---|---|---|---|---|
| montserrat-latin-600-normal.woff2 | Montserrat | 600 | npm `@fontsource/montserrat` 5.3.0 | SIL OFL 1.1 |
| montserrat-latin-700-normal.woff2 | Montserrat | 700 | npm `@fontsource/montserrat` 5.3.0 | SIL OFL 1.1 |
| roboto-latin-400-normal.woff2 | Roboto | 400 | npm `@fontsource/roboto` 5.3.0 | SIL OFL 1.1 |
| roboto-latin-500-normal.woff2 | Roboto | 500 | npm `@fontsource/roboto` 5.3.0 | SIL OFL 1.1 |
| roboto-mono-latin-400-normal.woff2 | Roboto Mono | 400 | npm `@fontsource/roboto-mono` 5.3.0 | SIL OFL 1.1 |

Latin subset only. The licence text is in `OFL.txt`. `DrDxMe_white.png` is the DrDx.Me wordmark for the header.
`style.css` is the dashboard theme (`sq-dashboard.css` tokens) mapped onto Streamlit's DOM.
