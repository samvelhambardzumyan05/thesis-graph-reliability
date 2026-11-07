# Գրաֆների կապակցվածության և ցանցի հուսալիության վերլուծություն

Դիպլոմային աշխատանք, Ս. Համբարձումյան. ղեկավար՝ Ն. Բադալյան.

## Կառուցվածք
- `thesis/` LaTeX աղբյուրներ (XeLaTeX; Overleaf-ում Menu → Compiler → XeLaTeX, main file՝ `thesis/main.tex`)
- `src/connectivity.py` BFS/DFS, Տարյանի ալգորիթմ, առավելագույն հոսք, λ և κ
- `src/reliability.py` լրիվ թվարկում, գործոնավորում, Monte Carlo
- `src/experiments.py` փորձեր → `results/results.json` + `thesis/figures/`
- `src/make_tables.py` գեներացնում է `thesis/tbl_*.tex` աղյուսակները results.json-ից
- `tests/` թեստեր (համեմատություն NetworkX-ի հետ)
- `slides/` ներկայացում

## Գործարկում
```
pip install networkx matplotlib pytest
python -m pytest -q tests
python src/experiments.py && python src/make_tables.py
```
Նշում. ժամանակների չափումները կախված են մեքենայից և փոքր-ինչ տարբերվում են գործարկումից գործարկում:
