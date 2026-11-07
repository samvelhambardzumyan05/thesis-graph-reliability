"""Generate LaTeX tables from results/results.json -> thesis/tbl_*.tex (so numbers are never typed by hand)."""
import json, os
ROOT = os.path.join(os.path.dirname(__file__), "..")
d = json.load(open(os.path.join(ROOT, "results", "results.json"), encoding="utf-8"))
T = os.path.join(ROOT, "thesis")

def f(x, k=4):
    return f"{x:.{k}f}".replace(".", "{,}")
def w(name, s):
    open(os.path.join(T, name), "w", encoding="utf-8").write(s)

# structure
rows = "".join(f"{r['name']} & {r['n']} & {r['m']} & {r['delta']} & {r['kappa']} & {r['lam']} & {r['bridges']} & {r['cut']}\\\\\n\\hline\n" for r in d["structure"])
w("tbl_structure.tex", "\\begin{tabular}{|l|c|c|c|c|c|c|c|}\n\\hline\nԳրաֆ & $n$ & $m$ & $\\delta$ & $\\kappa$ & $\\lambda$ & կամուրջներ & կտրող կետեր\\\\\n\\hline\n" + rows + "\\end{tabular}\n")

# reliability at selected p
ps = d["ps"]; sel = [0.5, 0.8, 0.9, 0.95]
hdr = " & ".join(f"$p={f(p,2)}$" for p in sel)
rows = "".join(f"{n} & " + " & ".join(f(v[ps.index(p)]) for p in sel) + "\\\\\n\\hline\n" for n, v in d["curves"].items())
w("tbl_curves.tex", "\\begin{tabular}{|l|c|c|c|c|}\n\\hline\nԳրաֆ & " + hdr + "\\\\\n\\hline\n" + rows + "\\end{tabular}\n")

# chords
rows = "".join(f"{n} & {f(v[0])} & {f(v[1])}\\\\\n\\hline\n" for n, v in d["chords"].items())
w("tbl_chords.tex", "\\begin{tabular}{|l|c|c|}\n\\hline\nԳրաֆ & $R$ ($p=0{,}7$) & $R$ ($p=0{,}9$)\\\\\n\\hline\n" + rows + "\\end{tabular}\n")

# bounds
rows = "".join(f"{b['name']} & {f(b['lower'])} & {f(b['exact'])} & {f(b['upper'])}\\\\\n\\hline\n" for b in d["bounds"])
w("tbl_bounds.tex", "\\begin{tabular}{|l|c|c|c|}\n\\hline\nԳրաֆ & ստորին & ճշգրիտ $R$ & վերին\\\\\n\\hline\n" + rows + "\\end{tabular}\n")

# timing
def tm(x): return "--" if x is None else f"{x:.4f}".replace(".", "{,}")
rows = "".join(f"{r['m']} & {tm(r.get('t_enum'))} & {tm(r['t_fact'])} & {tm(r['t_mc'])} & {f(r['R_exact'],5)} & {f(r['R_mc'],5)} & $\\pm${f(r['hw'],4)}\\\\\n\\hline\n" for r in d["timing"])
w("tbl_timing.tex", "\\begin{tabular}{|c|c|c|c|c|c|c|}\n\\hline\n$m$ & թվարկում, վրկ & գործոն., վրկ & MC, վրկ & $R$ ճշգրիտ & $R$ MC & 95\\% միջակայք\\\\\n\\hline\n" + rows + "\\end{tabular}\n")
inside = sum(abs(r['R_mc'] - r['R_exact']) <= r['hw'] for r in d["timing"])
w("stat_inside.tex", f"{inside}\n")

# MC error
rows = "".join(f"{N} & {f(e,5)}\\\\\n\\hline\n" for N, e in zip(d["mc_N"], d["mc_err"]))
w("tbl_mc.tex", "\\begin{tabular}{|r|c|}\n\\hline\n$N$ & միջին $|\\widehat R_N-R|$\\\\\n\\hline\n" + rows + "\\end{tabular}\n")
print("inside CI:", inside, "of", len(d["timing"]))
