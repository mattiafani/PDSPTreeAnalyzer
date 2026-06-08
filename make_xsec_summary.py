#!/usr/bin/env python3
"""Assemble the plots_xsec_* PDFs into one organized summary PDF.

Scans the plots_xsec_<channel>/ directories produced by plot_xsec.C and lays the
plots out in a single PDF via LaTeX (vector quality preserved).

Two layouts:
  --by type     (default) one section per observable, channels side by side in a
                grid. Best for comparing channels (mirrors the talk).
  --by channel  one section per channel, all of its observables together.

Usage:
  python3 make_xsec_summary.py                 # from the dir holding plots_xsec_*
  python3 make_xsec_summary.py --base . --by channel
  python3 make_xsec_summary.py --no-compile    # just write the .tex

Requires pdflatex/latexmk on PATH for the compile step.
"""
import argparse
import datetime
import os
import subprocess
import sys

# Display order + human names. Unknown keys fall back to the raw name.
CHANNEL_ORDER = ["inel", "abslike", "cexlike", "otherlike", "common"]
CHANNEL_NAME = {
    "inel": "Total inelastic",
    "abslike": "Absorption",
    "cexlike": "Charge exchange",
    "otherlike": "Quasi-elastic + other",
    "common": "Common (signal-independent)",
}
OBS_ORDER = [
    "xsec_total_inelastic",
    "xsec_Ninc_Nint_reco",
    "xsec_purity",
    "xsec_migration_matrix",
    "xsec_KE_reco_vs_true",
    "xsec_efficiency",
]
OBS_NAME = {
    "xsec_total_inelastic": "Cross section vs KE",
    "xsec_Ninc_Nint_reco": "Incident / interacting counts",
    "xsec_purity": "Purity",
    "xsec_migration_matrix": "Migration matrix (reco vs true slice)",
    "xsec_KE_reco_vs_true": "Reco vs true KE per slice",
    "xsec_efficiency": "Selection efficiency",
}
# Optional one-line explanatory note printed under a section header.
OBS_NOTE = {
    "xsec_total_inelastic": ("Red line = Geant4 Bertini input (the MC generator); "
                             "black = MC truth (thin-slice closure, must reproduce the line); "
                             "green open = MC reco, uncorrected."),
}


def tex_escape(s):
    for a, b in [("\\", r"\textbackslash{}"), ("_", r"\_"), ("&", r"\&"),
                 ("%", r"\%"), ("#", r"\#"), ("$", r"\$"), ("{", r"\{"), ("}", r"\}")]:
        s = s.replace(a, b)
    return s


def ascii_only(s):
    return s.encode("ascii", "ignore").decode("ascii")


def discover(base):
    """channel -> {obs_key -> relpath_from_base}."""
    out = {}
    for d in sorted(os.listdir(base)):
        full = os.path.join(base, d)
        if not (os.path.isdir(full) and d.startswith("plots_xsec_")):
            continue
        chan = d[len("plots_xsec_"):]
        files = {f[:-4]: os.path.join(d, f)
                 for f in sorted(os.listdir(full)) if f.lower().endswith(".pdf")}
        if files:
            out[chan] = files
    return out


def ordered(keys, order):
    known = [k for k in order if k in keys]
    rest = sorted(k for k in keys if k not in order)
    return known + rest


def grid(items, ncol=2):
    """items: list of (label, relpath) -> LaTeX minipage grid."""
    w = 0.98 / ncol
    lines = []
    for i, (label, rel) in enumerate(items):
        lines += [
            r"\begin{minipage}[t]{%.3f\textwidth}\centering" % w,
            r"{\small\bfseries %s}\par\vspace{2pt}" % label,
            r"\includegraphics[width=\linewidth]{%s}" % rel,
            r"\end{minipage}",
        ]
        lines.append(r"\par\vspace{10pt}" if (i + 1) % ncol == 0 else r"\hfill")
    return "\n".join(lines)


def name_of(key, names):
    return names.get(key, tex_escape(key))


def build_tex(data, by, base):
    body = []
    if by == "type":
        all_obs = ordered({o for f in data.values() for o in f}, OBS_ORDER)
        for obs in all_obs:
            items = [(name_of(ch, CHANNEL_NAME), data[ch][obs])
                     for ch in ordered(data.keys(), CHANNEL_ORDER) if obs in data[ch]]
            if not items:
                continue
            body.append(r"\section*{%s}" % name_of(obs, OBS_NAME))
            if obs in OBS_NOTE:
                body.append(r"{\small\itshape %s\par}\vspace{6pt}" % OBS_NOTE[obs])
            body.append(grid(items, ncol=2))
            body.append(r"\clearpage")
    else:  # by channel
        for ch in ordered(data.keys(), CHANNEL_ORDER):
            items = [(name_of(o, OBS_NAME), data[ch][o])
                     for o in ordered(data[ch].keys(), OBS_ORDER)]
            if not items:
                continue
            body.append(r"\section*{%s}" % name_of(ch, CHANNEL_NAME))
            body.append(grid(items, ncol=2))
            body.append(r"\clearpage")

    stamp = datetime.datetime.now().strftime("%Y-%m-%d %H:%M")
    base_disp = tex_escape(ascii_only(os.path.abspath(base))) or "."
    head = r"""\documentclass[10pt]{article}
\usepackage[a4paper,margin=1.4cm]{geometry}
\usepackage{graphicx}
\usepackage{float}
\setlength{\parindent}{0pt}
\begin{document}
\begin{center}
{\LARGE\bfseries 0.5 GeV/c $\pi^{+}$--Ar cross-section plots}\\[4pt]
{\large MC-only summary}\\[6pt]
{\small Generated %s \quad layout: by %s}\\
{\footnotesize\ttfamily %s}
\end{center}
\vspace{6pt}
""" % (stamp, by, base_disp)
    return head + "\n".join(body) + "\n\\end{document}\n"


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--base", default=".", help="dir containing plots_xsec_* (default: .)")
    ap.add_argument("--by", choices=["type", "channel"], default="type",
                    help="group by observable (default) or by channel")
    ap.add_argument("--out", default="xsec_summary.pdf", help="output PDF name")
    ap.add_argument("--no-compile", action="store_true", help="write .tex only")
    args = ap.parse_args()

    base = os.path.abspath(args.base)
    data = discover(base)
    if not data:
        sys.exit("No plots_xsec_* directories with PDFs found under %s" % base)
    print("Found channels: " + ", ".join(ordered(data.keys(), CHANNEL_ORDER)))

    stem = os.path.splitext(os.path.basename(args.out))[0]
    tex_path = os.path.join(base, stem + ".tex")
    with open(tex_path, "w") as fh:
        fh.write(build_tex(data, args.by, base))
    print("Wrote " + tex_path)

    if args.no_compile:
        return
    tool = next((t for t in ("latexmk", "pdflatex") if which(t)), None)
    if not tool:
        print("latexmk/pdflatex not found; .tex written, compile it yourself.")
        return
    cmd = (["latexmk", "-pdf", "-silent", "-interaction=nonstopmode", stem + ".tex"]
           if tool == "latexmk" else
           ["pdflatex", "-interaction=nonstopmode", stem + ".tex"])
    runs = 1 if tool == "latexmk" else 2
    for _ in range(runs):
        r = subprocess.run(cmd, cwd=base, stdout=subprocess.DEVNULL, stderr=subprocess.STDOUT)
    if which("latexmk"):
        subprocess.run(["latexmk", "-c", stem + ".tex"], cwd=base,
                       stdout=subprocess.DEVNULL, stderr=subprocess.STDOUT)
    pdf = os.path.join(base, stem + ".pdf")
    print(("Built " + pdf) if os.path.exists(pdf) else
          "Compile failed; inspect %s.log" % os.path.join(base, stem))


def which(prog):
    return any(os.access(os.path.join(p, prog), os.X_OK)
               for p in os.environ.get("PATH", "").split(os.pathsep) if p)


if __name__ == "__main__":
    main()
