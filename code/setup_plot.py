import matplotlib.pyplot as plt
import matplotlib as mpl

def on_press(event):
    """Print the mouse-button and data coordinates of a click event."""
    print("my position:", event.button, event.xdata, event.ydata)


def setup_plot(fontsize=8, fonttype="sans-serif", dpi=150):
    mpl.rcParams["font.size"] = fontsize
    mpl.rcParams["figure.dpi"] = dpi
    plt.rcParams["pdf.fonttype"] = 42

    if fonttype == "serif":
        mpl.rc("text", usetex=True)
        mpl.rc(
            "text.latex",
            preamble=r"\usepackage{amsmath, newtxtext, newtxmath}",
        )
        plt.rcParams["font.family"] = "serif"
    elif fonttype == "sans-serif":
        mpl.rc("text", usetex=True)
        plt.rcParams["font.family"] = "sans-serif"
        latex_preamble = r"""
        \usepackage[T1]{fontenc}
        \usepackage{bm, amsmath, sansmathfonts}

        \makeatletter
        \AtBeginDocument{
            \DeclareSymbolFont{sansextrabold}{T1}{cmss}{bx}{n}
            \DeclareMathSymbol{+}{\mathbin}{sansextrabold}{"2B}
            \DeclareMathSymbol{=}{\mathrel}{sansextrabold}{"3D}
            \DeclareMathSymbol{<}{\mathrel}{sansextrabold}{"3C}
            \DeclareMathSymbol{>}{\mathrel}{sansextrabold}{"3E}
            \DeclareMathSymbol{|}{\mathord}{sansextrabold}{"7C}
            \DeclareMathSymbol{/}{\mathord}{sansextrabold}{"2F}

            \DeclareSymbolFont{boldsanssymbols}{OMS}{cmsssy}{b}{n}
            \DeclareMathSymbol{-}{\mathbin}{boldsanssymbols}{"00}
            \DeclareMathSymbol{\times}{\mathbin}{boldsanssymbols}{"02}

            \DeclareMathSymbol{\oplus}{\mathbin}{boldsanssymbols}{"08}
            \DeclareMathSymbol{\otimes}{\mathbin}{boldsanssymbols}{"0A}
            \DeclareMathSymbol{\approx}{\mathrel}{boldsanssymbols}{"19}
            \DeclareMathSymbol{\perp}{\mathrel}{boldsanssymbols}{"3F}
            \DeclareMathSymbol{\sim}{\mathrel}{boldsanssymbols}{"18}
            \DeclareMathSymbol{\to}{\mathrel}{boldsanssymbols}{"21}
            \DeclareMathSymbol{\leftarrow}{\mathrel}{boldsanssymbols}{"20}
            \DeclareMathSymbol{\in}{\mathrel}{boldsanssymbols}{"32}

            \let\original@le\le
            \let\original@ge\ge
            \let\original@langle\langle
            \let\original@rangle\rangle

            \renewcommand{\le}{\bm{\original@le}}
            \renewcommand{\ge}{\bm{\original@ge}}
            \renewcommand{\langle}{\bm{\original@langle}}
            \renewcommand{\rangle}{\bm{\original@rangle}}
            \renewcommand{\|}{|\hspace{-1pt}|}
        }
        \makeatother

        \newcommand{\p}{\partial}
        \newcommand{\T}{^{\mathrm{T}}}
        \renewcommand{\vec}[1]{\boldsymbol{#1}}
        \newcommand{\bn}{\vec{\nabla}}
        \newcommand\ii{\mathrm{i}}
        \newcommand\ee{\mathrm{e}}
        \newcommand{\widebar}[1]{\mskip.5\thinmuskip\overline{\mskip-.5\thinmuskip {#1} \mskip-.5\thinmuskip}\mskip.5\thinmuskip}
        \newcommand{\ket}[1]{| #1 \rangle}
        \newcommand{\bra}[1]{\langle #1 |}
        """
        mpl.rc("text.latex", preamble=latex_preamble)
    else:
        raise ValueError("fonttype must be 'serif' or 'sans-serif'.")

