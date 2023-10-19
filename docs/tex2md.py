# tex2md.py
"""Translate LaTex style to Markdown style:
A tiny version of pandoc for sphinx with MyST.
"""

import re


FLAGS = re.DOTALL | re.MULTILINE


# Utilities -------------------------------------------------------------------
def _header(text):
    text = re.sub(r"^.+\\begin\{document\}\n", '', text, flags=FLAGS)
    return re.sub(r"\\end\{document\}", '', text, flags=FLAGS)


def _bold(text):
    return re.sub(r"\\textbf\{(.+?)\}", r"**\1**", text, flags=FLAGS)


def _italics(text):
    return re.sub(r"\\emph\{(.+?)\}", r"_\1_", text, flags=FLAGS)


def _texttt(text):
    return re.sub(
        r"\\texttt\{(.+?)\}", r"`\1`",
        re.sub(r"\\_", r"_", text, flags=FLAGS),
        flags=FLAGS
    )


def _eqref(text):
    return re.sub(r"\\eqref\{(.+?)\}", r"{eq}`\1`", text, flags=FLAGS)


def _refs(text):
    for prefix in ["Table", "Figure"]:
        text = re.sub(
            prefix + r"~\\ref\{(.+?)\}",
            f"[{prefix} TODO](\\1)",
            text,
            FLAGS,
        )
    return re.sub(r"Algorithm~\\ref\{(.+?)\}", r"{prf:ref}`\1`",
                  text, flags=FLAGS)


def _section(text):
    text = re.sub(r"\\section[\*]*\{(.+?)\}",
                  r"## \1", text, flags=FLAGS)
    text = re.sub(r"\\subsection[\*]*\{(.+?)\}",
                  r"### \1", text, flags=FLAGS)
    return re.sub(r"\\subsubsection[\*]*\{(.+?)\}",
                  r"#### \1", text, flags=FLAGS)


def _simbackslash(text):
    return re.sub(r"(\w)~\\", r"\1 \\", text, flags=FLAGS)


def _trp(text):
    text = re.sub(r"\\trp", r"^{\\mathsf{T}}", text, flags=FLAGS)
    return re.sub(r"\\invtrp", r"^{-\\mathsf{T}}", text, flags=FLAGS)


def _bbR(text):
    return re.sub(r"\\R", r"\\mathbb{R}", text, flags=FLAGS)


def _align(text):
    # Numberless equations
    text = re.sub(r"^\s*?\\begin\{align\*\}",
                  r"\n$$\n\\begin{aligned}", text, flags=FLAGS)
    text = re.sub(r"^\s*?\\end\{align\*\}",
                  r"\\end{aligned}\n$$\n", text, flags=FLAGS)
    # Numbered equations
    return re.sub(
        r"^\s*?\\begin\{align\}\s*?\\label\{(.+?)\}(.*?)^\s*?\\end\{align\}",
        r"\n$$\n\\begin{aligned}\2\n\\end{aligned}\n$$ (\1)\n",
        text, flags=FLAGS,
    )


def _listing(text):
    return re.sub(r"\\begin\{lstlisting\}\[.+?\](.+?)\\end\{lstlisting\}",
                  r"```matlab\1```", text, flags=FLAGS)


# Main routine ----------------------------------------------------------------
def main(infile: str, outfile: str = None) -> None:
    r"""Read a file, replace LaTex elements with Markdown elements, and export.

    Substitutions:

        \textbf{text}           ->      **text**
        \emph{text}             ->      _text_
        \texttt{a\_text}        ->      `a_text`
        \eqref{mylabel}         ->      {eq}`mylabel`
        Algorithm~\ref{mylabel} ->      {prf:ref}`mylabel`
        Table~\ref{mylabel}     ->      [Table X](mylabel)
        Figure~\ref{mylabel}    ->      [Figure X](mylabel)
        \section{TheSection}    ->      ## TheSection
        \subsection{TheSection} ->      ### TheSection
        ~\                      ->      \
        \trp                    ->      ^{\mathsf{T}}
        \invtrp                 ->      ^{-\mathsf{T}}
        \R                      ->      \mathbb{R}

        \begin{align}                   $$
            \label{mylabel}             \begin{aligned}
            SOME MATH           ->          SOME MATH
        \end{align}                     \end{aligned}
                                        $$ (mylabel)

        \begin{align*}                  $$
            SOME MATH                   \begin{aligned}
        \end{align*}           ->           SOME MATH
                                        \end{aligned}
                                        $$

        \begin{lstlisting}              ```matlab
            SOME CODE           ->      SOME CODE
        \end{listlisting}               ```

    Parameters
    ----------
    infile : str
        File to read from.
    outfile : str or None
        File to write to. If None (default), print the results.
    """
    # Read the input file.
    with open(infile, "r") as f:
        contents = f.read().strip()

    # Make the substitutions.
    for formatter in (
        _header,
        _bold,
        _italics,
        _texttt,
        _eqref,
        _refs,
        _section,
        _simbackslash,
        _trp,
        _bbR,
        _align,
        _listing,
    ):
        contents = formatter(contents)

    # If no outfile is provided, print the modified contents and quit.
    if outfile is None:
        return print(contents)

    # If an outfile is provided, write the modified contents to that file.
    with open(outfile, "w") as f:
        f.write(f"# {os.path.basename(outfile)}\n\n{contents}")


def _tests():
    if (out := _bold(r"\textbf{thetext}")) != "**thetext**":
        print(f"_bold() failed: {out}")
    if (out := _italics(r"\emph{thetext}")) != "_thetext_":
        print(f"_italics() failed: {out}")
    if (out := _texttt(r"\texttt{the\_text}")) != "`the_text`":
        print(f"_texttt() failed: {out}")
    if (out := _eqref(r"\eqref{myref}")) != "{eq}`myref`":
        print(f"_eqref() failed: {out}")
    if (out := _refs(r"Algorithm~\ref{myref}")) != "[Algorithm TODO](myref)":
        print(f"_refs() failed: {out}")
    if (out := _simbackslash(r"Equation~\ref{4}")) != r"Equation \ref{4}":
        print(f"_simbackslash() failed: {out}")
    out = _align(r"""
        \begin{align}
            \label{myequation}
            Here is some labeled math
            ...and some more
            And even more!
        \end{align}

        \begin{align*}
            Here is some unlabeled math.
            And some more!
        \end{align*}
    """)
    print(out)
    out = _listing(r"""
\begin{lstlisting}[style=Matlab-editor,frame=single,numbers=left]
classdef MyConstraint < Constraint

    methods (Access = public)

        function [u] = State_Solve(this, z)
            error('StateSolve() not implemented');
        end

        function [Mv] = c_u_Transpose_Inverse_Apply(this, v, u, z)
            error('c_u_Transpose_Inverse_Apply() not implemented');
        end

        function [Mv] = c_z_Transpose_Apply(this, v, u, z)
            error('c_z_Transpose_Apply() not implemented');
        end

        function [Mv] = c_u_Inverse_Apply(this, v, u, z)
           error('c_u_Inverse_Apply() not implemented');
        end

        function [Mv] = c_z_Apply(this, v, u, z)
            error('c_z_Apply() not implemented');
        end

        function [Mv] = c_uu_Apply(this, v, u, z, lambda)
            error('c_uu_Apply() not implemented');
        end

        function [Mv] = c_uz_Apply(this, v, u, z, lambda)
            error('c_uz_Apply() not implemented');
        end

        function [Mv] = c_zu_Apply(this, v, u, z, lambda)
            error('c_zu_Apply() not implemented');
        end

        function [Mv] = c_zz_Apply(this, v, u, z, lambda)
            error('c_zz_Apply() not implemented');
        end
    end
end
\end{lstlisting}
    """)
    print(out)


# =============================================================================
if __name__ == "__main__":
    import os
    import argparse

    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    fname = os.path.basename(__file__)
    parser.usage = f""" python3 {fname} --help
        python3 {fname} INFILE [-o OUTFILE]
    """

    parser.add_argument("infile", type=str, help="file to operator on")
    parser.add_argument("-o", type=str, default=None, help="file to write to")
    parser.add_argument("--test", action="store_true",
                        help="run internal unit tests")

    args = parser.parse_args()
    if args.test:
        _tests()
    else:
        main(args.infile, args.o)
