(TeX-add-style-hook
 "05-raspberry_electrophoresis"
 (lambda ()
   (TeX-add-to-alist 'LaTeX-provided-class-options
                     '(("scrartcl" "paper=a4" "fontsize=11pt" "twoside" "footsepline" "headsepline" "headinclude=false" "footinclude=false" "pagesize" "")))
   (add-to-list 'LaTeX-verbatim-environments-local "lstlisting")
   (add-to-list 'LaTeX-verbatim-macros-with-braces-local "lstinline")
   (add-to-list 'LaTeX-verbatim-macros-with-delims-local "lstinline")
   (TeX-run-style-hooks
    "latex2e"
    "common"
    "scrartcl"
    "scrartcl10"
    "graphicx"
    "siunitx"
    "verbatim"
    "listings")
   (TeX-add-symbols
    '("EScmd" 1))
   (LaTeX-add-labels
    "fig:rasp_snapshot"
    "sec:compiler_flags"
    "sec:espresso"
    "effectiveGammaEq")
   (LaTeX-add-environments
    "task"))
 :latex)

