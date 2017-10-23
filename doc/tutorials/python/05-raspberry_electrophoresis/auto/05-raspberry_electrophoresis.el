(TeX-add-style-hook
 "05-raspberry_electrophoresis"
 (lambda ()
   (TeX-add-to-alist 'LaTeX-provided-class-options
                     '(("scrartcl" "paper=a4" "fontsize=11pt" "twoside" "footsepline" "headsepline" "headinclude=false" "footinclude=false" "pagesize" "")))
   (TeX-run-style-hooks
    "latex2e"
    "common"
    "scrartcl"
    "scrartcl10"
    "siunitx")
   (LaTeX-add-labels
    "fig:rasp_snapshot"
    "sec:compiler_flags"
    "sec:espresso"
    "effectiveGammaEq")
   (LaTeX-add-environments
    "task"))
 :latex)

