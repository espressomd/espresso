configure_tutorial_target(
  TARGET tutorial_ek DEPENDS electrokinetics.ipynb figures/profiles.png
  figures/schlitzpore_3d.png scripts/eof_analytical.py scripts/plot.py)

nb_export(TARGET tutorial_ek SUFFIX "" FILE "electrokinetics.ipynb" HTML_RUN
          VAR_SUBST "integration_length=600;dt=0.5" ADD_SCRIPTS
          "scripts/eof_analytical.py;scripts/plot.py")
