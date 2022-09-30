configure_tutorial_target(TARGET tutorial_widom DEPENDS widom_insertion.ipynb)

nb_export(TARGET tutorial_widom SUFFIX "" FILE "widom_insertion.ipynb" HTML_RUN
          VAR_SUBST "\"ci_params={'mesh':(12,12,12),'cao':6}\"")
