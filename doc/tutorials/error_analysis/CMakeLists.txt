configure_tutorial_target(TARGET tutorial_err DEPENDS
                          error_analysis_part1.ipynb error_analysis_part2.ipynb)

nb_export(TARGET tutorial_err SUFFIX "1" FILE "error_analysis_part1.ipynb"
          HTML_RUN)
nb_export(TARGET tutorial_err SUFFIX "2" FILE "error_analysis_part2.ipynb"
          HTML_RUN)
