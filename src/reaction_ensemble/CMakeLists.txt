add_library(ReactionEnsemble SHARED "reaction_ensemble.cpp")
target_link_libraries(ReactionEnsemble PUBLIC EspressoCore)
target_include_directories(ReactionEnsemble PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
