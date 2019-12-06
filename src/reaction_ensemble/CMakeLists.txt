add_library(ReactionEnsemble SHARED "reaction_ensemble.cpp")
target_link_libraries(ReactionEnsemble PRIVATE EspressoCore PUBLIC utils)
target_include_directories(ReactionEnsemble PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
