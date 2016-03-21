check_PROGRAMS = Vector_test Variant_test Factory_test ObjectContainer_test MpiCallbacks_test

AM_CPPFLAGS=-I$(srcdir)/..
AM_LDFLAGS=$(BOOST_UNIT_TEST_FRAMEWORK_LIBS)

Vector_test_SOURCES = Vector_test.cpp
Variant_test_SOURCES = Variant_test.cpp
Factory_test_SOURCES = Factory_test.cpp testing_utils/HelloWorld.cpp testing_utils/HelloWorld.hpp
ObjectContainer_test_SOURCES = ObjectContainer_test.cpp ../ObjectContainer.hpp
MpiCallbacks_test_SOURCES = MpiCallbacks_test.cpp ../MpiCallbacks.hpp ../MpiCallbacks.cpp

TESTS = $(check_PROGRAMS)
