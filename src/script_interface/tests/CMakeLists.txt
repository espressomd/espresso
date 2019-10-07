include(unit_test)

unit_test(NAME ScriptInterface_test SRC ScriptInterface_test.cpp DEPENDS ScriptInterface)
unit_test(NAME AutoParameters_test SRC AutoParameters_test.cpp DEPENDS ScriptInterface)
unit_test(NAME AutoParameter_test SRC AutoParameter_test.cpp DEPENDS ScriptInterface)
unit_test(NAME Variant_test SRC Variant_test.cpp DEPENDS ScriptInterface)
unit_test(NAME get_value SRC get_value_test.cpp DEPENDS ScriptInterface)
unit_test(NAME None_test SRC None_test.cpp DEPENDS ScriptInterface)
unit_test(NAME LocalContext_test SRC LocalContext_test.cpp DEPENDS ScriptInterface)
