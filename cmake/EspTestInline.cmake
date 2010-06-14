# - Define macro to check inline keyword
#
#  GMX_TEST_INLINE(VARIABLE)
#
#  VARIABLE will be set to the keyword
#
#  Remember to have a cmakedefine for it too...

MACRO(ESP_TEST_INLINE VARIABLE)
    IF(NOT DEFINED TEST_${VARIABLE})

        MESSAGE(STATUS "Checking for inline keyword")

	FOREACH(KEYWORD "static inline" "inline static" "inline" "static")
            IF(NOT TEST_${VARIABLE})
                MESSAGE(STATUS "Checking for inline keyword - try ${KEYWORD}")
                TRY_COMPILE(TEST_${VARIABLE} "${CMAKE_BINARY_DIR}"    
                            "${CMAKE_SOURCE_DIR}/cmake/TestInline.c"
                            COMPILE_DEFINITIONS -D"INLINEDEF=${KEYWORD}" )
                SET(CHK_INLINE_KEYWORD ${KEYWORD})
            ENDIF(NOT TEST_${VARIABLE})
        ENDFOREACH(KEYWORD)
             
        IF(TEST_${VARIABLE})
            SET(${VARIABLE} ${CHK_INLINE_KEYWORD})
            MESSAGE(STATUS "Checking for inline keyword - using ${CHK_INLINE_KEYWORD}")
        ELSE(TEST_${VARIABLE})
	    SET(${VARIABLE} " ")
            MESSAGE(FATAL_ERROR "Checking for inline keyword - none found")
        ENDIF(TEST_${VARIABLE})

    ENDIF(NOT DEFINED TEST_${VARIABLE})        
ENDMACRO(ESP_TEST_INLINE VARIABLE)




