#define BOOST_TEST_MODULE variables
#include <boost/test/included/unit_test.hpp>
#include "exprtest.hpp"

double const x = 1;
double const y = 2;

static std::map<std::string, double> symtab() {
    std::map<std::string, double> symtab;
    symtab.insert(std::make_pair("x",  x));
    symtab.insert(std::make_pair("y",  y));
    symtab.insert(std::make_pair("e", -1));
    return symtab;
}

SYMEXPRTEST(var1, "x+y" , symtab(), x+y)
SYMEXPRTEST(var2, "x-y" , symtab(), x-y)
SYMEXPRTEST(var3, "x*y" , symtab(), x*y)
SYMEXPRTEST(var4, "x/y" , symtab(), x/y);
SYMEXPRTEST(var5, "x%y" , symtab(), std::fmod(x,y));
SYMEXPRTEST(var6, "x**y", symtab(), std::pow(x,y));

// Constants have higher priority than variables of the same name
SYMEXPRTEST(var7, "e", symtab(), boost::math::constants::e<double>())
