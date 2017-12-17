# Boost Matheval

[![Build status][travis-svg]][travis-link]
[![AppVeyor build status][appveyor-svg]][appveyor-link]
[![Code coverage report][codecov-svg]][codecov-link]
[Documentation][doxygen-link]

This header-only C++11 libary uses
[Boost.Spirit](http://www.boost.org/libs/spirit/index.html) and
[Boost.Phoenix](http://www.boost.org/libs/phoenix/index.html) to parse
and evaluate mathematical expressions.

### Motivating example 1

```cpp
#include <iostream>
#include "matheval.hpp"

int main()
{
    std::cout << matheval::parse<double>("1+1",{}) << '\n';
}
```
Outputs
```
2
```

### Motivating example 2

```cpp
#include <iostream>
#include <map>
#include "matheval.hpp"

int main()
{
    std::map<std::string,double> symtab;
    symtab.emplace(std::make_pair("x",  2));
    symtab.emplace(std::make_pair("y", -1));
    std::cout << matheval::parse<double>("cbrt(x/2 + sqrt(x**2/4 + y**3/24))",symtab) << '\n';
}
```
Outputs
```
1.25548
```

### Motivating example 3

We can also evaluate an expression multiple times without paying the
cost of parsing again.
```cpp
#include <iostream>
#include <map>
#include "matheval.hpp"

int main()
{
    matheval::Parser<double> parser;
    parser.parse("x + 1");
    std::cout << parser.evaluate({std::make_pair("x",1)}) << ' '
              << parser.evaluate({std::make_pair("x",2)}) << '\n';
}
```
Outputs
```
2 3
```

### Build instructions

Since Boost Matheval is header-only you can simply copy `matheval.hpp`
into your project and start using it.  If you want to build the
examples or run the tests, you can build them using CMake.
```bash
mkdir build
cd build
cmake ..
make         # build the examples
make check   # build and run the tests
```

### Requirements and Limitations

General:

* C++11 compatible compiler, i.e. GCC >= 4.8, Clang, Visual Studio 2015.
* Boost.Spirit, Boost.Phoenix, and Boost.MathConstants.
* No support for ternary functions (e.g. `if`).
* No support for complex numbers.

### Alternatives

* [GNU libmatheval](https://www.gnu.org/software/libmatheval/) is a C
  library built using the parser generator YACC with about the same
  scope as Boost Matheval.  It is a not header-only and does not have
  a C++ interface at the moment, but it is much faster in terms of
  compilation time.
* [ExprTk](http://www.partow.net/programming/exprtk/) is a true
  monster with almost 40000 lines in a single header file.  It
  implements a complete state machine including things like logical
  operations, control structures, and even file IO.  Compilation time
  is even longer than with Boost Matheval.

### License

Distributed under the [Boost Software License, Version 1.0](http://boost.org/LICENSE_1_0.txt).

[travis-svg]: https://travis-ci.org/hmenke/boost_matheval.svg?branch=master
[travis-link]: https://travis-ci.org/hmenke/boost_matheval
[appveyor-svg]: https://ci.appveyor.com/api/projects/status/bphe1739kownt81c/branch/master?svg=true
[appveyor-link]: https://ci.appveyor.com/project/hmenke/boost-matheval/branch/master
[codecov-svg]: https://codecov.io/gh/hmenke/boost_matheval/branch/master/graph/badge.svg
[codecov-link]: https://codecov.io/gh/hmenke/boost_matheval
[doxygen-link]: https://hmenke.github.io/boost_matheval
