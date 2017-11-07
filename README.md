# Boost Matheval

[![Build status][travis-svg]][travis-link]
[![AppVeyor build status][appveyor-svg]][appveyor-link]

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

* C++11 compatible compiler, i.e. GCC >= 5.0, Clang, Visual Studio 2015 (maybe).
* Boost.Spirit, Boost.Phoenix, and Boost.MathConstants (versions in Ubuntu 16.04 are known to work).
* No support for ternary functions (e.g. `if`).
* No support for complex numbers.

### License

Distributed under the [Boost Software License, Version 1.0](http://boost.org/LICENSE_1_0.txt).

[travis-svg]: https://travis-ci.org/hmenke/boost_matheval.svg?branch=master
[travis-link]: https://travis-ci.org/hmenke/boost_matheval
[appveyor-svg]: https://ci.appveyor.com/api/projects/status/bphe1739kownt81c/branch/master?svg=true
[appveyor-link]: https://ci.appveyor.com/project/hmenke/boost-matheval/branch/master
