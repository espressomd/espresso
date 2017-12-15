/** A little desktop calculator based on Boost Matheval
 *
 * This little calculator can be used for testing Boost Matheval or
 * just as a neat little calculator on the command line.  It supports
 * all the functions and constants from Boost Matheval but defines no
 * symbol table or variables.  GNU readline is used to keep a history of the inputs.
 *
 * c++ -Wall -Wextra -Wpedantic -std=c++11 -I/path/to/boost_matheval/include calc.cpp -lreadline
 */
#include <cstdlib>
#include <cstring>
#include <limits>
#include <iostream>

#include <boost/lexical_cast.hpp>
#include <boost/optional.hpp>

#include <readline/readline.h>
#include <readline/history.h>

#include <matheval.hpp>

template < typename Iterator >
boost::optional<double> evaluate(Iterator first, Iterator last)
{
    try {
        auto res = matheval::parse<double>(first, last, {});
        std::cout << res << '\n';
        return res;
    } catch (std::exception const &e) {
        std::cerr << "Error: " << e.what() << '\n';
        return boost::none;
    }
}

void commandline(int argc, char *argv[])
{
    std::string expr;
    for (int i = 1; i < argc; ++i)
        expr += argv[i];

    evaluate(expr.c_str(), expr.c_str() + expr.size());
}

void interactive()
{
    std::cout << "Type [q or Q] to quit\n";

    constexpr char const prompt[] = "> ";
    char * input = readline(prompt);
    double res = 0;

    while (input)
    {
        std::size_t len = std::strlen(input);

        if (len != 0)
        {
            if (input[0] == 'q' || input[0] == 'Q')
                break;

            add_history(input);

            // Replace % by the last result
            std::string line{input, input + len};
            std::size_t pos = line.find("%");
            if (pos != std::string::npos) {
                line.replace(pos, 1, boost::lexical_cast<std::string>(res));
            }

            // Update result on success
            res = evaluate(line.begin(), line.end()).value_or(res);
        }

        std::free(input);
        input = readline(prompt);
    }

    if (!input) // nullptr signals EOF
        std::cout << '\n';

    std::free(input);
}

int main(int argc, char *argv[])
{
    std::cout.precision(std::numeric_limits<double>::digits10);
    if (argc > 1)
        commandline(argc,argv);
    else
        interactive();
}
