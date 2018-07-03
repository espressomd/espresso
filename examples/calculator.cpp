#include <iostream>
#include <string>

#include "matheval.hpp"

int main()
{
    std::cout << "Type [q or Q] to quit\n";

    std::string str;
    while (std::getline(std::cin, str).good())
    {
        if (str.empty() || str[0] == 'q' || str[0] == 'Q') {
            break;
        }

        try {
            auto res = matheval::parse<double>(str, {});
            std::cout << res << '\n';
        } catch (std::exception const &e) {
            std::cout << "Error: " << e.what() << '\n';
        }
    }

    std::cout << "Bye... :-)\n";
}
