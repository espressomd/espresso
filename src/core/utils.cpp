#include <cstring>
#include "utils.hpp"

char *strcat_alloc(char *left, const char *right) {
  if (!left) {
    char *res = (char *)Utils::malloc(strlen(right) + 1);
    strcpy(res, right);
    return res;
  } else {
    char *res = Utils::realloc(left, strlen(left) + strlen(right) + 1);
    strcat(res, right);
    return res;
  }
}
