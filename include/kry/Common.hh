#ifndef KRY_INPUT_COMMON
#define KRY_INPUT_COMMON

#include <iostream>
#include <vector>
#include <sstream>

namespace kry { namespace input {

  using ulong = unsigned long;

struct SourceLocation
{
    SourceLocation(ulong r, ulong c) : row{r}, col{c} {}
    ulong row{}, col{};
};

std::ostream & operator << (std::ostream &o, const SourceLocation &s);

ulong rowIndex(const std::string &s, ulong r);
ulong nextRow(const std::string &s, ulong pos);
ulong previousRow(const std::string &s, ulong pos);
ulong startOfSection(const std::string &text, ulong start, std::string delim);
ulong endOfSection(const std::string &text, ulong start, std::string delim);
std::vector<ulong> rowIndices(const std::string &text, ulong start, ulong end);
ulong csvRowSize(const std::string &text, ulong rid);

std::vector<std::string> & split(const std::string &s, char delim, 
    std::vector<std::string> &elems);
std::vector<std::string> split(const std::string &s, char delim);

SourceLocation sourceCoordsFromAbsolute(const std::string &text, ulong loc);

}}

#endif
