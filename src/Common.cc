#include "kry/Common.hh"

using namespace kry::input;
using std::ostream;
using std::endl;

template <class T>
using vector = std::vector<T>;

using std::string;

ulong kry::input::rowIndex(const string &s, ulong r)
{
    ulong i{0};
    for(ulong rc=0; rc<r; ++rc )
    {
        while(s[i] != '\n') { ++i; } ++i;
    }
    return i;
}

ulong kry::input::nextRow(const string &s, ulong pos)
{
    while(s[pos] != '\n') { ++pos; } ++pos;
    return pos;
}

ulong kry::input::previousRow(const string &s, ulong pos)
{
    for(int i=0; i<2; ++i)
    {
        while(s[pos--] != '\n') {  }
    }
    return ++pos;
}

ulong kry::input::startOfSection(const string &text, ulong start, string delim)
{
    ulong result = text.find(delim, start);
    return nextRow(text, result);
}

ulong kry::input::endOfSection(const string &text, ulong start, string delim)
{
    ulong result = text.find(delim, start);
    return previousRow(text, result);
}

vector<ulong> kry::input::rowIndices(const string &text, ulong start, ulong end)
{
    vector<ulong> result;
    result.push_back(start);
    for(ulong i=start; i<=end; ++i)
    {
        if(text[i] == '\n') { result.push_back(i+1); }
    }
    return result;
}

ulong kry::input::csvRowSize(const string &text, ulong rid)
{
    ulong res{0}, i{rid};
    while(text[i] != '\n'){ if(text[i++] == ',') { ++res; } }
    return res;
}

vector<string> & kry::input::split(const string &s, char delim,
                                   vector<string> &elems) {
    std::stringstream ss(s);
    string item;
    while (std::getline(ss, item, delim)) {
        elems.push_back(item);
    }
    return elems;
}

ostream & kry::input::operator << (ostream &o, const SourceLocation &s)
{
    o << "[" << s.row << ":" << s.col << "]" << endl;
    return o;
}

vector<string> kry::input::split(const string &s, char delim) {
    vector<string> elems;
    split(s, delim, elems);
    return elems;
}

SourceLocation kry::input::sourceCoordsFromAbsolute(
    const string &text, ulong loc)
{
    ulong r{0};
    for(ulong i=0; i<loc; ++i) { if(text[i] == '\n') { ++r; } }
    ulong c = loc - previousRow(text, loc);
    return SourceLocation{r, c};
}
