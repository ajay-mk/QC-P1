#include "general.h"
#include <fstream>

using std::cout;
using std::endl;

int count_lines(std::string filename)
{
    int nlines;
    std::string line;
    std::ifstream input (filename);
    while(!input.eof())
        ++ nlines;
    return nlines;
}