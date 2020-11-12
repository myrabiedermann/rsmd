#pragma once

#include <string>
#include <vector>
#include <fstream>
#include <iostream>

#include <random>
#include <algorithm>



// 
// some useful functionality
//

namespace enhance
{
    // remove whitespaces before/after string
    std::string trimString(const std::string&);

    // split string at every occurence of char
    std::vector<std::string> splitString(const std::string&, char);

}


