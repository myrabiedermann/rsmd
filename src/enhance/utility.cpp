#include "enhance/utility.hpp"




std::string enhance::trimString(const std::string& input)
{
    auto start = input.begin();
    auto end   = input.end();

    // remove whitespaces from beginning
    while( start != end && std::isspace(*start) )
    {
        ++ start;
    }
    // remove whitespaces from end
    while( start != end && std::isspace( *(end-1) ) )
    {
        -- end;
    }

    return std::string(start, end);
}



std::vector<std::string> enhance::splitString(const std::string& input, char delimiter)
{
    std::vector<std::string> output {};
    std::size_t currentPosition = 0;
    std::size_t nextDelimiter = 0;
    std::size_t length = input.length();
    do
    {
        // find next position of delimiter character
        nextDelimiter = input.find(delimiter, currentPosition);
        // save substring
        output.emplace_back( input.substr(currentPosition, nextDelimiter-currentPosition) );
        // update currentPosition
        currentPosition = nextDelimiter < length ? nextDelimiter + 1 : length;
    }while( currentPosition < length || nextDelimiter < length);

    return output;
}



