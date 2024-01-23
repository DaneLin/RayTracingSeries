#pragma once

#include <cstring>

class Hello
{
    public:
    Hello(const std::string& word)
    {
        std::cout <<"Hello " << word << std::endl; 
    }
};