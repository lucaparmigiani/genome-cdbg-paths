#pragma once
#include <cstdlib>
#include <cstring>
#include <string>
#include <iostream>

struct str {
    char* c_str;
    size_t n;

    str() { n = 0; c_str = 0; }

    str(const char* char_ptr, size_t n) : n(n) {
        c_str = (char*) malloc(n*sizeof(char));
        memcpy(c_str, char_ptr, n);
    }

    str(std::string S) : n(S.length()) {
        c_str = (char*) malloc(n*sizeof(char));
        memcpy(c_str, S.c_str(), n);
    }

    //~str() { // no destructor. ownership is better
    //    if(c_str) free(c_str); 
    //}

    inline size_t length() const { return n; }
    inline char operator[](size_t i) const { return c_str[i]; }

    void clean() {
        n = 0;
        free(c_str);
    }

    std::string get_string() {
        return std::string(reinterpret_cast<char*>(c_str), n);
    }

    void print() {
        for (int i = 0; i < n; i++) {
            std::cout << c_str[i];
        }
        std::cout << std::flush;
    }
    void println() {
        print();
        std::cout << '\n' << std::flush;
    }

};
