//
// test the canonical code and decode calls work as expected
//
// g++ -O3 test_process_calls.cc popseq_utils.cc -std=c++11 -static -o test_process_calls
// process_calls input_file[.lzo|.gz|.bz] > output_file

#include <istream>
#include <iostream>
#include <fstream>
#include <cstring>
#include <unistd.h>
#include <vector>
#include <string>
#include <utility>
#include <algorithm>
#include <unordered_map>
#include <ctime>

#include "popseq_utils.h"

const int nprogeny = 20;

int main(int argc, char *argv[])
{
    srand48(time(NULL));
    
    std::vector<int> calls;
    calls.assign(nprogeny,0);
    std::string str(nprogeny,' ');
    std::string check(nprogeny,' ');
    
    for(auto i=0; i<100; i++)
    {
        for(auto j=0; j<nprogeny; j++)
        {
            if(drand48() < 0.5)
            {
                calls[j] = 0;
                str[j] = '0';
            }
            else
            {
                calls[j] = 1;
                str[j] = '1';
            }
        }
        
        uint32_t code = calls2code(calls);
        code2calls(code,nprogeny,check);
        
        std::cout << str << std::endl;
        std::cout << " " << check << std::endl << std::endl;
    }

    return 0;
}
