//
// deal only with RG-1 HA-0 kmers initially
// group calls (including repulsion phase)
// output call barcode patterns and counts of coupling and repulsion phase
// 
//
// g++ -O3 main_process_calls.cc popseq_utils.cc -std=c++11 -static -o process_calls
// process_calls input_file[.lzo|.gz|.bz] > output_file
//

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

int process_calls(std::ifstream&inp,std::unordered_map< uint64_t, std::pair<uint32_t,uint32_t> >&code2counts)
{
    //skip header line
    std::string header;
    getline(inp,header);
    
    int nprogeny = std::count(header.begin(), header.end(), ' ') - 2;
    
    if(nprogeny > 65)
    {
        std::cerr << "too many progeny to fit into a 64 bit integer" << std::endl;
        exit(1);
    }
    
    //vector of calls
    std::vector<int> calls;
    calls.assign(nprogeny,0);

    int ctr=0;
    while(1)
    {
        std::string seq;
        int parent1,parent2;

        //read the kmer
        inp >> seq;
        if(seq.length() != K || inp.eof()) break;
        
        //read calls of parents and progeny
        inp >> parent1 >> parent2;
        bool missing_flag = false;
        for(int i=0; i<nprogeny; i++)
        {
            inp >> calls[i];
            if(calls[i] == 3) missing_flag = true; //missing data is present
            //if(calls[i] != 0) calls[i] = 1; //convert to presence/absence
        }
        
        if(inp.eof())
        {
            std::cerr << "unexpected end of file" << std::endl;
            break;
        }
        
        ctr += 1;
        if(ctr % 1000000 == 0) std::cerr << "#"; //every million
        if(ctr % 50000000 == 0) std::cerr << std::endl; //every 50 million

        //ignore if not a Parent1-1copy Parent2-0copy kmer
        if(!(parent1 == 1 && parent2 == 0)) continue;
        
        //ignore if contains any missing calls
        if(missing_flag) continue;
        
        //convert to canonical binmap code encoded wrt progeny1's call
        uint64_t code = calls2code(calls);
        
        std::unordered_map< uint64_t, std::pair<uint32_t,uint32_t> >::iterator it = code2counts.find(code);
        
        //if code is new create new record otherwise increment existing count
        if(calls[0] == 0)
        {
            //progeny1 == 0
            if(it == code2counts.end()) code2counts[code] = std::make_pair(1,0);
            else                        it->second.first += 1;
        }
        else
        {
            //progeny1 == 1
            if(it == code2counts.end()) code2counts[code] = std::make_pair(0,1);
            else                        it->second.second += 1;
        }
    }
    
    std::cerr << std::endl << "codes " << code2counts.size() << std::endl;
    
    return nprogeny;
}

void dump_calls(std::unordered_map< uint64_t, std::pair<uint32_t,uint32_t> >&code2counts, int nprog)
{
    std::unordered_map< uint64_t, std::pair<uint32_t,uint32_t> >::iterator it;
    std::string calls(nprog-1,' ');
    
    for(it=code2counts.begin(); it!=code2counts.end(); it++)
    {
        code2calls(it->first,nprog,calls);
        //std::cout << it->first << " " << it->second.first << " " << it->second.second << std::endl;
        std::cout << calls << " " << it->second.first << " " << it->second.second << std::endl;
    }
}

int main(int argc, char *argv[])
{
    if(argc != 2)
    {
        std::cout << "usage: process_calls <input_file[.lzo|.gz|.bz]> > <output_file>"  << std::endl;
        exit(0);
    }

    std::string fname = argv[1];

    std::ifstream inp;
    open_compressed_file(fname,inp);

    if(!inp.good())
    {
        std::cerr << "failed to open file " << fname << std::endl;
        exit(1);
    }
    
    std::unordered_map< uint64_t, std::pair<uint32_t,uint32_t> > code2counts;
    
    int nprogeny = process_calls(inp,code2counts);
    
    dump_calls(code2counts,nprogeny);

    return 0;
}
