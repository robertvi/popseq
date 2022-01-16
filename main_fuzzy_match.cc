//
// load map bincode and generate all possible 1 and 2 error codes
// filter out collisions
// match kmer bincodes to 0-error or 1-error or 2-error code in that order of preference
// output consolidated bin info with each matched kmer
//
// g++ -O3 main_fuzzy_match.cc popseq_utils.cc -std=c++11 -static -o fuzzy_match
// fuzzy_match bin_specs input_file[.lzo|.gz|.bz] > output_file
//
// bin_specs contains canonicalcode chrm phase cm
// inputfile contains kmerseq rg ha p1...p19

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

void generate_errors(int nprogeny,
                     std::unordered_map< uint64_t, int64_t >&code2offset0,
                     std::unordered_map< uint64_t, int64_t >&code2offset1,
                     std::unordered_map< uint64_t, int64_t >&code2offset2)
{
    //generate all 1 error bitmasks
    std::vector< uint64_t > mlist;
    uint64_t mask=1;
    for(auto i=0; i<nprogeny-1; i++)
    {
        mlist.push_back(mask);
        mask <<= 1;
    }

    std::cerr << "1-error masks " << mlist.size() << std::endl;

    //for each existing code
    for(auto cit=code2offset0.begin(); cit!=code2offset0.end(); cit++)
    {
        uint64_t code = cit->first;
        int64_t offset = cit->second;

        //create all 1-error codes from existing codes
        for(auto mit=mlist.begin(); mit!=mlist.end(); mit++)
        {
            uint64_t error = code ^ (*mit);
            auto eit = code2offset1.find(error);

            if(eit != code2offset1.end())
            {
                //collision, flag as unusable
                code2offset1[error] = -1;
            }
            else
            {
                //store
                code2offset1[error] = offset;
            }
        }
    }

    std::cerr << "1-error codes " << code2offset1.size() << std::endl;

    //generate all 2 error bitmasks
    mlist.clear();
    for(auto i=0; i<nprogeny-2; i++)
    {
        for(auto j=i+1; j<nprogeny-1; j++)
        {
            mask = (1<<i) + (1<<j);
            mlist.push_back(mask);
        }
    }

    std::cerr << "2-error masks " << mlist.size() << std::endl;

    //for each existing code
    for(auto cit=code2offset0.begin(); cit!=code2offset0.end(); cit++)
    {
        uint64_t code = cit->first;
        int64_t offset = cit->second;

        //create all 2-error codes from existing codes
        for(auto mit=mlist.begin(); mit!=mlist.end(); mit++)
        {
            uint64_t error = code ^ (*mit);
            auto eit = code2offset2.find(error);

            if(eit != code2offset2.end())
            {
                //collision, flag as unusable
                code2offset2[error] = -1;
            }
            else
            {
                //store
                code2offset2[error] = offset;
            }
        }
    }

    std::cerr << "2-error codes " << code2offset2.size() << std::endl;
}

int load_bin_specs(std::ifstream&inp,
                   std::vector< std::vector<std::string> >&binspec,
                   std::unordered_map< uint64_t, int64_t >&code2offset0)
{
    int nprogeny = 0;

    std::vector< std::string > spec;

    int ctr=0;
    while(1)
    {
        std::string str;

        //code chrm count0 count1 phase0 mincm maxcm meancm
        spec.clear();
        for(unsigned i=0; i<8; i++)
        {
            inp >> str;
            spec.push_back(str);
        }

        if(inp.eof()) break;

        //canonical barcode
        str = spec[0];

        if(nprogeny == 0)
        {
            //get number of progeny from first bincode
            nprogeny = str.length() + 1;
            if(nprogeny > 65)
            {
                std::cerr << "too many progeny to fit into a 64 bit integer" << std::endl;
                exit(1);
            }
        }
        else
        {
            //ensure all subsequent codes are the correct length
            if(str.length() != nprogeny-1)
            {
                std::cerr << "unexpected bincode" << str << std::endl;
                exit(1);
            }
        }

        ctr += 1;

        //convert to canonical binmap code from string to uint64_t
        uint64_t code = str2code(str);

        //store new code and spec
        code2offset0[code] = binspec.size();
        binspec.push_back(spec);
    }

    std::cerr << std::endl << "codes " << code2offset0.size() << std::endl;

    return nprogeny;
}

void dump_output_line(int64_t offset,std::vector< std::vector<std::string> >&binspec,
                      const std::vector<int>&calls,const std::string&seq,const int parent1,const int parent2)
{
    std::vector< std::string >&spec = binspec[offset];
    std::string kmerphase;

    //do we need to take account of which progeny seems to have an error here?!
    //it could be the popseq or map which has the error
    if(calls[0] == 0)
    {
        kmerphase = spec[4]; //progeny1==0 therefore report bin's phase0 directly
    }
    else
    {
        if(spec[4] == "0") kmerphase = "1"; //progeny1==1 therefore report flip phase0 of bin
        else               kmerphase = "0";
    }

    //output: code chrm kmerphase mincm maxcm meancm
    std::cout << spec[0] << " " << spec[1] << " " << kmerphase << " ";
    std::cout << spec[5] << " " << spec[6] << " " << spec[7];

    //output: kmer parent1 parent2 progeny1... progenyN
    std::cout << " " << seq << " " << parent1 << " " << parent2;
    for(auto cit=calls.begin(); cit!=calls.end(); cit++) std::cout << " " << (*cit);

    std::cout << std::endl;
}

void process_calls_fuzzy(std::ifstream&inp,int nprogeny,
                         std::vector< std::vector<std::string> >&binspec,
                         std::unordered_map< uint64_t, int64_t >&code2offset0,
                         std::unordered_map< uint64_t, int64_t >&code2offset1,
                         std::unordered_map< uint64_t, int64_t >&code2offset2)
{
    //read header line
    std::string header;
    getline(inp,header);

    //check number of progeny agrees
    //kmer parent1 parent2 progeny1... progenyN
    int ncheck = std::count(header.begin(), header.end(), ' ') - 2;

    if(ncheck != nprogeny)
    {
        std::cerr << "wrong number of progeny" << std::endl;
        exit(1);
    }

    //vector of calls
    std::vector<int> calls;
    calls.assign(nprogeny,0);

    unsigned ctr=0;
    unsigned exact=0,error1=0,error2=0;
    unsigned ambig1=0,ambig2=0;

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

        std::unordered_map< uint64_t, int64_t >::iterator it;

        //if exact match
        it = code2offset0.find(code);
        if(it != code2offset0.end())
        {
            exact += 1;

            dump_output_line(it->second,binspec,calls,seq,parent1,parent2);

            continue;
        }

        //if match with single error
        it = code2offset1.find(code);
        if(it != code2offset1.end())
        {
            if(it->second == -1)
            {
                ambig1 += 1;
                continue;
            }
            else
            {
                error1 += 1;
                dump_output_line(it->second,binspec,calls,seq,parent1,parent2);
                continue;
            }
        }

        //if match with two errors
        it = code2offset2.find(code);
        if(it != code2offset2.end())
        {
            if(it->second == -1)
            {
                ambig2 += 1;
                continue;
            }
            else
            {
                error2 += 1;
                dump_output_line(it->second,binspec,calls,seq,parent1,parent2);
                continue;
            }
        }
    }

    std::cerr << std::endl;
    std::cerr << "unique matches: " << exact << " " << error1 << " " << error2 << std::endl;
    std::cerr << "ambiguous matches: " << ambig1 << " " << ambig2 << std::endl;
}

int main(int argc, char *argv[])
{
    if(argc != 3)
    {
        std::cout << "usage: fuzzy_match <bin_specs> <input_file[.lzo|.gz|.bz]> > <output_file>"  << std::endl;
        exit(0);
    }

    std::unordered_map< uint64_t, int64_t > code2offset0; //zero error bincodes
    std::vector< std::vector<std::string> > binspec;

    std::string binsfile = argv[1];
    std::ifstream bins(binsfile);

    if(!bins.good())
    {
        std::cerr << "failed to open file " << binsfile << std::endl;
        exit(1);
    }

    //load bin specs
    int nprogeny = load_bin_specs(bins,binspec,code2offset0);

    bins.close();

    //generate all possible 1 and 2 error bincodes, linked back to correct code
    //filter out any collisions
    std::unordered_map< uint64_t, int64_t > code2offset1; //one error bincodes
    std::unordered_map< uint64_t, int64_t > code2offset2; //two error bincodes

    generate_errors(nprogeny,code2offset0,code2offset1,code2offset2);

    //match kmer codes to map codes
    std::string fname = argv[2];
    std::ifstream inp;
    open_compressed_file(fname,inp);

    if(!inp.good())
    {
        std::cerr << "failed to open file " << fname << std::endl;
        exit(1);
    }

    process_calls_fuzzy(inp,nprogeny,binspec,code2offset0,code2offset1,code2offset2);

    //dump_calls(code2counts,nprogeny);

    return 0;
}
