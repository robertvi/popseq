//
// process kmer counts from both parents and progeny
// from first parent retain only high confidence 1-copy kmers in hash
// from remainder samples only consider kmers hq 1-copy in first parent
// classify remainder samples with respect to copy number
// output kmer strings with classification in all progeny
//
// g++ -O3 main_process_kmers.cc popseq_utils.cc -std=c++11 -static -o process_kmers
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

//argv[1] = config file

/* example config file (remove the comments)
subsample 0.01                   #process a random subsample of kmers, 1.0 to process all
seed 0                           #random number generator seed, 0 means use system time (in seconds)
parental_peak_fraction 0.6       #central fraction of parental kmer peaks to treat as confident copy number, max is 1.0
progeny_peak_fraction 1.0        #central fraction of progeny kmer peaks to treat as confident copy number, max is 1.0
filename_suffix _kmers           #suffix to apply to sample names to generate the kmer file name, - to use no suffix
parent1 5 10 15 20 25            #parent1 position of minimum1 peak1 minimum2 peak2 minimum3 in the kmer histogram plot as integer kmer counts
parent2 5 10 15 20 25              #same for parent2
progeny1 5 10 15 20 25              #same for progeny1
....                             #one line per progeny then file ends
*/

#include "popseq_utils.h"

struct config
{
    double subsample;
    uint64_t seed;
    double par_frac;
    double prog_frac;
    std::string suffix;
};

void load_config(std::ifstream&cfile,struct config&conf)
{
    std::string dummy;

    cfile >> dummy >> conf.subsample;
    cfile >> dummy >> conf.seed;
    cfile >> dummy >> conf.par_frac;
    cfile >> dummy >> conf.prog_frac;
    cfile >> dummy >> conf.suffix;

    if(conf.suffix == "-") conf.suffix = "";

    if(conf.seed == 0) srand48(time(NULL));
    else               srand48(conf.seed);
}

//input the location of the first two peaks and troughs (maxima and minima)
//in the kmer frequency plot, plus the fraction of each peak to treat as high confidence
//output the boundaries of the high confidence regions
void calc_interval(double trough1,double peak1,double trough2,double peak2,double trough3,
                   double frac,
                   uint64_t&max0,uint64_t&min1,uint64_t&max1,uint64_t&min2,uint64_t&max2)
{
    //calculate the kmer count intervals for high confidence 0,1,2-copy kmers
    max0 = round( (trough1 - 1.0)   * frac + 1.0 );
    min1 = round( (trough1 - peak1) * frac + peak1 );
    max1 = round( (trough2 - peak1) * frac + peak1 );
    min2 = round( (trough2 - peak2) * frac + peak2 );
    max2 = round( (trough3 - peak2) * frac + peak2 );
}

int process_sample(struct config&conf,
                   std::ifstream&cfile,
                   std::unordered_map<uint64_t, uint64_t>&kmer2copy,
                   int sampleno,
                   std::vector<std::string>&names)
{
    std::string sample;
    std::string fname;
    double trough1,peak1,trough2,peak2,trough3;
    uint64_t max0,min1,max1,min2,max2;

    //load sample name and minima and maxima positions
    //trough1 is the minimum between the error kmers and the first peak
    //peak1 is the peak on the 1-copy kmers (for a diploid heterozygote)
    //peak2 is the peak on the 2-copy kmers (for a diploid heterozygote)
    cfile >> sample >> trough1 >> peak1 >> trough2 >> peak2 >> trough3;

    if(cfile.eof()) return 1; //signal end of file, no more samples to process

    //record sample name
    names.push_back(sample);

    fname = sample + conf.suffix;

    //calculate the kmer count interval for the high confidence 0,1,2-copy kmers
    double frac;
    if(sampleno <= 2) frac = conf.par_frac;
    else              frac = conf.prog_frac;

    calc_interval(trough1,peak1,trough2,peak2,trough3,frac,max0,min1,max1,min2,max2);

    std::cerr << "------------------------" << std::endl;
    std::cerr << "sample number " << sampleno << std::endl;
    std::cerr << "name " << sample << std::endl;
    std::cerr << "file " << fname << std::endl;
    std::cerr << "max0 " << max0 << std::endl;
    std::cerr << "min1 " << min1 << std::endl;
    std::cerr << "max1 " << max1 << std::endl;
    std::cerr << "min2 " << min2 << std::endl;
    std::cerr << "max2 " << max2 << std::endl;

    std::ifstream inp;

    open_compressed_file(fname,inp);

    if(!inp.good())
    {
        std::cerr << "failed to open file " << fname << std::endl;
        exit(1);
    }

    int ctr=0;

    while(1)
    {
        std::string seq;
        uint64_t count,kmer;

        inp >> seq;
        if(seq.length() != K || inp.eof()) break;

        inp >> count;
        ctr += 1;

        if(ctr % 1000000 == 0) std::cerr << "#";
        if(ctr % 50000000 == 0) std::cerr << std::endl;

        //parent 1
        if(sampleno == 1)
        {
            //accept only a random subsample of kmers from parent1
            if(conf.subsample < 1.0) if(drand48() >= conf.subsample)
            {
                //std::cerr << "kmer not sampled" << std::endl;
                continue;
            }

            //reject any kmer that's not high confidence 1-copy
            if(count < min1 || count > max1)
            {
                //std::cerr << "reject kmer " << std::endl;
                continue;
            }

            //convert kmer from character string to 64 bit integer
            kmer = string2kmer(seq.c_str());

            kmer2copy[kmer] = 1;
            //std::cerr << "accept kmer" << std::endl;
        }
        //parent 2
        else if(sampleno == 2)
        {
            //if 0-copy (error kmer) nothing to do:
            //if already present in hash then parent 2 count is already initialised to zero
            //if not present then no need to add it as it is not present in parent 1
            if(count <= max0) continue;

            //convert kmer from character string to 64 bit integer
            kmer = string2kmer(seq.c_str());

            //ambiguous between error and 1-copy kmer
            //or ambiguous between 1-copy and >1 copy
            //or high confidence >1 copy
            //if present in hash delete it as parent 2 is not callable or >1 copy
            if( (count > max0 && count < min1) || count > max1)
            {
                kmer2copy.erase(kmer);
                continue;
            }

            //kmer is high confidence 1-copy in parent 2
            //determine if parent 1 is 1-copy as well
            std::unordered_map<uint64_t,uint64_t>::iterator it = kmer2copy.find(kmer);

            //if parent 1 is not 1-copy drop this kmer
            //because at this stage we are only interested in parent 1 kmers
            //however in the general case these could be retained if they are 0-copy in parent 1
            if(it == kmer2copy.end()) continue;

            //record parent 2 is 1-copy
            it->second += 1 << 2;
        }
        //simplified calling, setting to presence / absence or missing
        //ignores fraction parameter
        else
        {
            //if 0-copy (error kmer) nothing to do:
            //if already present in hash then count is already initialised to zero
            //if not present then do not add it
            if(count < trough1) continue;

            //convert kmer from character string to 64 bit integer
            kmer = string2kmer(seq.c_str());

            //find if kmer is in the hash, skip kmer if not present
            std::unordered_map<uint64_t,uint64_t>::iterator it = kmer2copy.find(kmer);
            if(it == kmer2copy.end()) continue;

            //set call to 1 or missing (3)
            uint64_t call;
            if(count > trough1)
            {
                call = 1; //1-copy
            }
            else //call == trough1 (ambiguous copy number)
            {
                call = 3;
            }

            it->second += call << (2*(sampleno-1));
        }
#if 0
        //original version using the fraction parameter
        else
        {
            //if 0-copy (error kmer) nothing to do:
            //if already present in hash then count is already initialised to zero
            //if not present then do not add it
            if(count <= max0) continue;

            //convert kmer from character string to 64 bit integer
            kmer = string2kmer(seq.c_str());

            //find if kmer is in the hash, skip kmer if not present
            std::unordered_map<uint64_t,uint64_t>::iterator it = kmer2copy.find(kmer);
            if(it == kmer2copy.end()) continue;

            //set call to 1,2 or missing
            uint64_t call=0;
            if(count >= min1 && count <= max1)
            {
                call = 1; //1-copy
            }
            else if(count >= min2 && count <= max2)
            {
                call = 2; //2-copy
            }
            else
            {
                call = 3; //missing (ambiguous copy number)
            }

            it->second += call << (2*(sampleno-1));
        }
#endif
    }

    std::cerr << std::endl;

    inp.close();

    std::cerr << fname << ' ' << ctr << " lines processed" << std::endl;
    std::cerr << kmer2copy.size() << " kmers in hash" << std::endl;

    return 0;
}

void output_results(std::unordered_map<uint64_t, uint64_t>&kmer2copy,std::vector< std::string >&names)
{
    //header
    std::cout << "kmer";
    for(auto it=names.begin(); it!=names.end(); it++) std::cout << " " << (*it);
    std::cout << std::endl;

    //kmer strings and calls
    std::string seq(31,' ');
    for(auto it=kmer2copy.begin(); it!=kmer2copy.end(); it++)
    {
        kmer2seq(it->first,seq);
        std::cout << seq;
        uint64_t call = it->second;

        for(int i=0; i<names.size(); i++)
        {
            std::cout << ' ' << (call & 0x3); //copy number 0,1,2, 3==missing
            call >>= 2;
        }
        std::cout << std::endl;
    }
}

int main(int argc, char *argv[])
{
    if(argc < 2)
    {
        std::cout << "usage: process_kmers <config_file> > <output_file>"  << std::endl;
        exit(0);
    }

    //load config from file
    std::ifstream cfile(argv[1]);
    struct config conf;
    load_config(cfile,conf);

    std::cerr << "subsample " << conf.subsample << std::endl;
    std::cerr << "seed " << conf.seed << std::endl;
    std::cerr << "parental peak fraction " << conf.par_frac << std::endl;
    std::cerr << "progeny peak fraction " << conf.prog_frac << std::endl;
    std::cerr << "filename suffix " << conf.suffix << std::endl;

    //key is 64 bit encoded canonical kmer
    //value is 2 bits per sample (two parents plus up to 31 progeny)
    std::unordered_map<uint64_t, uint64_t> kmer2copy;

    //record sample names
    std::vector< std::string > names;

    //process first parent, retain only high confidence 1-copy kmers
    //since this is the genome being assembled only haplotype specific kmers will be considered initially
    process_sample(conf,cfile,kmer2copy,1,names);

    //process second parent, retain kmers which are 1-copy in one or both parents
    process_sample(conf,cfile,kmer2copy,2,names);

    //process remaining samples until end of config file is reached
    int ctr=3;
    while(1)
    {
        //process progeny, classify each kmer as presence/absence/uncalled
        if(process_sample(conf,cfile,kmer2copy,ctr,names)) break;
        ctr += 1;
        if(ctr == 31)
        {
            std::cerr << "31 is maximum number of progeny due to using 2 bits per sample in 64 bit integer"
                      << std::endl;
            exit(1);
        }
    }

    cfile.close();

    //write results to stdout
    output_results(kmer2copy,names);

    return 0;
}
