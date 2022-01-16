#include "popseq_utils.h"

const uint64_t KMER_MASK=( ((uint64_t)1)<<(K*2) )-1;
const uint64_t KMER_FIRSTOFFSET=(K-1)*2;

#include <unistd.h>

//test if string ends with a given suffix
bool endswith(const std::string&str,const std::string&suf)
{
    if(str.length() >= suf.length())
    {
        if(str.compare(str.length()-suf.length(),suf.length(),suf) == 0)
        {
            return true;
        }
    }

    return false;
}

//test if string starts with a given prefix
bool startswith(const std::string&str,const std::string&pre)
{
    if(str.length() >= pre.length())
    {
        if(str.compare(0,pre.length(),pre) == 0)
        {
            return true;
        }
    }

    return false;
}

//return an ifstream to the decompressed version of a file
//decide decompression method based on the suffix
void open_compressed_file(const std::string&fname,std::ifstream&inp)
{
    //test if file exists
    std::ifstream tmp;
    tmp.open(fname);
    if(!tmp.good())
    {
        std::cerr << "unable to open file " << fname << std::endl;
        exit(1);
    }
    tmp.close();

    //decide if file needs decompressing
    std::string decomp = "";
    if(endswith(fname,".lzo")) decomp = "lzop -dcf";  //lzop
    else if(endswith(fname,".gz")) decomp = "zcat";   //gzip
    else if(endswith(fname,".bz2")) decomp = "bzcat"; //bzip2

    if(decomp == "")
    {
        //no compression detected, open file as normal
        inp.open(fname);
        return;
    }

    //pipe decompressed stream through a fifo to allow parallel decompression
    std::string fifo = fname + "_t_m_p_f_i_f_o";
    std::string cmd;
    cmd = "mkfifo " + fifo;
    system(cmd.c_str());

    if(fork() == 0)
    {
        //child process
        cmd = decomp + " " + fname + " > " + fifo + " ; rm " + fifo;
        system(cmd.c_str());
        exit(0);
    }

    inp.open(fifo);
}

//assumes string is 31mer kmer with no N's but need not be canonical
uint64_t string2kmer(const char*const seq)
{
    uint64_t fkmer=0,rkmer=0;

    for(auto p=0; p<K; p++)
    {
        //fkmer: grows from the right (LSB)
        //rkmer: grows from the left (MSB)
        switch(seq[p])
        {
            case 'A':
            case 'a':
                fkmer = ((fkmer << 2) +  uint64_t(0)) & KMER_MASK;
                rkmer =  (rkmer >> 2) + (uint64_t(3) << KMER_FIRSTOFFSET);
                break;
            case 'C':
            case 'c':
                fkmer = ((fkmer << 2) +  uint64_t(1)) & KMER_MASK;
                rkmer =  (rkmer >> 2) + (uint64_t(2) << KMER_FIRSTOFFSET);
                break;
            case 'G':
            case 'g':
                fkmer = ((fkmer << 2) +  uint64_t(2)) & KMER_MASK;
                rkmer =  (rkmer >> 2) + (uint64_t(1) << KMER_FIRSTOFFSET);
                break;
            case 'T':
            case 't':
                fkmer = ((fkmer << 2) +  uint64_t(3)) & KMER_MASK;
                rkmer =  (rkmer >> 2) + (uint64_t(0) << KMER_FIRSTOFFSET);
                break;
            default:
                std::cerr << "unexpected character: " << seq[p] << std::endl;
                exit(1);
        }
    }

    if (fkmer <= rkmer) return fkmer;
    else                return rkmer;
}

//convert canonical bincode from string to uint64_t
//string[0] => MSB
uint64_t str2code(const std::string&str)
{
    uint64_t code=0;

    for(auto p=0; p<str.size(); p++)
    {
        code <<= 1;
        if(str[p] == '1') code += 1;
    }

    return code;
}

//encode calls into a bitstring, bits are relative to call of progeny 1
//0==same progeny 1
//1==different to progeny 1
//all bits encoded wrt call[0]
//calls[1] => MSB
uint64_t calls2code(const std::vector<int>&calls)
{
    uint64_t code=0;

    for(auto p=1; p<calls.size(); p++)
    {
        code <<= 1;
        if(calls[p] != calls[0]) code += 1;
    }

    return code;
}

//canonical code is nprogeny-1 in length, preallocate the calls string
//MSB => calls[0]
void code2calls(uint64_t code,int progeny,std::string&calls)
{
    uint64_t mask = 1<<(progeny-2);

    for(auto p=0; p<progeny-1; p++)
    {
        if(mask & code) calls[p] = '1';
        else            calls[p] = '0';
        mask >>= 1;
    }
}

//convert kmer to string, string should be preallocated to length of 31
void kmer2seq(uint64_t kmer,std::string&seq)
{
    for(auto i=0; i<31; i++)
    {
        switch(kmer & uint64_t(3))
        {
            case 0:
                seq[30-i] = 'A';
                break;
            case 1:
                seq[30-i] = 'C';
                break;
            case 2:
                seq[30-i] = 'G';
                break;
            case 3:
                seq[30-i] = 'T';
                break;
        }

        kmer >>= 2;
    }
}
