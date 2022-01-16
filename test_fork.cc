//
// test decompressing .lzo files automatically
//
// g++ -O3 test_fork.cc -std=c++11 -static -o test_fork
//

#include <istream>
#include <iostream>
#include <fstream>
#include <cstring>
#include <unistd.h>
#include <vector>
#include <string>
#include <utility>
#include <sys/types.h>

//argv[1] = filename

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

void open_compressed_file(const std::string&fname,std::ifstream&inp)
{
    std::string decomp = "";
    
    if(endswith(fname,".lzo")) decomp = "lzop -dcf";
    else if(endswith(fname,".gz")) decomp = "zcat";
    else if(endswith(fname,".bz2")) decomp = "bzcat";

    if(decomp == "")
    {
        //no compression detected
        inp.open(fname);
        return;
    }

    //pipe decompressed stream through a fifo
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

int main(int argc, char *argv[])
{
    if(argc < 2)
    {
        std::cout << "usage: test_fork <filename>"  << std::endl;
        exit(0);
    }

    std::string fname = argv[1];
    
    std::ifstream inp;
    
    open_compressed_file(fname,inp);
    
    while(1)
    {
        std::string line;
        getline(inp,line);
        if(!inp.good()) break;
        std::cout << line;
    }
    
    inp.close();
    return 0;
}
