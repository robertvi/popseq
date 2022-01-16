#ifndef __RJV_POPSEQ_UTILS_H__
#define __RJV_POPSEQ_UTILS_H__

const int K=31; //kmer length (bases)

#include <string>
#include <fstream>
#include <iostream>
#include <vector>

//test if string ends with a given suffix
bool endswith(const std::string&str,const std::string&suf);

//test if string starts with a given prefix
bool startswith(const std::string&str,const std::string&pre);

//return an ifstream to the decompressed version of a file
//decide decompression method based on the suffix
void open_compressed_file(const std::string&fname,std::ifstream&inp);

//assumes string is 31mer kmer with no N's but need not be canonical
uint64_t string2kmer(const char*const seq);

//convert kmer to string, string should be preallocated to length of 31
void kmer2seq(uint64_t kmer,std::string&seq);

//convert progeny calls into "canonical" binmap code
uint64_t calls2code(const std::vector<int>&calls);

//convert "canonical" binmap code to progeny calls
void code2calls(uint64_t code,int progeny,std::string&calls);

//convert canonical bincode from string to uint64_t
uint64_t str2code(const std::string&str);
#endif
