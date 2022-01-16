'''
helper functions relating to kmers
'''

def generate_sequences(filename):
    '''
    generator to yield sequences from a 2 or 3 column flat file
    '''

    with open(filename) as f:
        for line in f:
            column = line.strip().split()
            identifier = column[0]
            sequence = column[1]

            if len(column) > 2: quality = column[2]
            else:               quality = None

            yield identifier,sequence,quality

def generate_kmers(seq,kmer_size):
    '''
    generator to yield kmers from a sequence
    '''

    seqlen = len(seq)

    for i in range(0,seqlen-kmer_size+1):
        yield seq[i:i+kmer_size]

def reverse_compliment(kmer):
    '''
    return reverse complement of the kmer
    only 'ATCGN' are considered legitimate bases
    '''

    rev = ''

    for ch in kmer:
        if   ch == 'A': rev = 'T' + rev
        elif ch == 'T': rev = 'A' + rev
        elif ch == 'C': rev = 'G' + rev
        elif ch == 'G': rev = 'C' + rev
        elif ch == 'N': rev = 'N' + rev
        else:           raise Exception

    return rev
