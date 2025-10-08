#include <iostream>
#include <string.h>
#include <vector>
#include <cstring>
#include <stdint.h>
#include <fstream>
#include <sstream>
#include <zlib.h>
#include <filesystem>
#include <fstream>
#include "vec.hpp"
#include "str.hpp"
#include "kmer.hpp"
#include "kseq.h"

/* Assumptions:
 * - GFA comes from a compacted de Bruijn graph (eg, Bifrost)
 * - this means that every k-mer appears *once*
 */

#define MCALLOC(ptr, len) ((ptr) = (__typeof__(ptr))calloc((len), sizeof(*(ptr))))
#define MMALLOC(ptr, len) ((ptr) = (__typeof__(ptr))malloc((len) * sizeof(*(ptr))))
#define MREALLOC(ptr, len) ((ptr) = (__typeof__(ptr))realloc((ptr), (len) * sizeof(*(ptr))))

KSEQ_INIT(gzFile, gzread)

typedef struct fasta_t {
    size_t nk;
    std::vector<size_t> lens;
    std::vector<char*> sequences;
    std::vector<std::string> names;
} fasta_t;

fasta_t read_fasta(const char *fn) {
    kseq_t *ks;
	gzFile fp;
    fasta_t fa;
    memset(&fa, 0, sizeof(fasta_t));

	if ((fp = gzopen(fn, "r")) == 0) return fa;
    ks = kseq_init(fp);

    int ret;
    while ((ret = kseq_read(ks)) >= 0) {
        size_t l = ret; // get len sequence
        char* seq = (char*) malloc(l);
        memcpy(seq, ks->seq.s, l);
        fa.sequences.push_back(seq);
        fa.names.push_back(std::string(ks->name.s, ks->name.s+ks->name.l));
        fa.lens.push_back(l);
        fa.nk += l;
    }
    return fa;
}

void error(const char *msg) {
    std::cerr << "Error: " << msg << '\n' << std::flush;
    exit(1);
}

std::string base_name(std::string const & path) {
    return path.substr(path.find_last_of("/\\") + 1);
}

std::string remove_extension(std::string const & filename) {
    typename std::string::size_type const p(filename.find_last_of('.'));
    return p > 0 && p != std::string::npos ? filename.substr(0, p) : filename;
}

void read_gfa(const char* file, Vec<str>& nodes){
    std::cerr << "reading file..." << std::flush;
    std::ifstream fin(file);
    if (!fin.is_open()) {
        std::cerr << "opening file: " << file << '\n';
        std::exit(EXIT_FAILURE);
    }

    std::string line;
    std::stringstream ss; 
    nodes.push({"empty", 5});

    while (std::getline(fin, line, '\n')) {
        std::string type;
        ss.str(line);

        ss >> type;
        if (type[0] == 'S') {
            std::string node, _;
            ss >> _ >> node;
            str node_str(node);
            nodes.push(node_str);
        }

        ss.clear();
    }
    std::cerr << "ok\n" << std::flush;
}

void gfa2gff(kmertable_t *kmer_table, std::string filepath, int k, Vec<str>& nodes) {
    std::string filename = remove_extension(base_name(filepath));
    //std::string filename = base_name(filepath);
    std::cerr << "printing " << filename << '\n' << std::flush;

    fasta_t fa = read_fasta(filepath.c_str());

    kmer_t mask = (1ULL << (2*k)) - 1;

    for (int z = 0; z < fa.lens.size(); z++) {
        kmer_t kmer = 0;
        int len = 0;
        uint64_t last_rid = UINT64_MAX;
        uint64_t last_pos = UINT64_MAX;
        int last_strand = 2;
        uint64_t start = 0;
        uint64_t end = 0;
        char* sequence = fa.sequences[z];
        size_t len_sequence = fa.lens[z];
        std::string seqname = fa.names[z];
        int64_t i = 0;
        int64_t j = 0;

        std::cout << "P\t"<< filename << "#" << seqname << '\t' << std::flush;

        bool is_first_node = true;
        while (i < len_sequence) {
            unsigned char c = nt_2_bits[sequence[i]];
            //std::cerr << "out" << sequence[i] << " " << i << "\n";
            if (c < 4) {
                kmer = (kmer << 2 | c) & mask;
                if (len+1 < k) len++;
                else {
                    auto alg = kmer_table->align(kmer);
                    if (alg.pos == UINT32_MAX) {
                        std::cerr << "\nError: kmer " << bits2kmer(kmer, k) << " not found" << '\n' << std::flush;
                        exit(1);
                    }
                    str node = nodes[alg.rid];
                    size_t len_node = node.n;
                    start = i - k + 1;
                    if (!alg.strand) {
                        j = alg.pos;

                        while (i+1 < len_sequence && j+1 < len_node && 
                                nt_2_bits[sequence[i+1]] == nt_2_bits[node[j+1]]){
                            i++;
                            j++;
                            c = nt_2_bits[sequence[i]];
                            kmer = (kmer << 2 | c) & mask;
                        }
                        end = i;
                    } else {
                        j = alg.pos - k + 1;

                        while (i+1 < len_sequence && j-1 >= 0 && nt_2_bits[sequence[i+1]] < 4 &&
                                ((uint64_t)(3-nt_2_bits[sequence[i+1]])&3ULL) == (uint64_t)nt_2_bits[node[j-1]]){
                            i++;
                            j--;
                            c = nt_2_bits[sequence[i]];
                            kmer = (kmer << 2 | c) & mask;
                        }
                        end = i;
                    }
                    if (!is_first_node) std::cout << ',';
                    is_first_node = false;
                    std::cout << alg.rid << (alg.strand? '-' : '+');
                }
            } 
            if (nt_2_bits[sequence[i]] == 4) {
                start = 0;
                len = 0, kmer = 0;
            }
            i++;
        }
        std::cout << "\t*\n" << std::flush;
    }
}

void read_file(const char* file, std::string& fa){
    std::cerr << "reading file..." << std::flush;
    std::ifstream fin(file);
    if (!fin.is_open()) {
        std::cout << "Error opening file: " << file << '\n';
        std::exit(EXIT_FAILURE);
    }

    std::string line;
    std::stringstream ss; 
    while (std::getline(fin, line, '\n')) {
        std::string tmp;
        ss.str(line);
        ss >> tmp;
        fa += tmp;
        ss.clear();
    }
    std::cerr << "complete\n" << std::flush;
}

int main(int argc, char **argv) {
    Vec<str> nodes;

    if (argc < 4) {
        std::cerr << "Usage:\n"
                  << "  " << argv[0] << " <k> <gfa> <fasta> [more fasta ...] [options]\n\n"
                  << "Required arguments:\n"
                  << "  <k>       Integer, k-mer size.\n"
                  << "  <gfa>     Path to GFA file.\n"
                  << "  <fasta>   Path to at least one FASTA file.\n"
                  << "            You may provide multiple FASTA files.\n\n"
                  << "Options:\n"
                  << "  -t, --threads <num>   Number of threads to use (default: number of cores).\n"
                  << "  -h, --help            Show this help message.\n";
        exit(1);
    }

    Vec<char*> args;
    args.push(argv[0]); 

    NUM_THREADS = omp_get_max_threads(); 

    for (int i = 1; i < argc; i++) {
        std::string arg = argv[i];

        if (arg == "--threads" || arg == "-t") {
            if (i + 1 < argc) {
                NUM_THREADS = std::stoi(argv[++i]); // consume next argument
            } else {
                std::cerr << "Error: --threads requires a number" << std::endl;
                return 1;
            }
        } else {
            args.push(argv[i]);
        }
    }

    int k = std::stoi(args[1]);
    std::cerr << "Running on " << NUM_THREADS << " threads" << '\n';

    std::string gfa_file = args[2];
    read_gfa(gfa_file.c_str(), nodes);
    kmertable_t *kmer_table = count_kmers(nodes, k);
    std::cerr << "Finished creating hashtable of kmers: " << kmer_table->num_kmers << " kmers found" << '\n';
    
    for (int i = 3; i < args.size(); i++) {
        std::cerr << "["<<i-2 << "/" << args.size()-3<<"] " << std::flush;
        gfa2gff(kmer_table, args[i], k, nodes);
    }

    for (int i = 0; i < nodes.size(); i++) nodes[i].clean();
    nodes.clean();

    return 0;
}
