#include <iostream>
#include <stdint.h>
#include <fstream>
#include <sstream>
#include <zlib.h>
#include <filesystem>
#include <fstream>
#include "vec.hpp"
#include "str.hpp"
#include "kmer.hpp"

/* Assumptions:
 * - GFA comes from a compacted de Bruijn graph (eg, Bifrost)
 * - this means that every k-mer appears *once*
 * - filename of fasta are different, since I use that as the name 
 *   of the path
 */

#define BUFLEN      16384

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

void reconstruct_genome(kmertable_t *kmer_table, std::string dna, int k, Vec<str>& nodes) {
    char bases[] = {'T','G','C','A'};
    kmer_t kmer = 0;
    int len = 0;
    kmer_t mask = (1ULL << (2*k)) - 1;

    uint64_t last_rid = UINT64_MAX;
    int last_strand = 2;
    Vec<char> strands;
    Vec<char> small_dna;
    for (uint64_t i = 0; i < dna.size(); i++) {
        unsigned char c = nt_2_bits[dna[i]];
        if (c < 4) {
            kmer = (kmer << 2 | c) & mask;
            if (len+1 < k) { len++; std::cout << dna[i]; }
            else {
                small_dna.clean();
                auto alg = kmer_table->align(kmer);
                if (!alg.strand) {
                    std::cout << nodes[alg.rid][alg.pos];
                    strands.push('+');
                }
                else {
                    //You take the start of the kmer not the end
                    //query TAATCCTGCGAATGGCAACACAAAGCCGATG
                    //                             pos-k+1                       pos
                    //                             *                             *
                    //rc                           CATCGGCTTTGTGTTGCCATTCGCAGGATTA
                    //   AATATCTGATCGCCCTGACCAACAAACATCGGCTTTGTGTTGCCATTCGCAGGATTAA
                    std::cout << bases[nt_2_bits[(char)nodes[alg.rid][alg.pos-k+1]]];
                    strands.push('-');
                }
            }
        } else {
            std::cout << dna[i];
            len = 0, kmer = 0;
        }
    }
    std::cout << '\n' << std::flush;

}

void print_path_string(kmertable_t *kmer_table, std::string dna, int k) {
    kmer_t kmer = 0;
    int len = 0;
    kmer_t mask = (1ULL << (2*k)) - 1;

    std::cout << "path\t" << std::flush;
    uint64_t last_rid = UINT64_MAX;
    int last_strand = 2;
    for (uint64_t i = 0; i < dna.size(); i++) {
        unsigned char c = nt_2_bits[dna[i]];
        if (c < 4) {
            kmer = (kmer << 2 | c) & mask;
            if (len+1 < k) len++;
            else {
                auto alg = kmer_table->align(kmer);
                if (alg.rid != last_rid || alg.strand != last_strand) {
                    last_rid = alg.rid;
                    last_strand = alg.strand;
                    std::cout << alg.rid << (alg.strand ? '-' : '+') << ',' << std::flush;
                }
            }
        } else len = 0, kmer = 0;
    }
    std::cout << "\t*\n" << std::flush;
}

void print_path_fasta(kmertable_t *kmer_table, std::string filepath, int k, bool break_path_at_non_ACGT_character) {
    std::string filename = remove_extension(base_name(filepath));
    std::cerr << "printing " << filename << '\n' << std::flush;

    kmer_t kmer = 0;
    int len = 0;
    kmer_t mask = (1ULL << (2*k)) - 1;
    uint64_t last_rid = UINT64_MAX;
    uint64_t last_pos = UINT64_MAX;
    bool empty_path = true;
    int last_strand = 2;
    int counter_comments = -1;
    int contiguous_string = 0;

    gzFile file;
    static char buf[BUFLEN];
    if ((file = gzopen(filepath.c_str(), "rb")) == 0) error("opening file");
    
    bool is_first_node = false;
    for (;;) {
        int buf_len = gzread(file, buf, sizeof(buf));

        if (buf_len < 0) error("reading gzip file");
        else if (buf_len == 0) goto finish;

        int z = 0;
        while (z < buf_len) {
            if (buf[z] == '>' || buf[z] == '@') {
                do {
                    z++;
                    if (z >= buf_len) {
                        buf_len = gzread(file, buf, sizeof(buf));
                        if (buf_len < 0) error("reading gzip file");
                        else if (buf_len == 0) goto finish;
                        z = 0;
                    }
                } while (buf[z] != '\n');
                counter_comments++;
                contiguous_string = 0;
                if (counter_comments) std::cout << "\t*\n" << std::flush;
                std::cout << "P\t"<< filename << "#" << counter_comments << "#" << contiguous_string << '\t' << std::flush;
                is_first_node = true;
                z++;
                len = 0, kmer = 0;
            }

            if (z < buf_len && buf[z] != '\n') {
                unsigned char c = nt_2_bits[buf[z]];
                if (c < 4) {
                    kmer = (kmer << 2 | c) & mask;
                    if (len+1 < k) len++;
                    else {
                        auto alg = kmer_table->align(kmer);
                        if (alg.pos == UINT32_MAX) {
                            std::cerr << "\nError: kmer " << bits2kmer(kmer, k) << " not found" << '\n' << std::flush;
                            exit(1);
                        }
                        if (alg.rid != last_rid || alg.strand != last_strand || (!alg.strand && alg.pos <= last_pos) || (alg.strand && alg.pos >= last_pos)) {
                            last_rid = alg.rid;
                            last_strand = alg.strand;
                            last_pos = alg.pos;
                            if (!is_first_node) std::cout << ',';
                            is_first_node = false;
                            std::cout << alg.rid << (alg.strand ? '-' : '+');
                            empty_path = false;
                        }
                    }
                } else { 
                    if (!empty_path && break_path_at_non_ACGT_character) { 
                        contiguous_string++;
                        std::cout << "\t*\n" << std::flush;
                        std::cout << "P\t"<< filename << "#" << counter_comments << "#" << contiguous_string << '\t' << std::flush;
                        is_first_node = true;
                        empty_path = true;
                    }
                    len = 0, kmer = 0;
                }
            }
            z++;
        }
    }

finish:
    if (gzclose(file) != Z_OK) error("closing file");

    std::cout << "\t*\n" << std::flush;
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

void local_test() {
    std::string gfa_file = "/home/luca/@focus/util/print_paths_cdbg/test/graph_k31.gfa";
    Vec<str> nodes;
    read_gfa(gfa_file.c_str(), nodes);
    std::filesystem::path folder_path = "/home/luca/@focus/util/print_paths_cdbg/test/oneline";
    int k = 31;
    kmertable_t *kmer_table = count_kmers(nodes, k);
    std::cerr << "Finished creating hashtable of kmers: " << kmer_table->num_kmers << " kmers found" << '\n';
    for (const auto& entry : std::filesystem::directory_iterator(folder_path)) {
        if (!std::filesystem::is_regular_file(entry)) continue;
        std::filesystem::path input_file = entry.path();
        std::filesystem::path output_file = input_file;
        if (input_file.extension() != ".fa") continue;
        output_file.replace_extension(".res");
        std::streambuf* orig_buf = std::cout.rdbuf();
        std::ofstream outfile(output_file);
        std::cout.rdbuf(outfile.rdbuf());

        std::string dna;
        read_file(input_file.c_str(), dna);
        reconstruct_genome(kmer_table, dna, k, nodes);

        std::cout.rdbuf(orig_buf);
    }
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
                  << "  --break               Split paths at non-ACGT characters.\n"
                  << "  -t, --threads <num>   Number of threads to use (default: number of cores).\n"
                  << "  -h, --help            Show this help message.\n";
        exit(1);
    }

    bool break_path_at_non_ACGT_character = false;

    Vec<char*> args;
    args.push(argv[0]); 

    NUM_THREADS = omp_get_max_threads(); 

    for (int i = 1; i < argc; i++) {
        std::string arg = argv[i];

        if (arg == "--break") {
            break_path_at_non_ACGT_character = true;
        } else if (arg == "--threads" || arg == "-t") {
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
    
    for (int i = 0; i < nodes.size(); i++) nodes[i].clean();
    nodes.clean();

    for (int i = 3; i < args.size(); i++) {
        std::cerr << "["<<i-2 << "/" << args.size()-3<<"] " << std::flush;
        print_path_fasta(kmer_table, args[i], k, break_path_at_non_ACGT_character);
    }

    return 0;
}
