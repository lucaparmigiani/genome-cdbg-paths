#pragma once
#include <algorithm>
#include <iostream>
#include <stdio.h>
#include <cstring>
#include <fstream>
#include <omp.h>
#include "vec.hpp"
#include "hash.hpp"
#include "str.hpp"

#define kmer_t          uint64_t
#define value_t         uint64_t

int NUM_THREADS = 1;

uint64_t lookup_rc[256] = {
    0xff,0xbf,0x7f,0x3f, 0xef,0xaf,0x6f,0x2f, 0xdf,0x9f,0x5f,0x1f, 0xcf,0x8f,0x4f,0x0f,
    0xfb,0xbb,0x7b,0x3b, 0xeb,0xab,0x6b,0x2b, 0xdb,0x9b,0x5b,0x1b, 0xcb,0x8b,0x4b,0x0b,
    0xf7,0xb7,0x77,0x37, 0xe7,0xa7,0x67,0x27, 0xd7,0x97,0x57,0x17, 0xc7,0x87,0x47,0x07,
    0xf3,0xb3,0x73,0x33, 0xe3,0xa3,0x63,0x23, 0xd3,0x93,0x53,0x13, 0xc3,0x83,0x43,0x03,
    0xfe,0xbe,0x7e,0x3e, 0xee,0xae,0x6e,0x2e, 0xde,0x9e,0x5e,0x1e, 0xce,0x8e,0x4e,0x0e,
    0xfa,0xba,0x7a,0x3a, 0xea,0xaa,0x6a,0x2a, 0xda,0x9a,0x5a,0x1a, 0xca,0x8a,0x4a,0x0a,
    0xf6,0xb6,0x76,0x36, 0xe6,0xa6,0x66,0x26, 0xd6,0x96,0x56,0x16, 0xc6,0x86,0x46,0x06,
    0xf2,0xb2,0x72,0x32, 0xe2,0xa2,0x62,0x22, 0xd2,0x92,0x52,0x12, 0xc2,0x82,0x42,0x02,
    0xfd,0xbd,0x7d,0x3d, 0xed,0xad,0x6d,0x2d, 0xdd,0x9d,0x5d,0x1d, 0xcd,0x8d,0x4d,0x0d,
    0xf9,0xb9,0x79,0x39, 0xe9,0xa9,0x69,0x29, 0xd9,0x99,0x59,0x19, 0xc9,0x89,0x49,0x09,
    0xf5,0xb5,0x75,0x35, 0xe5,0xa5,0x65,0x25, 0xd5,0x95,0x55,0x15, 0xc5,0x85,0x45,0x05,
    0xf1,0xb1,0x71,0x31, 0xe1,0xa1,0x61,0x21, 0xd1,0x91,0x51,0x11, 0xc1,0x81,0x41,0x01,
    0xfc,0xbc,0x7c,0x3c, 0xec,0xac,0x6c,0x2c, 0xdc,0x9c,0x5c,0x1c, 0xcc,0x8c,0x4c,0x0c,
    0xf8,0xb8,0x78,0x38, 0xe8,0xa8,0x68,0x28, 0xd8,0x98,0x58,0x18, 0xc8,0x88,0x48,0x08,
    0xf4,0xb4,0x74,0x34, 0xe4,0xa4,0x64,0x24, 0xd4,0x94,0x54,0x14, 0xc4,0x84,0x44,0x04,
    0xf0,0xb0,0x70,0x30, 0xe0,0xa0,0x60,0x20, 0xd0,0x90,0x50,0x10, 0xc0,0x80,0x40,0x00
};

uint64_t revcmp(uint64_t kmer, int k) {
    return (lookup_rc[ (kmer)       & 0xffULL ]<<56 |
            lookup_rc[ (kmer >>8)   & 0xffULL ]<<48 | 
            lookup_rc[ (kmer >>16 ) & 0xffULL ]<<40 |
            lookup_rc[ (kmer >>24 ) & 0xffULL ]<<32 | 
            lookup_rc[ (kmer >>32 ) & 0xffULL ]<<24 | 
            lookup_rc[ (kmer >>40 ) & 0xffULL ]<<16 | 
            lookup_rc[ (kmer >>48 ) & 0xffULL ]<<8  | 
            lookup_rc[ (kmer >>56 ) & 0xffULL ] ) >> (64 - k*2);
}

unsigned char nt_2_bits[256] = { // translate ACGT to 0123
    0,1,2,3, 4,4,4,4, 4,4,4,4, 4,4,4,4,
    4,4,4,4, 4,4,4,4, 4,4,4,4, 4,4,4,4,
    4,4,4,4, 4,4,4,4, 4,4,4,4, 4,4,4,4,
    4,4,4,4, 4,4,4,4, 4,4,4,4, 4,4,4,4,
    4,0,4,1, 4,4,4,2, 4,4,4,4, 4,4,4,4,
    4,4,4,4, 3,3,4,4, 4,4,4,4, 4,4,4,4,
    4,0,4,1, 4,4,4,2, 4,4,4,4, 4,4,4,4,
    4,4,4,4, 3,3,4,4, 4,4,4,4, 4,4,4,4,
    4,4,4,4, 4,4,4,4, 4,4,4,4, 4,4,4,4,
    4,4,4,4, 4,4,4,4, 4,4,4,4, 4,4,4,4,
    4,4,4,4, 4,4,4,4, 4,4,4,4, 4,4,4,4,
    4,4,4,4, 4,4,4,4, 4,4,4,4, 4,4,4,4,
    4,4,4,4, 4,4,4,4, 4,4,4,4, 4,4,4,4,
    4,4,4,4, 4,4,4,4, 4,4,4,4, 4,4,4,4,
    4,4,4,4, 4,4,4,4, 4,4,4,4, 4,4,4,4,
    4,4,4,4, 4,4,4,4, 4,4,4,4, 4,4,4,4
};

uint64_t kmer2bits(const char *kmer_char, int k, bool canonical=true) {
    kmer_t kmer[2] = {0,0};
    kmer_t shift = (k - 1) * 2;

    for (int i = 0; i < k; i++) {
        unsigned char c = nt_2_bits[kmer_char[i]];
        if (c < 4) {
            kmer[0] = (kmer[0] << 2 | c);
            kmer[1] = kmer[1] >> 2 | (uint64_t)(3 - c) << shift;  
        }
        else {
            std::cout << kmer[i] << " is not a nucleotide base" << '\n' << std::flush;
            return SIZE_MAX;
        }
    }
    return canonical ? (kmer[0] < kmer[1] ? kmer[0] : kmer[1]) : kmer[0];
}

std::string bits2kmer(uint64_t kmer_bits, int k) {
    char nucleotides[4] = {'A','C','G','T'};
    char kmer_char[k+1];
    for (int i = 0; i < k; i++) {
        kmer_char[i] = nucleotides[(kmer_bits>>(2*k-(i+1)*2)) & 3];
    }
    kmer_char[k] = '\0';
    return kmer_char;
}

typedef struct {
    kmer_t kmer;
    value_t value;
} kmer_value_t;

typedef struct {
    uint32_t rid, pos; /* position corresponds to the last nt */
    bool strand;
} res_align_t;

typedef struct {
    HashMap<kmer_t, value_t> h;
} hashmap_t;

typedef struct {
    uint64_t k, num_kmers, suf, mask_suf ;
    hashmap_t *hm;

    res_align_t align(kmer_t kmer, bool canonical=true) {
        int strand = 0;
        if (canonical) {
            kmer_t kmer_rc = revcmp(kmer, k);
            strand = kmer > kmer_rc;
            kmer = kmer < kmer_rc ? kmer : kmer_rc;
        }
        kmer_t kmer_suf = kmer & mask_suf;
        size_t idx = hm[kmer_suf].h.find(kmer);
        if (idx != hm[kmer_suf].h.end()) {
            value_t val = hm[kmer_suf].h.values[idx];
            return {
                (uint32_t) (val >> 31), 
                (uint32_t) ((val & ((1ULL << 31) - 1))>>1),
                (bool)     ((val & 1ULL) ^ strand)
            };
        }
        return {UINT32_MAX, UINT32_MAX, 0};
    }

    /* print table for debugging */
    void print() {
        for (int i = 0; i < 1<<suf; i++) {
            if (hm[i].h.num_elements) {
                std::cout << i << ": { " << std::flush;
                for (int j = 0; j < hm[i].h.end(); j++) {
                    if (hm[i].h.is_used(j)) {
                        std::cout << hm[i].h.keys[j]   << ":"
                                  << (hm[i].h.values[j]>>32) << ' ' << std::flush;
                    }
                }
                std::cout << "}" << '\n' << std::flush;
            }
        }
    }

    /* returns SIZE_MAX in case it is not there */
    size_t find(kmer_t kmer, bool canonical=true) {
        if (canonical) {
            kmer_t kmer_rc = revcmp(kmer, k);
            kmer = kmer < kmer_rc ? kmer : kmer_rc;
        }
        kmer_t kmer_suf = kmer & mask_suf;
        size_t idx = hm[kmer_suf].h.find(kmer);
        if (idx != hm[kmer_suf].h.end()) return idx;
        return SIZE_MAX;
    }
} kmertable_t;


void insert_to_hashtable(kmertable_t*& kmer_table, int size_suf, uint64_t* n_ins, Vec<kmer_value_t>* kmer_buffer, int k) {
    std::cerr << "Adding kmers to hashtable..." << std::flush;
    #pragma omp parallel for num_threads(NUM_THREADS)
    for (int i = 0; i < size_suf; i++) {
        for (int j = 0; j < kmer_buffer[i].size(); j++) {
            bool absent;
            size_t idx = kmer_table->hm[i].h.put(kmer_buffer[i][j].kmer, &absent);
            if (absent) n_ins[omp_get_thread_num()]++;
            else std::cerr << "kmer " << bits2kmer(kmer_buffer[i][j].kmer, k) << " appears twice\n" << std::flush;

            kmer_table->hm[i].h.values[idx] = kmer_buffer[i][j].value;
        }
        kmer_buffer[i].clean();
    }
    for (int i = 0; i < NUM_THREADS; i++) kmer_table->num_kmers += n_ins[i];
    std::cerr << "ok\n" << std::flush;
}

//kmertable_t* count_kmers(Vec<std::string> dna, int k) {
kmertable_t* count_kmers(Vec<str> &dna, int k) {
    kmer_t kmer[2] = {0,0};
    int len = 0;
    int suf = 10;
    int size_suf = 1<<suf;
    int added_kmers = 0;
    Vec<kmer_value_t> kmer_buffer[size_suf];

    kmer_t mask = (1ULL << (2*k)) - 1;
    kmer_t shift = (k - 1) * 2;
    kmer_t mask_suf = (1 << suf) - 1;
    
    //initialize kmer_table
    kmertable_t *kmer_table = (kmertable_t*) calloc(1, sizeof(kmertable_t));
    kmer_table->k = k;
    kmer_table->suf = suf;
    kmer_table->mask_suf = mask_suf;
    kmer_table->hm = (hashmap_t*) calloc(1 << kmer_table->suf, sizeof(hashmap_t));
    for (int i = 0; i < 1<<kmer_table->suf; ++i)
        kmer_table->hm[i].h = HashMap<kmer_t, value_t>(); 

    uint64_t n_ins[NUM_THREADS];
    memset(n_ins, 0, sizeof(n_ins));

    int dbc = 0;
    for (uint64_t i = 0; i < dna.size(); i++) {
        for (uint64_t j = 0; j < dna[i].length(); j++) {
            unsigned char c = nt_2_bits[dna[i][j]];
            if (c < 4) {
                kmer[0] = (kmer[0] << 2 | c) & mask;
                kmer[1] = kmer[1] >> 2 | (uint64_t)(3 - c) << shift;  
                if (len+1 < k) len++;
                else {
                    dbc++;
                    uint64_t strand = kmer[0] < kmer[1] ? 0 : 1;
                    kmer_buffer[kmer[strand] & mask_suf].push({kmer[strand], (i<<31) | (j<<1) | strand });
                    if(++added_kmers > 1024*1024*1024) {
                        insert_to_hashtable(kmer_table, size_suf, n_ins, kmer_buffer, k);
                        added_kmers = 0;
                    }
                }
            } else len = 0, kmer[0] = kmer[1] = 0;
        }
        len = 0, kmer[0] = kmer[1] = 0;
    }
    if (added_kmers>0) insert_to_hashtable(kmer_table, size_suf, n_ins, kmer_buffer, k);
   
    return kmer_table;
}

/*
int main(int argc, char *argv[])
{
    Vec<std::string> read;
    read.push("ACGTAGCTACTGATCG");
    read.push("ACGTAGCTACTGATGC");
    kmertable_t *kmer_table = count_kmers(read, 11);
    std::cout << kmer_table->num_kmers << '\n' << std::flush;
    return 0;
}
*/
