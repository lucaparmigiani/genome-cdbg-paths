#pragma once
#include<iostream>
#include<cstdint>
#include<cstdlib> //free-malloc
#include<type_traits> //is_same

#define hroundup(x) ((x) == 1 ? 1 : 1 << (64-__builtin_clzl((x)-1)))
#define hlog(x) (64-__builtin_clzl((x-1)))

#define is_set(flag, i)    (flag[i>>5] >> (i&0x1fU) & 1U)
#define set_bit(flag, i)   (flag[i>>5] |= 1U<<(i&0x1fU))
#define unset_bit(flag, i) (flag[i>>5] &= ~(1U<<(i&0x1fU)))

//const uint32_t __bitmask_true[] = { 1, 2, 4, 8, 16, 32, 64, 128 };
//const uint32_t __bitmask_false[] = {254, 253, 251, 247, 239, 223, 191, 127};
//const uint32_t __mask3bits = (1<<3) - 1;
//#define is_set(A, i) (A[i>>3] & __bitmask_true[i & __mask3bits])
//#define set_bit(A, i) (A[i>>3] |= __bitmask_true[i & __mask3bits])
//#define unset_bit(A, i) (A[i>>3] &= __bitmask_false[i & __mask3bits])

template<typename key_t>
struct HashSet {
    size_t num_buckets, num_elements, mask;
    uint32_t* used;
    key_t* keys;

    inline size_t capacity() const  { return num_buckets; }
    inline size_t size() const { return num_elements; }
    inline size_t begin() const { return 0; }
    inline size_t end() const { return num_buckets; }
    inline bool is_used(const size_t i) const { return is_set(used, i); }

    HashSet() : num_buckets(0), num_elements(0), mask(0), used(0), keys(0) {};

    static inline size_t hashx(key_t key) {
        if (std::is_same<size_t, uint32_t>::value) return hash32(key);
        else if (std::is_same<size_t,uint64_t>::value) return hash64(key);
    }
    static inline uint32_t hash32(uint32_t key) { return key * 2654435769U; }
	static inline uint64_t hash64(uint64_t key) { return key * 11400714819323198485ULL; }

    bool contains(const key_t &key) const {
        if (!num_buckets) return false;
        size_t i, first;
        i = first = hashx(key) & mask;
        while (is_set(used, i) && key != keys[i]) {
            i = (i+1) & mask;
            if (i == first) break;
        }
        return is_set(used, i) && key == keys[i];
    }

    void resize(size_t new_num_buckets) {
        if (num_buckets >= new_num_buckets) return; // No shrink support for the moment

        new_num_buckets = hroundup(new_num_buckets);
        uint32_t* new_used = (__typeof__(used)) std::calloc(new_num_buckets, sizeof(*(used))); 
        keys = (key_t*) std::realloc(keys, new_num_buckets * sizeof(key_t));

        size_t new_mask = new_num_buckets - 1;
        for (size_t i = 0; i < num_buckets; i++) {
            if (!is_set(used, i)) continue;

            key_t key = keys[i];
            unset_bit(used, i);
            while (1) {
                size_t j = hashx(key) & new_mask;
                while (is_set(new_used, j)) j = (j+1) & new_mask;
                set_bit(new_used, j);
                if (j < num_buckets && is_set(used, j)) {
                    key_t tmp = keys[j];
                    keys[j] = key;
                    key = tmp;
                    unset_bit(used, j);
                }
                else {
                    keys[j] = key;
                    break;
                }
            }
        }
        num_buckets = new_num_buckets;
        if (used) { free(used); }
        used = new_used;
        mask = new_mask;
    }

    //inline size_t const put(const key_t& key, bool *absent=0) {
    //    if (num_elements >= (num_buckets>>1) + (num_buckets>>2)) {
    //        resize(num_buckets ? num_buckets << 1 : 4);
    //    }
    //    size_t i, first;
    //    i = first = hashx(key);
    //    while (is_set(used,i) && key != keys[i]) {
    //        i = (i+1) & mask;
    //        if (i == first) break;
    //    }
    //    if (!is_set(used,i)) {
    //        keys[i] = key;
    //        num_elements++;
    //        if(absent) *absent = 1;
    //        set_bit(used,i);
    //    } else if(absent) *absent = 0;

    //    return i;
    //}

    size_t put(const key_t& key) { // return num_buckets if it was already there
        if (num_elements >= (num_buckets>>1) + (num_buckets>>2)) {
            resize(num_buckets ? num_buckets << 1 : 4);
        }
        size_t i, first;
        i = first = hashx(key) & mask;
        while (is_set(used,i) && key != keys[i]) {
            i = (i+1) & mask;
            if (i == first) break;
        }
        if (!is_set(used,i)) {
            keys[i] = key;
            num_elements++;
            set_bit(used, i);
        } else return num_buckets;

        return i;
    }

    void clean() {
        if (used) free(used);
        if (keys) free(keys);
    }
};

template<typename key_t, typename val_t>
struct HashMap {
    size_t num_buckets, num_elements, mask, log_mask;
    uint32_t* used;
    key_t* keys;
    val_t* values;

    inline size_t capacity() const  { return num_buckets; }
    inline size_t size() const { return num_elements; }
    inline size_t begin() const { return 0; }
    inline size_t end() const { return num_buckets; }
    inline bool is_used(const size_t i) const { return is_set(used, i); }

    HashMap() : num_buckets(0), num_elements(0), mask(0), log_mask(0),
                used(0), keys(0), values(0) {};

    inline size_t hashx(key_t key) const {
        if (std::is_same<size_t, uint32_t>::value) return hash32(key);
        else if (std::is_same<size_t,uint64_t>::value) return hash64(key);
    }

    inline uint32_t hash32(size_t key) const { return key * 2654435769U >> (32-log_mask); }
	inline uint64_t hash64(size_t key) const { return key * 11400714819323198485ULL>>(64-log_mask); }

    bool contains(const key_t &key) const {
        if (!num_buckets) return false;
        size_t i, first;
        i = first = hashx(key);
        while (is_set(used, i) && key != keys[i]) {
            i = (i+1) & mask;
            if (i == first) return false;
        }
        return is_set(used,i) && key == keys[i];
    }

    /* Find id of the key if present, otherwise get()==end() */
    size_t find(const key_t &key) const {
        if (!num_buckets) return num_buckets;
        size_t i, first;
        i = first = hashx(key);
        while (is_set(used, i) && key != keys[i]) {
            i = (i+1) & mask;
            if (i == first) return num_buckets;
        }
        return is_set(used,i) && key == keys[i] ? i : num_buckets;
    }

    inline void const resize(size_t new_num_buckets) {
        if (num_buckets >= new_num_buckets) return; // No shrink support for the moment

        new_num_buckets = hroundup(new_num_buckets);
        uint32_t* new_used = (uint32_t*) std::calloc(new_num_buckets, sizeof(*(used))); 
        keys = (key_t*) std::realloc(keys, new_num_buckets * sizeof(key_t));
        values = (val_t*) std::realloc(values, new_num_buckets * sizeof(val_t));

        size_t new_mask = new_num_buckets - 1;
        log_mask = hlog(new_num_buckets);
        for (size_t i = 0; i < num_buckets; i++) {
            if (!is_set(used, i)) continue;

            key_t key = keys[i];
            val_t val = values[i];
            unset_bit(used, i);
            while (1) {
                size_t j = hashx(key);
                while (is_set(new_used, j)) j = (j+1) & new_mask;
                set_bit(new_used, j);
                if (j < num_buckets && is_set(used, j)) {
                    key_t tmp_key = keys[j];
                    keys[j] = key;
                    key = tmp_key;

                    val_t tmp_val = values[j];
                    values[j] = val;
                    val = tmp_val;

                    unset_bit(used, j);
                }
                else {
                    keys[j] = key;
                    values[j] = val;
                    break;
                }
            }
        }
        num_buckets = new_num_buckets;
        if (used) { free(used); }
        used = new_used;
        mask = new_mask;
    }

    inline size_t const put(const key_t& key, bool *absent=0) {
        if (num_elements >= (num_buckets>>1) + (num_buckets>>2)) {
            resize(num_buckets ? num_buckets << 1 : 4);
        }
        size_t i, first;
        i = first = hashx(key);
        while (is_set(used,i) && key != keys[i]) {
            i = (i+1) & mask;
            if (i == first) break;
        }
        if (!is_set(used,i)) {
            keys[i] = key;
            num_elements++;
            if(absent) *absent = 1;
            set_bit(used,i);
        } else if(absent) *absent = 0;

        return i;
    }

    void clean() {
        if (used) free(used);
        if (keys) free(keys);
        if (values) free(values);
    }
};
