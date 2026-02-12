#pragma once

#include <vector>

namespace Lloyd {
    class BitArray {
    public:
        int size;
        std::vector<char> bits;
        BitArray(int size){
            this->size = size;
            this->bits.resize((size + 7) / 8);
        }
        BitArray() = default;
        ~BitArray() = default;

        void set(int index){
            bits[index / 8] |= (1 << (index % 8));
        }

        bool get(int index){
            return (bits[index / 8] & (1 << (index % 8))) != 0;
        }

        void clear(int index){
            bits[index / 8] &= ~(1 << (index % 8));
        }

        void resize(int new_size){
            this->size = new_size;
            this->bits.resize((new_size + 7) / 8);
        }
    };
}