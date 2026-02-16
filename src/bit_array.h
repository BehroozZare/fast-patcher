#pragma once

#include <vector>
#include <algorithm>

namespace Lloyd {

/**
 * @brief Compact bit-packed boolean array.
 *
 * Stores one bit per element using a byte-backed vector, providing O(1)
 * per-element set/get/clear operations.  Used as an active-vertex mask
 * during BFS and Lloyd iterations.
 */
class BitArray {
public:
    int size;                  ///< Number of logical bits in the array.
    std::vector<char> bits;    ///< Underlying byte storage (ceil(size/8) bytes).

    /**
     * @brief Construct a BitArray with a given number of bits.
     *
     * All bits are initially zero (cleared).
     *
     * @param size  Number of bits to allocate.
     */
    BitArray(int size){
        this->size = size;
        this->bits.resize((size + 7) / 8);
    }

    /** @brief Default constructor (empty array). */
    BitArray() = default;

    /** @brief Default destructor. */
    ~BitArray() = default;

    /**
     * @brief Set the bit at the given index to 1.
     * @param index  Zero-based bit index (must be in [0, size)).
     */
    void set(int index){
        bits[index / 8] |= (1 << (index % 8));
    }

    /**
     * @brief Query the bit at the given index.
     * @param index  Zero-based bit index (must be in [0, size)).
     * @return True if the bit is set (1), false otherwise.
     */
    bool get(int index) const {
        return (bits[index / 8] & (1 << (index % 8))) != 0;
    }

    /**
     * @brief Clear the bit at the given index to 0.
     * @param index  Zero-based bit index (must be in [0, size)).
     */
    void clear(int index){
        bits[index / 8] &= ~(1 << (index % 8));
    }

    /**
     * @brief Set every bit in the array to 1.
     *
     * @note Trailing bits beyond @c size in the last byte are also set,
     *       but they are never queried by valid indices.
     */
    void set_all(){
        std::fill(bits.begin(), bits.end(), static_cast<char>(0xFF));
    }

    /**
     * @brief Clear every bit in the array to 0.
     */
    void clear_all(){
        std::fill(bits.begin(), bits.end(), 0);
    }

    /**
     * @brief Resize the bit array to hold @p new_size bits.
     *
     * Existing bits within the new range are preserved; any newly added
     * bytes are zero-initialized by the underlying std::vector.
     *
     * @param new_size  New number of logical bits.
     */
    void resize(int new_size){
        this->size = new_size;
        this->bits.resize((new_size + 7) / 8);
    }
};

} // namespace Lloyd
