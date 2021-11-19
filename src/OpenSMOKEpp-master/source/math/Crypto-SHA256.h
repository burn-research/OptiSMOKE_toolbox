#ifndef SHA_256_H_INCLUDED
#define SHA_256_H_INCLUDED

// This is a relatively straightforward implementation of SHA-256. It makes no particular
// attempt at optimization, instead aiming toward easy verification against the standard.
// To that end, many of the variable names are identical to those used in FIPS 180-2 and
// FIPS 180-3. 
//
// The code should be fairly portable, within a few limitations:
// 1. It requires that 'char' have 8 bits. In theory this is avoidable, but I don't think
// it's worth the bother.
// 2. It only deals with inputs in (8-bit) bytes. In theory, SHA-256 can deal with a number of 
// bits that's not a multiple of 8, but I've never needed it. Since the padding always results
// in a byte-sized stream, the only parts that would need changing would be reading and padding
// the input. The main hashing portion would be unaffected.
//
// Originally written in February 2008 for SHA-1.
// Converted to SHA-256 sometime later (sorry, I don't remember exactly when).
// 
// You can use this software any way you want to, with following limitations
// (shamelessly stolen from the Boost software license):
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE, TITLE AND NON-INFRINGEMENT. IN NO EVENT
// SHALL THE COPYRIGHT HOLDERS OR ANYONE DISTRIBUTING THE SOFTWARE BE LIABLE
// FOR ANY DAMAGES OR OTHER LIABILITY, WHETHER IN CONTRACT, TORT OR OTHERWISE,
// ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
// DEALINGS IN THE SOFTWARE.
// 
// If you put this to real use, I'd be happy to hear about it. If you find a bug, 
// I'd be interested in hearing about that too. There's even a pretty good chance 
// that I'll try to fix it, though I certainly can't guarantee that.
// 

#include <algorithm>
#include <vector>
#include <string>
#include <assert.h>
#include <iostream>
#include <sstream>
#include <iomanip>

#if defined(_MSC_VER) && _MSC_VER < 1600
typedef unsigned int uint32_t;
typedef unsigned __int64 uint64_t;
#else
#include <stdint.h>
#endif

namespace crypto 
{ 
	struct ternary_operator 
	{ 
		virtual uint32_t operator()(uint32_t x, uint32_t y, uint32_t z) = 0;
	};

class sha256 { 
    static const size_t hash_size = 8;
    static const size_t min_pad = 64;
    static const size_t block_bits = 512;
    static const size_t block_bytes = block_bits / 8;
    static const size_t block_words = block_bytes / 4;

    std::vector<uint32_t> K;
    std::vector<uint32_t> H;
    std::vector<uint32_t> W;
    std::vector<ternary_operator *> fs;
    std::vector<uint32_t> temp;

    static const size_t block_size = 16;
    static const size_t bytes_per_word = 4;
    size_t total_size;

    // hash a 512-bit block of input.
    //
    void hash_block(std::vector<uint32_t> const &block);

    // Pad the input to a multiple of 512 bits, and add the length
    // in binary to the end.
    static std::string pad(std::string const &input);

    // Turn 64 bytes into a block of 16 uint32_t's.
    std::vector<uint32_t> make_block(std::string const &in);

public:
    // Construct a SHA-256 object. More expensive that typical 
    // ctor, but not expected to be copied a lot or anything
    // like that, so it should be fairly harmless.
    sha256();

    // The two ways to provide input for hashing: as a stream or a string.
    // Either way, you get the result as a vector<uint32_t>. It's a fairly
    // small vector, so even if your compiler doesn't do return-value 
    // optimization, the time for copying isn't like to be significant.
    // 
    std::vector<uint32_t> operator()(std::string const &input);

    friend std::ostream &operator<<(std::ostream &os, sha256 const &s);
};

// A minimal test harness to check that it's working correctly. Strictly black-box
// testing, with no attempt at things like coverage analysis. Nonetheless, I believe
// it should cover most of the code -- the core hashing code all gets used for every
// possible value. The padding code should be tested fairly thoroughly as well -- the
// first test is a fairly simple case, and the second the more complex one (where the 
// padding requires adding another block).
class tester {
    bool verify(uint32_t *test_val, std::vector<uint32_t> const &hash, std::ostream &os) {
        // Verify that a result matches a test value and report result.
        for (size_t i=0; i<hash.size(); i++)
            if (hash[i] != test_val[i]) {
                os << "Mismatch. Expected: " << test_val[i] << ", but found: " << hash[i] << "\n";
                return false;
            }
            os << "Message digest Verified.\n\n";
            return true;
    }

public:

    bool operator()(uint32_t *test_val, std::string const &input) {
        std::cout << "Testing hashing from string:\n\"" << input << "\"\n";
        crypto::sha256 hasher1;
        std::vector<uint32_t> hash = hasher1(input);
        std::cout << "Message digest is:\n\t" << hasher1;
        return verify(test_val, hash, std::cerr);
    }
};
}

namespace crypto {
namespace {
    uint32_t word(int a, int b, int c, int d) {
        a &= 0xff;
        b &= 0xff;
        c &= 0xff;
        d &= 0xff;
        int val =  a << 24 | b << 16 | c << 8 | d;
        return val;
    }

    uint32_t ROTR(uint32_t number, unsigned bits) { 
        return (number >> bits) | (number << (32-bits));
    }

    uint32_t f1(uint32_t x, uint32_t y, uint32_t z) { 
        return (x & y) ^ (~x & z);
    }
    uint32_t f2(uint32_t x, uint32_t y, uint32_t z) { 
        return (x & y) ^ (x&z) ^ (y&z);
    }
    uint32_t f3(uint32_t x) { 
        return ROTR(x, 2) ^ ROTR(x, 13) ^ ROTR(x, 22);  
    }
    uint32_t f4(uint32_t x) { 
        return ROTR(x, 6) ^ ROTR(x, 11) ^ ROTR(x, 25);
    }
    uint32_t f5(uint32_t x) { 
        return ROTR(x, 7) ^ ROTR(x, 18) ^ (x >> 3);
    }
    uint32_t f6(uint32_t x) { 
        return ROTR(x, 17) ^ ROTR(x, 19) ^ (x >> 10);
    }

    uint32_t add(uint32_t a, uint32_t b) {
        return a+b;
    }
}

// Pad the input to a multiple of 512 bits, and add the length
// in binary to the end.
std::string sha256::pad(std::string const &input) {
    uint64_t length = input.size() * 8 + 1;
    size_t remainder = length % block_bits;
    size_t k = (remainder <= 448) ? 448 - remainder : 960 - remainder;

    std::string padding("\x80");
    padding.append(std::string(k/8, '\0'));
    --length;

    for (int i=sizeof(length)-1; i>-1; i--) {
        unsigned char byte = length >> (i*8) & 0xff;
        padding.push_back(byte);
    }

    std::string ret(input+padding);
    return ret;
}

// Turn 64 bytes into a vector of 16 uint32_t's.
std::vector<uint32_t> sha256::make_block(std::string const &in) { 
    assert(in.size() >= block_bytes);

    std::vector<uint32_t> ret(block_words);

    for (size_t i=0; i<block_words; i++) {
        size_t s = i*4;
        ret[i] = word(in[s], in[s+1], in[s+2], in[s+3]);
    }
    return ret;
}

sha256::sha256() : H(hash_size), W(64), temp(10) { 
    static const uint32_t H0[hash_size] = {
        0x6a09e667, 0xbb67ae85, 0x3c6ef372, 0xa54ff53a, 
        0x510e527f, 0x9b05688c, 0x1f83d9ab, 0x5be0cd19
    };

    std::copy(H0, H0+hash_size, H.begin());
}

void sha256::hash_block(std::vector<uint32_t> const &block) {
    static const uint32_t K[] = {
        0x428a2f98, 0x71374491, 0xb5c0fbcf, 0xe9b5dba5,
        0x3956c25b, 0x59f111f1, 0x923f82a4, 0xab1c5ed5, 
        0xd807aa98, 0x12835b01, 0x243185be, 0x550c7dc3,
        0x72be5d74, 0x80deb1fe, 0x9bdc06a7, 0xc19bf174,
        0xe49b69c1, 0xefbe4786, 0x0fc19dc6, 0x240ca1cc,
        0x2de92c6f, 0x4a7484aa, 0x5cb0a9dc, 0x76f988da,
        0x983e5152, 0xa831c66d, 0xb00327c8, 0xbf597fc7,
        0xc6e00bf3, 0xd5a79147, 0x06ca6351, 0x14292967,
        0x27b70a85, 0x2e1b2138, 0x4d2c6dfc, 0x53380d13,
        0x650a7354, 0x766a0abb, 0x81c2c92e, 0x92722c85,
        0xa2bfe8a1, 0xa81a664b, 0xc24b8b70, 0xc76c51a3,
        0xd192e819, 0xd6990624, 0xf40e3585, 0x106aa070, 
        0x19a4c116, 0x1e376c08, 0x2748774c, 0x34b0bcb5,
        0x391c0cb3, 0x4ed8aa4a, 0x5b9cca4f, 0x682e6ff3,
        0x748f82ee, 0x78a5636f, 0x84c87814, 0x8cc70208,
        0x90befffa, 0xa4506ceb, 0xbef9a3f7, 0xc67178f2
    };

    assert(block.size() == 16);

    std::copy(block.begin(), block.end(), W.begin());
    for (int t=16; t<64; ++t) 
        W[t] = f6(W[t-2]) + W[t-7] + f5(W[t-15]) + W[t-16];
    std::copy(H.begin(), H.end(), temp.begin());

    for (int t=0; t<64; ++t) {
        temp[8] = temp[7]+f4(temp[4]) + f1(temp[4],temp[5],temp[6])+K[t]+W[t];
        temp[9] = f3(temp[0]) + f2(temp[0], temp[1], temp[2]);
        temp[7] = temp[6];
        temp[6] = temp[5];
        temp[5] = temp[4];
        temp[4] = temp[3] + temp[8];
        temp[3] = temp[2];
        temp[2] = temp[1];
        temp[1] = temp[0];
        temp[0] = temp[8] + temp[9];
    }
    std::transform(H.begin(), H.end(), temp.begin(), H.begin(), add);
}

// Take a `std::string` as input, produce a SHA-256 hash as a vector of 16 uint32_ts'.
//
std::vector<uint32_t> sha256::operator()(std::string const &input) { 
    std::string temp(pad(input));
    std::vector<uint32_t> block(block_size);

    size_t num = temp.size()/block_bytes;

    for (unsigned block_num=0; block_num<num; block_num++) {
        size_t s;
        for (size_t i=0; i<block_size; i++) {
            s = block_num*block_bytes+i*4;
            block[i] = word(temp[s], temp[s+1], temp[s+2], temp[s+3]);
        }
        hash_block(block);  
    }
    return H;
}

std::ostream &operator<<(std::ostream &os, sha256 const &s) { 
    // Display hash result in hex.
    for (size_t i=0; i<(s.H).size(); i++)
        os << std::fixed << std::setprecision(8) << std::hex << std::setfill('0') << (s.H)[i] << " ";
    return os << std::dec << std::setfill(' ');
}
}

#endif
