#include <cstddef>
#include <cstdint>

typedef int16_t daqint;
#define bswap(sample) { bswap_16(sample) }

typedef uint16_t bbsize;

typedef float bbreal;
