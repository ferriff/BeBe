#include <cstddef>
#include <cstdint>

#define bswap(sample) { bswap_16(sample) }

namespace bb {

        typedef int16_t  daqint_t;
        //typedef uint16_t bbsize;
        typedef float    real_t;

}
