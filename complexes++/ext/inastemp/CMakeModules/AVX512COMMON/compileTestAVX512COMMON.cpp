
#include <x86intrin.h>
#include <xmmintrin.h> // SSE
#include <emmintrin.h> // SSE2
#include <pmmintrin.h> // SSE3
#include <tmmintrin.h> // SSSE3
#include <smmintrin.h> // SSE4

#include <immintrin.h> // AVX

int main() {
    __m512d res0d, res1d;
    res0d = _mm512_add_pd(res0d, res1d);

    __m512 res0, res1;
    res0 = _mm512_add_ps(res0, res1);

    // er
    {
        __m512d src;
        __mmask8 k;
        __m512d a;
        _mm512_mask_rcp28_pd(src, k, a);
    }
    // cd
    {
        __mmask8 k;
        __m512i tmp = _mm512_broadcastmb_epi64(k);
    }

    return 0;
}
