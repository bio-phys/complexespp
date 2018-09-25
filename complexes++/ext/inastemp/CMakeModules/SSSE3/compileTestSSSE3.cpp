
#include <x86intrin.h>
#include <xmmintrin.h> // SSE
#include <emmintrin.h> // SSE2
#include <pmmintrin.h> // SSE3

#include <tmmintrin.h>  //SSSE3

#ifndef __SSSE3__
#error __SSSE3__ must be defined
#endif

int main(){
	__m128d res0d, res1d;
	res0d = _mm_hadd_pd(res0d, res1d);

	__m128 res0, res1;
	res0 = _mm_hadd_ps(res0, res1);

	return 0;
}
