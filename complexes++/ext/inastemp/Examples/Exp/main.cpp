///////////////////////////////////////////////////////////////////////////
// Inastemp - Berenger Bramas MPCDF - 2016
// Under MIT Licence, please you must read the LICENCE file.
///////////////////////////////////////////////////////////////////////////

// In this example we time the duration to compute a given number of Exponential
// We compare the scalar version and the BestType

#include "InastempConfig.h"
#include "SCALAR/InaVecSCALARDouble.hpp"
#include "SCALAR/InaVecSCALARFloat.hpp"
#include "Common/InaTimer.hpp"

#include <memory>
#include <iostream>
#include <cmath>


#ifdef INASTEMP_USE_SSE41

#include "SSE41/InaVecSSE41Double.hpp"
#include "SSE41/InaVecSSE41Float.hpp"

inline void InaVecSSE41_exp(const float inVal[], float outVal[]) {
    __m128 vec = _mm_load_ps(inVal);
#ifdef __INTEL_COMPILER
    vec = _mm_exp_ps(vec);
#else
    const __m128 COEFF_LOG2E = _mm_set1_ps(float(InaFastExp::CoeffLog2E()));
    const __m128 COEFF_A     = _mm_set1_ps(float(InaFastExp::CoeffA32()));
    const __m128 COEFF_B     = _mm_set1_ps(float(InaFastExp::CoeffB32()));
    const __m128 COEFF_P5_A  = _mm_set1_ps(float(InaFastExp::GetCoefficient6_5()));
    const __m128 COEFF_P5_B  = _mm_set1_ps(float(InaFastExp::GetCoefficient6_4()));
    const __m128 COEFF_P5_C  = _mm_set1_ps(float(InaFastExp::GetCoefficient6_3()));
    const __m128 COEFF_P5_D  = _mm_set1_ps(float(InaFastExp::GetCoefficient6_2()));
    const __m128 COEFF_P5_E  = _mm_set1_ps(float(InaFastExp::GetCoefficient6_1()));
    const __m128 COEFF_P5_F  = _mm_set1_ps(float(InaFastExp::GetCoefficient6_0()));

    __m128 x = _mm_mul_ps(vec, COEFF_LOG2E);

    const __m128 fractional_part = _mm_sub_ps(x, _mm_floor_ps(x));

    __m128 factor = _mm_add_ps(_mm_mul_ps(_mm_add_ps(_mm_mul_ps(_mm_add_ps(
                                                                    _mm_mul_ps(_mm_add_ps(_mm_mul_ps(_mm_add_ps(_mm_mul_ps(
                                                                                                                    COEFF_P5_A, fractional_part),
                                                                                                                COEFF_P5_B),
                                                                                                     fractional_part),
                                                                                          COEFF_P5_C),
                                                                               fractional_part),
                                                                    COEFF_P5_D),
                                                                fractional_part),
                                                     COEFF_P5_E),
                                          fractional_part),
                               COEFF_P5_F);

    x = _mm_sub_ps(x, factor);

    __m128i castedInteger = _mm_cvtps_epi32(_mm_add_ps(_mm_mul_ps(COEFF_A, x), COEFF_B));

    vec                       = _mm_castsi128_ps(castedInteger);
#endif
    _mm_storeu_ps(outVal, vec);
}

inline void InaVecSSE41_exp(const double inVal[], double outVal[]) {
    __m128d vec = _mm_load_pd(inVal);
#ifdef __INTEL_COMPILER
    vec = _mm_exp_pd(vec);
#else
    const __m128d COEFF_LOG2E = _mm_set1_pd(double(InaFastExp::CoeffLog2E()));
    const __m128d COEFF_A     = _mm_set1_pd(double(InaFastExp::CoeffA64()));
    const __m128d COEFF_B     = _mm_set1_pd(double(InaFastExp::CoeffB64()));
    const __m128d COEFF_P5_X  = _mm_set1_pd(double(InaFastExp::GetCoefficient9_8()));
    const __m128d COEFF_P5_Y  = _mm_set1_pd(double(InaFastExp::GetCoefficient9_7()));
    const __m128d COEFF_P5_Z  = _mm_set1_pd(double(InaFastExp::GetCoefficient9_6()));
    const __m128d COEFF_P5_A  = _mm_set1_pd(double(InaFastExp::GetCoefficient9_5()));
    const __m128d COEFF_P5_B  = _mm_set1_pd(double(InaFastExp::GetCoefficient9_4()));
    const __m128d COEFF_P5_C  = _mm_set1_pd(double(InaFastExp::GetCoefficient9_3()));
    const __m128d COEFF_P5_D  = _mm_set1_pd(double(InaFastExp::GetCoefficient9_2()));
    const __m128d COEFF_P5_E  = _mm_set1_pd(double(InaFastExp::GetCoefficient9_1()));
    const __m128d COEFF_P5_F  = _mm_set1_pd(double(InaFastExp::GetCoefficient9_0()));

    __m128d x = _mm_mul_pd(vec, COEFF_LOG2E);

    const __m128d fractional_part = _mm_sub_pd(x, _mm_floor_pd(x));

    __m128d factor = _mm_add_pd(_mm_mul_pd(_mm_add_pd(
                                               _mm_mul_pd(_mm_add_pd(_mm_mul_pd(_mm_add_pd(
                                                                                    _mm_mul_pd(_mm_add_pd(_mm_mul_pd(_mm_add_pd(
                                                                                                                         _mm_mul_pd(_mm_add_pd(_mm_mul_pd(_mm_add_pd(_mm_mul_pd(
                                                                                                                                                                         COEFF_P5_X, fractional_part),
                                                                                                                                                                     COEFF_P5_Y),
                                                                                                                                                          fractional_part),
                                                                                                                                               COEFF_P5_Z),
                                                                                                                                    fractional_part),
                                                                                                                         COEFF_P5_A),
                                                                                                                     fractional_part),
                                                                                                          COEFF_P5_B),
                                                                                               fractional_part),
                                                                                    COEFF_P5_C),
                                                                                fractional_part),
                                                                     COEFF_P5_D),
                                                          fractional_part),
                                               COEFF_P5_E),
                                           fractional_part),
                                COEFF_P5_F);

    x = _mm_sub_pd(x, factor);

    x = _mm_add_pd(_mm_mul_pd(COEFF_A, x), COEFF_B);

    alignas(64) long int allvalint[2] = { _mm_cvtsd_si64(x),
                                          _mm_cvtsd_si64(_mm_shuffle_pd(x, x, 1)) };

    vec = _mm_castsi128_pd(_mm_set_epi64x(allvalint[1], allvalint[0]));
#endif
    _mm_storeu_pd(outVal, vec);
}

#endif

#ifdef INASTEMP_USE_AVX

#include "AVX/InaVecAVXDouble.hpp"
#include "AVX/InaVecAVXFloat.hpp"

inline void InaVecAVX_exp(const float inVal[], float outVal[]) {
    __m256 vec = _mm256_load_ps(inVal);
#ifdef __INTEL_COMPILER
    vec = _mm256_exp_ps(vec);
#else
    const __m256 COEFF_LOG2E = _mm256_set1_ps(float(InaFastExp::CoeffLog2E()));
    const __m256 COEFF_A     = _mm256_set1_ps(float(InaFastExp::CoeffA32()));
    const __m256 COEFF_B     = _mm256_set1_ps(float(InaFastExp::CoeffB32()));
    const __m256 COEFF_P5_A  = _mm256_set1_ps(float(InaFastExp::GetCoefficient6_5()));
    const __m256 COEFF_P5_B  = _mm256_set1_ps(float(InaFastExp::GetCoefficient6_4()));
    const __m256 COEFF_P5_C  = _mm256_set1_ps(float(InaFastExp::GetCoefficient6_3()));
    const __m256 COEFF_P5_D  = _mm256_set1_ps(float(InaFastExp::GetCoefficient6_2()));
    const __m256 COEFF_P5_E  = _mm256_set1_ps(float(InaFastExp::GetCoefficient6_1()));
    const __m256 COEFF_P5_F  = _mm256_set1_ps(float(InaFastExp::GetCoefficient6_0()));

    __m256 x = _mm256_mul_ps(vec, COEFF_LOG2E);

    const __m256 fractional_part = _mm256_sub_ps(x, _mm256_floor_ps(x));

    __m256 factor = _mm256_add_ps(_mm256_mul_ps(_mm256_add_ps(_mm256_mul_ps(_mm256_add_ps(
                                                                                _mm256_mul_ps(_mm256_add_ps(_mm256_mul_ps(_mm256_add_ps(_mm256_mul_ps(
                                                                                                                                            COEFF_P5_A, fractional_part),
                                                                                                                                        COEFF_P5_B),
                                                                                                                          fractional_part),
                                                                                                            COEFF_P5_C),
                                                                                              fractional_part),
                                                                                COEFF_P5_D),
                                                                            fractional_part),
                                                              COEFF_P5_E),
                                                fractional_part),
                                  COEFF_P5_F);

    x = _mm256_sub_ps(x, factor);

    __m256i castedInteger = _mm256_cvtps_epi32(_mm256_add_ps(_mm256_mul_ps(COEFF_A, x), COEFF_B));

    vec                       = _mm256_castsi256_ps(castedInteger);
#endif
    _mm256_storeu_ps(outVal, vec);
}

inline void InaVecAVX_exp(const double inVal[], double outVal[]) {
    __m256d vec = _mm256_load_pd(inVal);
#ifdef __INTEL_COMPILER
    vec = _mm256_exp_pd(vec);
#else
    const __m256d COEFF_LOG2E = _mm256_set1_pd(double(InaFastExp::CoeffLog2E()));
    const __m256d COEFF_A     = _mm256_set1_pd(double(InaFastExp::CoeffA64()));
    const __m256d COEFF_B     = _mm256_set1_pd(double(InaFastExp::CoeffB64()));
    const __m256d COEFF_P5_X  = _mm256_set1_pd(double(InaFastExp::GetCoefficient9_8()));
    const __m256d COEFF_P5_Y  = _mm256_set1_pd(double(InaFastExp::GetCoefficient9_7()));
    const __m256d COEFF_P5_Z  = _mm256_set1_pd(double(InaFastExp::GetCoefficient9_6()));
    const __m256d COEFF_P5_A  = _mm256_set1_pd(double(InaFastExp::GetCoefficient9_5()));
    const __m256d COEFF_P5_B  = _mm256_set1_pd(double(InaFastExp::GetCoefficient9_4()));
    const __m256d COEFF_P5_C  = _mm256_set1_pd(double(InaFastExp::GetCoefficient9_3()));
    const __m256d COEFF_P5_D  = _mm256_set1_pd(double(InaFastExp::GetCoefficient9_2()));
    const __m256d COEFF_P5_E  = _mm256_set1_pd(double(InaFastExp::GetCoefficient9_1()));
    const __m256d COEFF_P5_F  = _mm256_set1_pd(double(InaFastExp::GetCoefficient9_0()));

    __m256d x = _mm256_mul_pd(vec, COEFF_LOG2E);

    const __m256d fractional_part = _mm256_sub_pd(x, _mm256_floor_pd(x));

    __m256d factor = _mm256_add_pd(_mm256_mul_pd(_mm256_add_pd(
                                                     _mm256_mul_pd(_mm256_add_pd(_mm256_mul_pd(_mm256_add_pd(
                                                                                                   _mm256_mul_pd(_mm256_add_pd(_mm256_mul_pd(_mm256_add_pd(
                                                                                                                                                 _mm256_mul_pd(_mm256_add_pd(_mm256_mul_pd(_mm256_add_pd(_mm256_mul_pd(
                                                                                                                                                                                                             COEFF_P5_X, fractional_part),
                                                                                                                                                                                                         COEFF_P5_Y),
                                                                                                                                                                                           fractional_part),
                                                                                                                                                                             COEFF_P5_Z),
                                                                                                                                                               fractional_part),
                                                                                                                                                 COEFF_P5_A),
                                                                                                                                             fractional_part),
                                                                                                                               COEFF_P5_B),
                                                                                                                 fractional_part),
                                                                                                   COEFF_P5_C),
                                                                                               fractional_part),
                                                                                 COEFF_P5_D),
                                                                   fractional_part),
                                                     COEFF_P5_E),
                                                 fractional_part),
                                   COEFF_P5_F);

    x = _mm256_sub_pd(x, factor);

    x = _mm256_add_pd(_mm256_mul_pd(COEFF_A, x), COEFF_B);

    __m128d valupper = _mm256_extractf128_pd(x, 1);
    __m128d vallower = _mm256_castpd256_pd128(x);
    // Removed because it makes GCC bugging:_mm256_zeroupper();

    alignas(64) long long int allvalint[4] = { _mm_cvtsd_si64(vallower),
                                               _mm_cvtsd_si64(_mm_shuffle_pd(vallower, vallower, 1)),
                                               _mm_cvtsd_si64(valupper),
                                               _mm_cvtsd_si64(_mm_shuffle_pd(valupper, valupper, 1)) };

    vec = _mm256_castsi256_pd(_mm256_load_si256(reinterpret_cast< const __m256i* >(allvalint)));
#endif
    _mm256_storeu_pd(outVal, vec);
}


#endif

#ifdef INASTEMP_USE_AVX512KNL

#include "AVX512KNL/InaVecAVX512KNLDouble.hpp"
#include "AVX512KNL/InaVecAVX512KNLFloat.hpp"

inline void InaVecAVX512KNL_exp(const float inVal[], float outVal[]) {
    __m512 vec = _mm512_load_ps(inVal);

    const __m512 COEFF_LOG2E = _mm512_set1_ps(float(InaFastExp::CoeffLog2E()));
    const __m512 COEFF_A     = _mm512_set1_ps(float(InaFastExp::CoeffA32()));
    const __m512 COEFF_B     = _mm512_set1_ps(float(InaFastExp::CoeffB32()));
    const __m512 COEFF_P5_A  = _mm512_set1_ps(float(InaFastExp::GetCoefficient6_5()));
    const __m512 COEFF_P5_B  = _mm512_set1_ps(float(InaFastExp::GetCoefficient6_4()));
    const __m512 COEFF_P5_C  = _mm512_set1_ps(float(InaFastExp::GetCoefficient6_3()));
    const __m512 COEFF_P5_D  = _mm512_set1_ps(float(InaFastExp::GetCoefficient6_2()));
    const __m512 COEFF_P5_E  = _mm512_set1_ps(float(InaFastExp::GetCoefficient6_1()));
    const __m512 COEFF_P5_F  = _mm512_set1_ps(float(InaFastExp::GetCoefficient6_0()));

    __m512 x = _mm512_mul_ps(vec, COEFF_LOG2E);

    const __m512 fractional_part = _mm512_sub_ps(x, _mm512_floor_ps(x));

    __m512 factor = _mm512_add_ps(_mm512_mul_ps(_mm512_add_ps(_mm512_mul_ps(_mm512_add_ps(
                                                                                _mm512_mul_ps(_mm512_add_ps(_mm512_mul_ps(_mm512_add_ps(_mm512_mul_ps(
                                                                                                                                            COEFF_P5_A, fractional_part),
                                                                                                                                        COEFF_P5_B),
                                                                                                                          fractional_part),
                                                                                                            COEFF_P5_C),
                                                                                              fractional_part),
                                                                                COEFF_P5_D),
                                                                            fractional_part),
                                                              COEFF_P5_E),
                                                fractional_part),
                                  COEFF_P5_F);

    x = _mm512_sub_ps(x, factor);

    __m512i castedInteger = _mm512_cvtps_epi32(_mm512_add_ps(_mm512_mul_ps(COEFF_A, x), COEFF_B));

    vec = _mm512_castsi512_ps(castedInteger);
    _mm512_storeu_ps(outVal, vec);
}

inline void InaVecAVX512KNL_exp(const double inVal[], double outVal[]) {
    __m512d vec = _mm512_load_pd(inVal);

    const __m512d COEFF_LOG2E = _mm512_set1_pd(double(InaFastExp::CoeffLog2E()));
    const __m512d COEFF_A     = _mm512_set1_pd(double(InaFastExp::CoeffA64()));
    const __m512d COEFF_B     = _mm512_set1_pd(double(InaFastExp::CoeffB64()));
    const __m512d COEFF_P5_X  = _mm512_set1_pd(double(InaFastExp::GetCoefficient9_8()));
    const __m512d COEFF_P5_Y  = _mm512_set1_pd(double(InaFastExp::GetCoefficient9_7()));
    const __m512d COEFF_P5_Z  = _mm512_set1_pd(double(InaFastExp::GetCoefficient9_6()));
    const __m512d COEFF_P5_A  = _mm512_set1_pd(double(InaFastExp::GetCoefficient9_5()));
    const __m512d COEFF_P5_B  = _mm512_set1_pd(double(InaFastExp::GetCoefficient9_4()));
    const __m512d COEFF_P5_C  = _mm512_set1_pd(double(InaFastExp::GetCoefficient9_3()));
    const __m512d COEFF_P5_D  = _mm512_set1_pd(double(InaFastExp::GetCoefficient9_2()));
    const __m512d COEFF_P5_E  = _mm512_set1_pd(double(InaFastExp::GetCoefficient9_1()));
    const __m512d COEFF_P5_F  = _mm512_set1_pd(double(InaFastExp::GetCoefficient9_0()));

    __m512d x = _mm512_mul_pd(vec, COEFF_LOG2E);

    const __m512d fractional_part = _mm512_sub_pd(x, _mm512_floor_pd(x));

    __m512d factor = _mm512_add_pd(_mm512_mul_pd(_mm512_add_pd(
                                                     _mm512_mul_pd(_mm512_add_pd(_mm512_mul_pd(_mm512_add_pd(
                                                                                                   _mm512_mul_pd(_mm512_add_pd(_mm512_mul_pd(_mm512_add_pd(
                                                                                                                                                 _mm512_mul_pd(_mm512_add_pd(_mm512_mul_pd(_mm512_add_pd(_mm512_mul_pd(
                                                                                                                                                                                                             COEFF_P5_X, fractional_part),
                                                                                                                                                                                                         COEFF_P5_Y),
                                                                                                                                                                                           fractional_part),
                                                                                                                                                                             COEFF_P5_Z),
                                                                                                                                                               fractional_part),
                                                                                                                                                 COEFF_P5_A),
                                                                                                                                             fractional_part),
                                                                                                                               COEFF_P5_B),
                                                                                                                 fractional_part),
                                                                                                   COEFF_P5_C),
                                                                                               fractional_part),
                                                                                 COEFF_P5_D),
                                                                   fractional_part),
                                                     COEFF_P5_E),
                                                 fractional_part),
                                   COEFF_P5_F);

    x = _mm512_sub_pd(x, factor);

    x = _mm512_add_pd(_mm512_mul_pd(COEFF_A, x), COEFF_B);

    alignas(64) double allvalreal[8];
    _mm512_store_pd(allvalreal, x);

    alignas(64) long long int allvalint[8] = { static_cast< long long int >(allvalreal[0]), static_cast< long long int >(allvalreal[1]),
                                               static_cast< long long int >(allvalreal[2]), static_cast< long long int >(allvalreal[3]),
                                               static_cast< long long int >(allvalreal[4]), static_cast< long long int >(allvalreal[5]),
                                               static_cast< long long int >(allvalreal[6]), static_cast< long long int >(allvalreal[7]) };

    vec = _mm512_castsi512_pd(_mm512_load_epi64(reinterpret_cast< const __m512i* >(allvalint)));
    _mm512_storeu_pd(outVal, vec);
}
#endif


#ifdef INASTEMP_USE_ALTIVEC

#include "ALTIVEC/InaVecALTIVECDouble.hpp"
#include "ALTIVEC/InaVecALTIVECFloat.hpp"

inline void InaVecALTIVEC_exp(const float inVal[], float outVal[]) {
    __vector float vec = vec_xl(0, const_cast< float* >(inVal));

    const __vector float COEFF_LOG2E = vec_splats(float(InaFastExp::CoeffLog2E()));
    const __vector float COEFF_A     = vec_splats(float(InaFastExp::CoeffA32()));
    const __vector float COEFF_B     = vec_splats(float(InaFastExp::CoeffB32()));
    const __vector float COEFF_P5_A  = vec_splats(float(InaFastExp::GetCoefficient6_5()));
    const __vector float COEFF_P5_B  = vec_splats(float(InaFastExp::GetCoefficient6_4()));
    const __vector float COEFF_P5_C  = vec_splats(float(InaFastExp::GetCoefficient6_3()));
    const __vector float COEFF_P5_D  = vec_splats(float(InaFastExp::GetCoefficient6_2()));
    const __vector float COEFF_P5_E  = vec_splats(float(InaFastExp::GetCoefficient6_1()));
    const __vector float COEFF_P5_F  = vec_splats(float(InaFastExp::GetCoefficient6_0()));

    __vector float x = vec * COEFF_LOG2E;

    const __vector float fractional_part = x - vec_floor(x);

    __vector float factor = (((((COEFF_P5_A * fractional_part + COEFF_P5_B) * fractional_part + COEFF_P5_C) * fractional_part + COEFF_P5_D) * fractional_part + COEFF_P5_E) * fractional_part + COEFF_P5_F);

    x -= factor;

    __vector int castedInteger = vec_cts(COEFF_A * x + COEFF_B, 0);

    vec = reinterpret_cast< __vector float >(castedInteger);
    vec_xst(vec, 0, outVal);
}

inline void InaVecALTIVEC_exp(const double inVal[], double outVal[]) {
    __vector double vec = vec_xl(0, const_cast< double* >(&inVal[0]));

    const __vector double COEFF_LOG2E = vec_splats(double(InaFastExp::CoeffLog2E()));
    const __vector double COEFF_A     = vec_splats(double(InaFastExp::CoeffA64()));
    const __vector double COEFF_B     = vec_splats(double(InaFastExp::CoeffB64()));
    const __vector double COEFF_P5_C  = vec_splats(double(InaFastExp::GetCoefficient4_3()));
    const __vector double COEFF_P5_D  = vec_splats(double(InaFastExp::GetCoefficient4_2()));
    const __vector double COEFF_P5_E  = vec_splats(double(InaFastExp::GetCoefficient4_1()));
    const __vector double COEFF_P5_F  = vec_splats(double(InaFastExp::GetCoefficient4_0()));

    __vector double x = vec * COEFF_LOG2E;

    const __vector double fractional_part = x - vec_floor(x);

    __vector double factor = (((COEFF_P5_C * fractional_part + COEFF_P5_D) * fractional_part + COEFF_P5_E) * fractional_part + COEFF_P5_F);

    x -= factor;

    x = COEFF_A * x + COEFF_B;

    // TODO find conversion function
    //__vector long castedInteger = vec_cts(x, 0);
    //return reinterpret_cast<__vector double>(castedInteger);
    alignas(16) double tmpptr[2];
    vec_st(reinterpret_cast< __vector unsigned int >(x), 0, reinterpret_cast< unsigned int* >(tmpptr));

#ifdef INASTEMP_USE_XL
    alignas(16) long long ltmpptr[2];
#else
    alignas(16) long ltmpptr[2];
#endif
    ltmpptr[0] = long(tmpptr[0]);
    ltmpptr[1] = long(tmpptr[1]);
    vec        = reinterpret_cast< __vector double >(vec_xl(0, ltmpptr));
    vec_xst(reinterpret_cast< __vector unsigned int >(vec), 0, reinterpret_cast< unsigned int* >(outVal));
}
#endif

template < class VecType >
void GenericExpInavec(const size_t NbOverLoop, const size_t NbExp) {
    using RealType = typename VecType::RealType;
    // Note : we increase the length of the vector to avoid checking the loop size
    std::unique_ptr< RealType[] > resIna(new RealType[NbExp + VecType::VecLength]);
    InaTimer timer;

    for(size_t idxLoop = 0; idxLoop < NbOverLoop; ++idxLoop) {
        for(size_t idx = 0; idx < NbExp; idx += VecType::VecLength) {
            alignas(64) RealType bufferX[VecType::VecLength];
            // Copy value into a buffer since we do it on the fly
            for(size_t idxX = 0; idxX < VecType::VecLength; ++idxX) {
                bufferX[idxX] = static_cast< RealType >((idx + idxX) % 200);
            }
            VecType().setFromAlignedArray(bufferX).exp().storeInArray(&resIna[idx]);
        }
    }

    timer.stop();
    std::cout << "Vector " << VecType::GetName() << " for " << NbExp * NbOverLoop
              << " exp took " << timer.getElapsed() << "s (" << timer.getElapsed() / double(NbExp * NbOverLoop) << "s per exp)\n";
}


template < class RealType >
void compareExpTime(const size_t NbOverLoop, const size_t NbExp) {

    /////////////////////////////////////////////////////////////

    std::unique_ptr< RealType[] > resScalar(new RealType[NbExp]);
    {
        InaTimer timer;

        for(size_t idxLoop = 0; idxLoop < NbOverLoop; ++idxLoop) {
            for(size_t idx = 0; idx < NbExp; ++idx) {
                resScalar[idx] = static_cast< RealType >(std::exp(RealType(idx % 200)));
            }
        }

        timer.stop();
        std::cout << "Scalar for " << NbExp * NbOverLoop
                  << " exp took " << timer.getElapsed() << "s (" << timer.getElapsed() / double(NbExp * NbOverLoop) << "s per exp)\n";
        // Ensure that optimization compute for real
        volatile RealType tmp;
        tmp = resScalar[0];
    }
    std::cout << "\n";
    /////////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////

    {
        InaTimer timer;

        for(size_t idxLoop = 0; idxLoop < NbOverLoop; ++idxLoop) {
            for(size_t idx = 0; idx < NbExp; ++idx) {
                resScalar[idx] = static_cast< RealType >(InaVecSCALAR< RealType >(RealType(idx % 200)).exp());
            }
        }

        timer.stop();
        std::cout << "Inastemp Scalar for " << NbExp * NbOverLoop
                  << " exp took " << timer.getElapsed() << "s (" << timer.getElapsed() / double(NbExp * NbOverLoop) << "s per exp)\n";
        // Ensure that optimization compute for real
        volatile RealType tmp;
        tmp = resScalar[0];
    }
    std::cout << "\n";
/////////////////////////////////////////////////////////////

#ifdef INASTEMP_USE_SSE41
    GenericExpInavec< InaVecSSE41< RealType > >(NbOverLoop, NbExp);

    {
        // Raw SIMD
        const int VecLength = InaVecSSE41< RealType >::VecLength;
        // Note : we increase the length of the vector to avoid checking the loop size
        std::unique_ptr< RealType[] > resSimd(new RealType[NbExp + VecLength]);
        InaTimer timer;

        for(size_t idxLoop = 0; idxLoop < NbOverLoop; ++idxLoop) {
            for(size_t idx = 0; idx < NbExp; idx += VecLength) {
                alignas(64) RealType bufferX[VecLength];
                // Copy value into a buffer since we do it on the fly
                for(size_t idxX = 0; idxX < VecLength; ++idxX) {
                    bufferX[idxX] = static_cast< RealType >((idx + idxX) % 200);
                }
                InaVecSSE41_exp(bufferX, &resSimd[idx]);
            }
        }

        timer.stop();
        std::cout << "Vector "
                  << "SSE41"
                  << " for " << NbExp * NbOverLoop
                  << " exp took " << timer.getElapsed() << "s (" << timer.getElapsed() / double(NbExp * NbOverLoop) << "s per exp)\n";
    }
    std::cout << "\n";
#endif
/////////////////////////////////////////////////////////////

#ifdef INASTEMP_USE_AVX
    GenericExpInavec< InaVecAVX< RealType > >(NbOverLoop, NbExp);

    {
        // Raw SIMD
        const int VecLength = InaVecAVX< RealType >::VecLength;
        // Note : we increase the length of the vector to avoid checking the loop size
        std::unique_ptr< RealType[] > resSimd(new RealType[NbExp + VecLength]);
        InaTimer timer;

        for(size_t idxLoop = 0; idxLoop < NbOverLoop; ++idxLoop) {
            for(size_t idx = 0; idx < NbExp; idx += VecLength) {
                alignas(64) RealType bufferX[VecLength];
                // Copy value into a buffer since we do it on the fly
                for(size_t idxX = 0; idxX < VecLength; ++idxX) {
                    bufferX[idxX] = static_cast< RealType >((idx + idxX) % 200);
                }
                InaVecAVX_exp(bufferX, &resSimd[idx]);
            }
        }

        timer.stop();
        std::cout << "Vector "
                  << "AVX"
                  << " for " << NbExp * NbOverLoop
                  << " exp took " << timer.getElapsed() << "s (" << timer.getElapsed() / double(NbExp * NbOverLoop) << "s per exp)\n";
    }
    std::cout << "\n";
#endif
/////////////////////////////////////////////////////////////

#ifdef INASTEMP_USE_AVX512KNL
    GenericExpInavec< InaVecAVX512KNL< RealType > >(NbOverLoop, NbExp);

    {
        // Raw SIMD
        const int VecLength = InaVecAVX512KNL< RealType >::VecLength;
        // Note : we increase the length of the vector to avoid checking the loop size
        std::unique_ptr< RealType[] > resSimd(new RealType[NbExp + VecLength]);
        InaTimer timer;

        for(size_t idxLoop = 0; idxLoop < NbOverLoop; ++idxLoop) {
            for(size_t idx = 0; idx < NbExp; idx += VecLength) {
                alignas(64) RealType bufferX[VecLength];
                // Copy value into a buffer since we do it on the fly
                for(size_t idxX = 0; idxX < VecLength; ++idxX) {
                    bufferX[idxX] = static_cast< RealType >((idx + idxX) % 200);
                }
                InaVecAVX512KNL_exp(bufferX, &resSimd[idx]);
            }
        }

        timer.stop();
        std::cout << "Vector "
                  << "AVX512KNL"
                  << " for " << NbExp * NbOverLoop
                  << " exp took " << timer.getElapsed() << "s (" << timer.getElapsed() / double(NbExp * NbOverLoop) << "s per exp)\n";
    }
    std::cout << "\n";
#endif

#ifdef INASTEMP_USE_ALTIVEC
    GenericExpInavec< InaVecALTIVEC< RealType > >(NbOverLoop, NbExp);

    {
        // Raw SIMD
        const int VecLength = InaVecALTIVEC< RealType >::VecLength;
        // Note : we increase the length of the vector to avoid checking the loop size
        std::unique_ptr< RealType[] > resSimd(new RealType[NbExp + VecLength]);
        InaTimer timer;

        for(size_t idxLoop = 0; idxLoop < NbOverLoop; ++idxLoop) {
            for(size_t idx = 0; idx < NbExp; idx += VecLength) {
                alignas(64) RealType bufferX[VecLength];
                // Copy value into a buffer since we do it on the fly
                for(size_t idxX = 0; idxX < VecLength; ++idxX) {
                    bufferX[idxX] = static_cast< RealType >((idx + idxX) % 200);
                }
                InaVecALTIVEC_exp(bufferX, &resSimd[idx]);
            }
        }

        timer.stop();
        std::cout << "Vector "
                  << "ALTIVEC"
                  << " for " << NbExp * NbOverLoop
                  << " exp took " << timer.getElapsed() << "s (" << timer.getElapsed() / double(NbExp * NbOverLoop) << "s per exp)\n";
    }
    std::cout << "\n";
#endif
}

int main(int /*argc*/, char* /*argv*/ []) {
    std::cout << "[INFO] This program runs the computation of exp() using scalar, intrinsic vectors or inastemp vectors. \n";

    const size_t NbOverLoop = 7;
    const size_t NbExp      = 10240000;
    std::cout << "[INFO] It will compute " << NbExp << " consecutive exp, and store them in an array. \n";
    std::cout << "[INFO] This process will be done " << NbOverLoop << " times. \n";

    std::cout << "[INFO] In Float:" << std::endl;
    compareExpTime< float >(NbOverLoop, NbExp);

    std::cout << "[INFO] In Double:" << std::endl;
    compareExpTime< double >(NbOverLoop, NbExp);

    return 0;
}
