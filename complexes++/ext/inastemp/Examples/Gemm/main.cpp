///////////////////////////////////////////////////////////////////////////
// Inastemp - Berenger Bramas MPCDF - 2016
// Under MIT Licence, please you must read the LICENCE file.
///////////////////////////////////////////////////////////////////////////

// In this example we time the duration to compute a gemm
// to simplify the example, the matrices are square
// and the size is a factor of the blocking coefficient

#include "InastempConfig.h"
#include "Common/InaTimer.hpp"


#include "SCALAR/InaVecSCALARDouble.hpp"
#include "SCALAR/InaVecSCALARFloat.hpp"

#include <memory>
#include <iostream>
#include <cstring>
#include <cassert>

///////////////////////////////////////////////////////////////////////////////////
///
/// Utils
///
///////////////////////////////////////////////////////////////////////////////////

template < class RealType, size_t BlockSize >
void CopyMat(RealType* __restrict__ dest, const size_t leadingDest,
             const RealType* __restrict__ src, const size_t leadingSrc) {

    for(size_t idxl = 0; idxl < leadingDest; ++idxl) {
        for(size_t idxb = 0; idxb < BlockSize; ++idxb) {
            dest[idxb * leadingDest + idxl] = src[idxb * leadingSrc + idxl];
        }
    }
}

template < class RealType, size_t BlockSize >
void CopyMatT(RealType* __restrict__ dest, const size_t leadingDest, const size_t lengthCopy,
              const RealType* __restrict__ src, const size_t leadingSrc) {

    for(size_t idxl = 0; idxl < lengthCopy; ++idxl) {
        for(size_t idxb = 0; idxb < BlockSize; ++idxb) {
            dest[idxl * leadingDest + idxb] = src[idxb * leadingSrc + idxl];
        }
    }
}

template < class RealType >
void CheckEquality(const RealType goodMat[], const RealType testMat[], const size_t matDim) {
    for(size_t idxRow = 0; idxRow < matDim; ++idxRow) {
        for(size_t idxCol = 0; idxCol < matDim; ++idxCol) {
            if(std::abs((testMat[idxCol * matDim + idxRow] - goodMat[idxCol * matDim + idxRow]) / goodMat[idxCol * matDim + idxRow]) > (sizeof(RealType) == 4 ? 1E-6 : 1E-10)) {
                std::cout << "Error at pos " << idxRow << " " << idxCol << "\n";
                std::cout << "    Should be " << goodMat[idxCol * matDim + idxRow] << "\n";
                std::cout << "    but is " << testMat[idxCol * matDim + idxRow] << "\n";
            }
        }
    }
}


///////////////////////////////////////////////////////////////////////////////////
///
/// Scalar basic implementation (very slow)
///
///////////////////////////////////////////////////////////////////////////////////


template < class RealType >
void ScalarGemmNoBlock(const RealType* __restrict__ A, const RealType* __restrict__ B,
                       RealType* __restrict__ C, const size_t size) {

    for(size_t j = 0; j < size; ++j) {
        for(size_t i = 0; i < size; ++i) {
            for(size_t k = 0; k < size; ++k) {
                C[k * size + j] += A[k * size + i] * B[j * size + i];
            }
        }
    }
}


///////////////////////////////////////////////////////////////////////////////////
///
/// Scalar with a blocking scheme
/// (similar as http://www.netlib.org/utk/papers/autoblock/node2.html)
///
///////////////////////////////////////////////////////////////////////////////////

template < class RealType, size_t BlockSize >
void ScalarGemm(const RealType* __restrict__ A, const RealType* __restrict__ B,
                RealType* __restrict__ C, const size_t size) {

    for(size_t j_outer = 0; j_outer < size; j_outer += BlockSize) {
        for(size_t i_outer = 0; i_outer < size; i_outer += BlockSize) {
            for(size_t k = 0; k < size; ++k) {

                const RealType* __restrict__ ptrA = &A[k * size + i_outer];
                const RealType* __restrict__ ptrB = &B[j_outer * size + i_outer];
                RealType* __restrict__ ptrC       = &C[k * size + j_outer];

                for(size_t idxBlockCol = 0; idxBlockCol < BlockSize; ++idxBlockCol) {
                    RealType v = 0;
                    for(size_t idxVal = 0; idxVal < BlockSize; ++idxVal) {
                        v += ptrA[idxVal] * ptrB[idxBlockCol * size + idxVal];
                    }
                    ptrC[idxBlockCol] += v;
                }
            }
        }
    }
}


///////////////////////////////////////////////////////////////////////////////////
///
/// Scalar with a blocking scheme more close Goto-blas
///
///////////////////////////////////////////////////////////////////////////////////

template < class RealType, size_t PanelSizeK, size_t PanelSizeiA,
           size_t PanelSizejB, size_t BlockSize >
void ScalarGemmV2(const RealType* __restrict__ A, const RealType* __restrict__ B,
                  RealType* __restrict__ C, const size_t size) {

    static_assert(PanelSizeK >= BlockSize, "PanelSizeK must be greater than block");
    static_assert(PanelSizeiA >= BlockSize, "PanelSizeiA must be greater than block");
    static_assert(PanelSizejB >= BlockSize, "PanelSizejB must be greater than block");
    static_assert((PanelSizeK / BlockSize) * BlockSize == PanelSizeK, "PanelSizeK must be a multiple of block");
    static_assert((PanelSizeiA / BlockSize) * BlockSize == PanelSizeiA, "PanelSizeiA must be a multiple of block");
    static_assert((PanelSizejB / BlockSize) * BlockSize == PanelSizejB, "PanelSizejB must be a multiple of block");
    // Restrict to a multiple of panelsize for simplcity
    assert((size / PanelSizeK) * PanelSizeK == size);
    assert((size / PanelSizeiA) * PanelSizeiA == size);
    assert((size / PanelSizejB) * PanelSizejB == size);

    for(size_t ip = 0; ip < size; ip += PanelSizeiA) {
        for(size_t jp = 0; jp < size; jp += PanelSizejB) {

            for(size_t kp = 0; kp < size; kp += PanelSizeK) {

                alignas(64) RealType panelA[PanelSizeiA * PanelSizeK];
                alignas(64) RealType panelB[PanelSizeK * BlockSize];

                for(size_t jb = 0; jb < PanelSizejB; jb += BlockSize) {

                    CopyMat< RealType, BlockSize >(panelB, PanelSizeK, &B[jp * size + kp], size);

                    for(size_t ib = 0; ib < PanelSizeiA; ib += BlockSize) {

                        if(jb == 0) {
                            CopyMat< RealType, BlockSize >(&panelA[PanelSizeK * ib], PanelSizeK, &A[(ib + ip) * size + kp], size);
                        }

                        for(size_t idxRow = 0; idxRow < BlockSize; ++idxRow) {
                            for(size_t idxCol = 0; idxCol < BlockSize; ++idxCol) {
                                RealType sum = 0;
                                for(size_t idxK = 0; idxK < PanelSizeK; ++idxK) {
                                    sum += panelA[(idxRow + ib) * PanelSizeK + idxK] * panelB[idxCol * PanelSizeK + idxK];
                                }
                                C[(jp + jb + idxCol) * size + ip + ib + idxRow] += sum;
                            }
                        }
                    }
                }
            }
        }
    }
}


///////////////////////////////////////////////////////////////////////////////////
///
/// Native implementation
///
///////////////////////////////////////////////////////////////////////////////////


#ifdef INASTEMP_USE_SSE41

#include "SSE41/InaVecSSE41Double.hpp"
#include "SSE41/InaVecSSE41Float.hpp"

template < size_t PanelSizeK, size_t PanelSizeiA,
           size_t PanelSizejB, size_t VecTypeLength >
void InaVecSSE41_ScalarGemmInaV2(const float* __restrict__ A, const float* __restrict__ B,
                                 float* __restrict__ C, const size_t size) {

    const int BlockSize = VecTypeLength;

    static_assert(PanelSizeK >= BlockSize, "PanelSizeK must be greater than block");
    static_assert(PanelSizeiA >= BlockSize, "PanelSizeiA must be greater than block");
    static_assert(PanelSizejB >= BlockSize, "PanelSizejB must be greater than block");
    static_assert((PanelSizeK / BlockSize) * BlockSize == PanelSizeK, "PanelSizeK must be a multiple of block");
    static_assert((PanelSizeiA / BlockSize) * BlockSize == PanelSizeiA, "PanelSizeiA must be a multiple of block");
    static_assert((PanelSizejB / BlockSize) * BlockSize == PanelSizejB, "PanelSizejB must be a multiple of block");
    // Restrict to a multiple of panelsize for simplcity
    assert((size / PanelSizeK) * PanelSizeK == size);
    assert((size / PanelSizeiA) * PanelSizeiA == size);
    assert((size / PanelSizejB) * PanelSizejB == size);

    for(size_t ip = 0; ip < size; ip += PanelSizeiA) {
        for(size_t jp = 0; jp < size; jp += PanelSizejB) {

            for(size_t kp = 0; kp < size; kp += PanelSizeK) {

                alignas(64) float panelA[PanelSizeiA * PanelSizeK];
                alignas(64) float panelB[PanelSizeK * BlockSize];

                for(size_t jb = 0; jb < PanelSizejB; jb += BlockSize) {

                    CopyMat< float, BlockSize >(panelB, PanelSizeK, &B[jp * size + kp], size);

                    for(size_t ib = 0; ib < PanelSizeiA; ib += BlockSize) {

                        if(jb == 0) {
                            CopyMatT< float, BlockSize >(&panelA[ib], PanelSizeiA, PanelSizeK,
                                                         &A[(ib + ip) * size + kp], size);
                        }

                        __m128 sum[BlockSize];

                        for(size_t idxCol = 0; idxCol < BlockSize; ++idxCol) {
                            sum[idxCol] = _mm_set1_ps(0.);
                        }

                        for(size_t idxK = 0; idxK < PanelSizeK; ++idxK) {
                            const __m128 valA = _mm_loadu_ps(&panelA[idxK * PanelSizeiA + ib]);
                            for(size_t idxCol = 0; idxCol < BlockSize; ++idxCol) {
                                sum[idxCol] = _mm_add_ps(sum[idxCol], _mm_mul_ps(valA, _mm_loadu_ps(&panelB[idxCol * PanelSizeK + idxK])));
                            }
                        }

                        float* __restrict__ ptrC = &C[(jp + jb) * size + ip + ib];
                        for(size_t idxCol = 0; idxCol < BlockSize; ++idxCol) {
                            __m128 res = _mm_add_ps(sum[idxCol], _mm_loadu_ps(&ptrC[idxCol * size]));
                            _mm_storeu_ps(&ptrC[idxCol * size], res);
                        }
                    }
                }
            }
        }
    }
}


template < size_t PanelSizeK, size_t PanelSizeiA,
           size_t PanelSizejB, size_t VecTypeLength >
void InaVecSSE41_ScalarGemmInaV2(const double* __restrict__ A, const double* __restrict__ B,
                                 double* __restrict__ C, const size_t size) {

    const int BlockSize = VecTypeLength;

    static_assert(PanelSizeK >= BlockSize, "PanelSizeK must be greater than block");
    static_assert(PanelSizeiA >= BlockSize, "PanelSizeiA must be greater than block");
    static_assert(PanelSizejB >= BlockSize, "PanelSizejB must be greater than block");
    static_assert((PanelSizeK / BlockSize) * BlockSize == PanelSizeK, "PanelSizeK must be a multiple of block");
    static_assert((PanelSizeiA / BlockSize) * BlockSize == PanelSizeiA, "PanelSizeiA must be a multiple of block");
    static_assert((PanelSizejB / BlockSize) * BlockSize == PanelSizejB, "PanelSizejB must be a multiple of block");
    // Restrict to a multiple of panelsize for simplcity
    assert((size / PanelSizeK) * PanelSizeK == size);
    assert((size / PanelSizeiA) * PanelSizeiA == size);
    assert((size / PanelSizejB) * PanelSizejB == size);

    for(size_t ip = 0; ip < size; ip += PanelSizeiA) {
        for(size_t jp = 0; jp < size; jp += PanelSizejB) {

            for(size_t kp = 0; kp < size; kp += PanelSizeK) {

                alignas(64) double panelA[PanelSizeiA * PanelSizeK];
                alignas(64) double panelB[PanelSizeK * BlockSize];

                for(size_t jb = 0; jb < PanelSizejB; jb += BlockSize) {

                    CopyMat< double, BlockSize >(panelB, PanelSizeK, &B[jp * size + kp], size);

                    for(size_t ib = 0; ib < PanelSizeiA; ib += BlockSize) {

                        if(jb == 0) {
                            CopyMatT< double, BlockSize >(&panelA[ib], PanelSizeiA, PanelSizeK,
                                                          &A[(ib + ip) * size + kp], size);
                        }

                        __m128d sum[BlockSize];

                        for(size_t idxCol = 0; idxCol < BlockSize; ++idxCol) {
                            sum[idxCol] = _mm_set1_pd(0.);
                        }

                        for(size_t idxK = 0; idxK < PanelSizeK; ++idxK) {
                            const __m128d valA = _mm_loadu_pd(&panelA[idxK * PanelSizeiA + ib]);
                            for(size_t idxCol = 0; idxCol < BlockSize; ++idxCol) {
                                sum[idxCol] = _mm_add_pd(sum[idxCol], _mm_mul_pd(valA, _mm_loadu_pd(&panelB[idxCol * PanelSizeK + idxK])));
                            }
                        }

                        double* __restrict__ ptrC = &C[(jp + jb) * size + ip + ib];
                        for(size_t idxCol = 0; idxCol < BlockSize; ++idxCol) {
                            __m128d res = _mm_add_pd(sum[idxCol], _mm_loadu_pd(&ptrC[idxCol * size]));
                            _mm_storeu_pd(&ptrC[idxCol * size], res);
                        }
                    }
                }
            }
        }
    }
}

#endif

#ifdef INASTEMP_USE_AVX

#include "AVX/InaVecAVXDouble.hpp"
#include "AVX/InaVecAVXFloat.hpp"

template < size_t PanelSizeK, size_t PanelSizeiA,
           size_t PanelSizejB, size_t VecTypeLength >
void InaVecAVX_ScalarGemmInaV2(const float* __restrict__ A, const float* __restrict__ B,
                               float* __restrict__ C, const size_t size) {

    const int BlockSize = VecTypeLength;

    static_assert(PanelSizeK >= BlockSize, "PanelSizeK must be greater than block");
    static_assert(PanelSizeiA >= BlockSize, "PanelSizeiA must be greater than block");
    static_assert(PanelSizejB >= BlockSize, "PanelSizejB must be greater than block");
    static_assert((PanelSizeK / BlockSize) * BlockSize == PanelSizeK, "PanelSizeK must be a multiple of block");
    static_assert((PanelSizeiA / BlockSize) * BlockSize == PanelSizeiA, "PanelSizeiA must be a multiple of block");
    static_assert((PanelSizejB / BlockSize) * BlockSize == PanelSizejB, "PanelSizejB must be a multiple of block");
    // Restrict to a multiple of panelsize for simplcity
    assert((size / PanelSizeK) * PanelSizeK == size);
    assert((size / PanelSizeiA) * PanelSizeiA == size);
    assert((size / PanelSizejB) * PanelSizejB == size);

    for(size_t ip = 0; ip < size; ip += PanelSizeiA) {
        for(size_t jp = 0; jp < size; jp += PanelSizejB) {

            for(size_t kp = 0; kp < size; kp += PanelSizeK) {

                alignas(64) float panelA[PanelSizeiA * PanelSizeK];
                alignas(64) float panelB[PanelSizeK * BlockSize];

                for(size_t jb = 0; jb < PanelSizejB; jb += BlockSize) {

                    CopyMat< float, BlockSize >(panelB, PanelSizeK, &B[jp * size + kp], size);

                    for(size_t ib = 0; ib < PanelSizeiA; ib += BlockSize) {

                        if(jb == 0) {
                            CopyMatT< float, BlockSize >(&panelA[ib], PanelSizeiA, PanelSizeK,
                                                         &A[(ib + ip) * size + kp], size);
                        }

                        __m256 sum[BlockSize];

                        for(size_t idxCol = 0; idxCol < BlockSize; ++idxCol) {
                            sum[idxCol] = _mm256_set1_ps(0.);
                        }

                        for(size_t idxK = 0; idxK < PanelSizeK; ++idxK) {
                            const __m256 valA = _mm256_loadu_ps(&panelA[idxK * PanelSizeiA + ib]);
                            for(size_t idxCol = 0; idxCol < BlockSize; ++idxCol) {
                                sum[idxCol] = _mm256_add_ps(sum[idxCol], _mm256_mul_ps(valA, _mm256_loadu_ps(&panelB[idxCol * PanelSizeK + idxK])));
                            }
                        }

                        float* __restrict__ ptrC = &C[(jp + jb) * size + ip + ib];
                        for(size_t idxCol = 0; idxCol < BlockSize; ++idxCol) {
                            __m256 res = _mm256_add_ps(sum[idxCol], _mm256_loadu_ps(&ptrC[idxCol * size]));
                            _mm256_storeu_ps(&ptrC[idxCol * size], res);
                        }
                    }
                }
            }
        }
    }
}


template < size_t PanelSizeK, size_t PanelSizeiA,
           size_t PanelSizejB, size_t VecTypeLength >
void InaVecAVX_ScalarGemmInaV2(const double* __restrict__ A, const double* __restrict__ B,
                               double* __restrict__ C, const size_t size) {

    const int BlockSize = VecTypeLength;

    static_assert(PanelSizeK >= BlockSize, "PanelSizeK must be greater than block");
    static_assert(PanelSizeiA >= BlockSize, "PanelSizeiA must be greater than block");
    static_assert(PanelSizejB >= BlockSize, "PanelSizejB must be greater than block");
    static_assert((PanelSizeK / BlockSize) * BlockSize == PanelSizeK, "PanelSizeK must be a multiple of block");
    static_assert((PanelSizeiA / BlockSize) * BlockSize == PanelSizeiA, "PanelSizeiA must be a multiple of block");
    static_assert((PanelSizejB / BlockSize) * BlockSize == PanelSizejB, "PanelSizejB must be a multiple of block");
    // Restrict to a multiple of panelsize for simplcity
    assert((size / PanelSizeK) * PanelSizeK == size);
    assert((size / PanelSizeiA) * PanelSizeiA == size);
    assert((size / PanelSizejB) * PanelSizejB == size);

    for(size_t ip = 0; ip < size; ip += PanelSizeiA) {
        for(size_t jp = 0; jp < size; jp += PanelSizejB) {

            for(size_t kp = 0; kp < size; kp += PanelSizeK) {

                alignas(64) double panelA[PanelSizeiA * PanelSizeK];
                alignas(64) double panelB[PanelSizeK * BlockSize];

                for(size_t jb = 0; jb < PanelSizejB; jb += BlockSize) {

                    CopyMat< double, BlockSize >(panelB, PanelSizeK, &B[jp * size + kp], size);

                    for(size_t ib = 0; ib < PanelSizeiA; ib += BlockSize) {

                        if(jb == 0) {
                            CopyMatT< double, BlockSize >(&panelA[ib], PanelSizeiA, PanelSizeK,
                                                          &A[(ib + ip) * size + kp], size);
                        }

                        __m256d sum[BlockSize];

                        for(size_t idxCol = 0; idxCol < BlockSize; ++idxCol) {
                            sum[idxCol] = _mm256_set1_pd(0.);
                        }

                        for(size_t idxK = 0; idxK < PanelSizeK; ++idxK) {
                            const __m256d valA = _mm256_loadu_pd(&panelA[idxK * PanelSizeiA + ib]);
                            for(size_t idxCol = 0; idxCol < BlockSize; ++idxCol) {
                                sum[idxCol] = _mm256_add_pd(sum[idxCol], _mm256_mul_pd(valA, _mm256_loadu_pd(&panelB[idxCol * PanelSizeK + idxK])));
                            }
                        }

                        double* __restrict__ ptrC = &C[(jp + jb) * size + ip + ib];
                        for(size_t idxCol = 0; idxCol < BlockSize; ++idxCol) {
                            __m256d res = _mm256_add_pd(sum[idxCol], _mm256_loadu_pd(&ptrC[idxCol * size]));
                            _mm256_storeu_pd(&ptrC[idxCol * size], res);
                        }
                    }
                }
            }
        }
    }
}


#endif

#ifdef INASTEMP_USE_AVX512KNL

#include "AVX512KNL/InaVecAVX512KNLDouble.hpp"
#include "AVX512KNL/InaVecAVX512KNLFloat.hpp"

template < size_t PanelSizeK, size_t PanelSizeiA,
           size_t PanelSizejB, size_t VecTypeLength >
void InaVecAVX512KNL_ScalarGemmInaV2(const float* __restrict__ A, const float* __restrict__ B,
                                     float* __restrict__ C, const size_t size) {

    const int BlockSize = VecTypeLength;

    static_assert(PanelSizeK >= BlockSize, "PanelSizeK must be greater than block");
    static_assert(PanelSizeiA >= BlockSize, "PanelSizeiA must be greater than block");
    static_assert(PanelSizejB >= BlockSize, "PanelSizejB must be greater than block");
    static_assert((PanelSizeK / BlockSize) * BlockSize == PanelSizeK, "PanelSizeK must be a multiple of block");
    static_assert((PanelSizeiA / BlockSize) * BlockSize == PanelSizeiA, "PanelSizeiA must be a multiple of block");
    static_assert((PanelSizejB / BlockSize) * BlockSize == PanelSizejB, "PanelSizejB must be a multiple of block");
    // Restrict to a multiple of panelsize for simplcity
    assert((size / PanelSizeK) * PanelSizeK == size);
    assert((size / PanelSizeiA) * PanelSizeiA == size);
    assert((size / PanelSizejB) * PanelSizejB == size);

    for(size_t ip = 0; ip < size; ip += PanelSizeiA) {
        for(size_t jp = 0; jp < size; jp += PanelSizejB) {

            for(size_t kp = 0; kp < size; kp += PanelSizeK) {

                alignas(64) float panelA[PanelSizeiA * PanelSizeK];
                alignas(64) float panelB[PanelSizeK * BlockSize];

                for(size_t jb = 0; jb < PanelSizejB; jb += BlockSize) {

                    CopyMat< float, BlockSize >(panelB, PanelSizeK, &B[jp * size + kp], size);

                    for(size_t ib = 0; ib < PanelSizeiA; ib += BlockSize) {

                        if(jb == 0) {
                            CopyMatT< float, BlockSize >(&panelA[ib], PanelSizeiA, PanelSizeK,
                                                         &A[(ib + ip) * size + kp], size);
                        }

                        __m512 sum[BlockSize];

                        for(size_t idxCol = 0; idxCol < BlockSize; ++idxCol) {
                            sum[idxCol] = _mm512_set1_ps(0.);
                        }

                        for(size_t idxK = 0; idxK < PanelSizeK; ++idxK) {
                            const __m512 valA = _mm512_loadu_ps(&panelA[idxK * PanelSizeiA + ib]);
                            for(size_t idxCol = 0; idxCol < BlockSize; ++idxCol) {
                                sum[idxCol] = _mm512_add_ps(sum[idxCol], _mm512_mul_ps(valA, _mm512_loadu_ps(&panelB[idxCol * PanelSizeK + idxK])));
                            }
                        }

                        float* __restrict__ ptrC = &C[(jp + jb) * size + ip + ib];
                        for(size_t idxCol = 0; idxCol < BlockSize; ++idxCol) {
                            __m512 res = _mm512_add_ps(sum[idxCol], _mm512_loadu_ps(&ptrC[idxCol * size]));
                            _mm512_storeu_ps(&ptrC[idxCol * size], res);
                        }
                    }
                }
            }
        }
    }
}


template < size_t PanelSizeK, size_t PanelSizeiA,
           size_t PanelSizejB, size_t VecTypeLength >
void InaVecAVX512KNL_ScalarGemmInaV2(const double* __restrict__ A, const double* __restrict__ B,
                                     double* __restrict__ C, const size_t size) {

    const int BlockSize = VecTypeLength;

    static_assert(PanelSizeK >= BlockSize, "PanelSizeK must be greater than block");
    static_assert(PanelSizeiA >= BlockSize, "PanelSizeiA must be greater than block");
    static_assert(PanelSizejB >= BlockSize, "PanelSizejB must be greater than block");
    static_assert((PanelSizeK / BlockSize) * BlockSize == PanelSizeK, "PanelSizeK must be a multiple of block");
    static_assert((PanelSizeiA / BlockSize) * BlockSize == PanelSizeiA, "PanelSizeiA must be a multiple of block");
    static_assert((PanelSizejB / BlockSize) * BlockSize == PanelSizejB, "PanelSizejB must be a multiple of block");
    // Restrict to a multiple of panelsize for simplcity
    assert((size / PanelSizeK) * PanelSizeK == size);
    assert((size / PanelSizeiA) * PanelSizeiA == size);
    assert((size / PanelSizejB) * PanelSizejB == size);

    for(size_t ip = 0; ip < size; ip += PanelSizeiA) {
        for(size_t jp = 0; jp < size; jp += PanelSizejB) {

            for(size_t kp = 0; kp < size; kp += PanelSizeK) {

                alignas(64) double panelA[PanelSizeiA * PanelSizeK];
                alignas(64) double panelB[PanelSizeK * BlockSize];

                for(size_t jb = 0; jb < PanelSizejB; jb += BlockSize) {

                    CopyMat< double, BlockSize >(panelB, PanelSizeK, &B[jp * size + kp], size);

                    for(size_t ib = 0; ib < PanelSizeiA; ib += BlockSize) {

                        if(jb == 0) {
                            CopyMatT< double, BlockSize >(&panelA[ib], PanelSizeiA, PanelSizeK,
                                                          &A[(ib + ip) * size + kp], size);
                        }

                        __m512d sum[BlockSize];

                        for(size_t idxCol = 0; idxCol < BlockSize; ++idxCol) {
                            sum[idxCol] = _mm512_set1_pd(0.);
                        }

                        for(size_t idxK = 0; idxK < PanelSizeK; ++idxK) {
                            const __m512d valA = _mm512_loadu_pd(&panelA[idxK * PanelSizeiA + ib]);
                            for(size_t idxCol = 0; idxCol < BlockSize; ++idxCol) {
                                sum[idxCol] = _mm512_add_pd(sum[idxCol], _mm512_mul_pd(valA, _mm512_loadu_pd(&panelB[idxCol * PanelSizeK + idxK])));
                            }
                        }

                        double* __restrict__ ptrC = &C[(jp + jb) * size + ip + ib];
                        for(size_t idxCol = 0; idxCol < BlockSize; ++idxCol) {
                            __m512d res = _mm512_add_pd(sum[idxCol], _mm512_loadu_pd(&ptrC[idxCol * size]));
                            _mm512_storeu_pd(&ptrC[idxCol * size], res);
                        }
                    }
                }
            }
        }
    }
}
#endif

#ifdef INASTEMP_USE_ALTIVEC

#include "ALTIVEC/InaVecALTIVECDouble.hpp"
#include "ALTIVEC/InaVecALTIVECFloat.hpp"

template < size_t PanelSizeK, size_t PanelSizeiA,
           size_t PanelSizejB, size_t VecTypeLength >
void InaVecALTIVEC_ScalarGemmInaV2(const float* __restrict__ A, const float* __restrict__ B,
                                   float* __restrict__ C, const size_t size) {

    const int BlockSize = VecTypeLength;

    static_assert(PanelSizeK >= BlockSize, "PanelSizeK must be greater than block");
    static_assert(PanelSizeiA >= BlockSize, "PanelSizeiA must be greater than block");
    static_assert(PanelSizejB >= BlockSize, "PanelSizejB must be greater than block");
    static_assert((PanelSizeK / BlockSize) * BlockSize == PanelSizeK, "PanelSizeK must be a multiple of block");
    static_assert((PanelSizeiA / BlockSize) * BlockSize == PanelSizeiA, "PanelSizeiA must be a multiple of block");
    static_assert((PanelSizejB / BlockSize) * BlockSize == PanelSizejB, "PanelSizejB must be a multiple of block");
    // Restrict to a multiple of panelsize for simplcity
    assert((size / PanelSizeK) * PanelSizeK == size);
    assert((size / PanelSizeiA) * PanelSizeiA == size);
    assert((size / PanelSizejB) * PanelSizejB == size);

    for(size_t ip = 0; ip < size; ip += PanelSizeiA) {
        for(size_t jp = 0; jp < size; jp += PanelSizejB) {

            for(size_t kp = 0; kp < size; kp += PanelSizeK) {

                alignas(64) float panelA[PanelSizeiA * PanelSizeK];
                alignas(64) float panelB[PanelSizeK * BlockSize];

                for(size_t jb = 0; jb < PanelSizejB; jb += BlockSize) {

                    CopyMat< float, BlockSize >(panelB, PanelSizeK, &B[jp * size + kp], size);

                    for(size_t ib = 0; ib < PanelSizeiA; ib += BlockSize) {

                        if(jb == 0) {
                            CopyMatT< float, BlockSize >(&panelA[ib], PanelSizeiA, PanelSizeK,
                                                         &A[(ib + ip) * size + kp], size);
                        }

                        __vector float sum[BlockSize];

                        for(size_t idxCol = 0; idxCol < BlockSize; ++idxCol) {
                            sum[idxCol] = vec_splats(0.f);
                        }

                        for(size_t idxK = 0; idxK < PanelSizeK; ++idxK) {
                            const __vector float valA = vec_xl(0, &panelA[idxK * PanelSizeiA + ib]);
                            for(size_t idxCol = 0; idxCol < BlockSize; ++idxCol) {
                                sum[idxCol] += valA * vec_xl(0, &panelB[idxCol * PanelSizeK + idxK]);
                            }
                        }

                        float* __restrict__ ptrC = &C[(jp + jb) * size + ip + ib];
                        for(size_t idxCol = 0; idxCol < BlockSize; ++idxCol) {
                            __vector float res = sum[idxCol] + vec_xl(0, &ptrC[idxCol * size]);
                            vec_xst(res, 0, &ptrC[idxCol * size]);
                        }
                    }
                }
            }
        }
    }
}


template < size_t PanelSizeK, size_t PanelSizeiA,
           size_t PanelSizejB, size_t VecTypeLength >
void InaVecALTIVEC_ScalarGemmInaV2(const double* __restrict__ A, const double* __restrict__ B,
                                   double* __restrict__ C, const size_t size) {

    const int BlockSize = VecTypeLength;

    static_assert(PanelSizeK >= BlockSize, "PanelSizeK must be greater than block");
    static_assert(PanelSizeiA >= BlockSize, "PanelSizeiA must be greater than block");
    static_assert(PanelSizejB >= BlockSize, "PanelSizejB must be greater than block");
    static_assert((PanelSizeK / BlockSize) * BlockSize == PanelSizeK, "PanelSizeK must be a multiple of block");
    static_assert((PanelSizeiA / BlockSize) * BlockSize == PanelSizeiA, "PanelSizeiA must be a multiple of block");
    static_assert((PanelSizejB / BlockSize) * BlockSize == PanelSizejB, "PanelSizejB must be a multiple of block");
    // Restrict to a multiple of panelsize for simplcity
    assert((size / PanelSizeK) * PanelSizeK == size);
    assert((size / PanelSizeiA) * PanelSizeiA == size);
    assert((size / PanelSizejB) * PanelSizejB == size);

    for(size_t ip = 0; ip < size; ip += PanelSizeiA) {
        for(size_t jp = 0; jp < size; jp += PanelSizejB) {

            for(size_t kp = 0; kp < size; kp += PanelSizeK) {

                alignas(64) double panelA[PanelSizeiA * PanelSizeK];
                alignas(64) double panelB[PanelSizeK * BlockSize];

                for(size_t jb = 0; jb < PanelSizejB; jb += BlockSize) {

                    CopyMat< double, BlockSize >(panelB, PanelSizeK, &B[jp * size + kp], size);

                    for(size_t ib = 0; ib < PanelSizeiA; ib += BlockSize) {

                        if(jb == 0) {
                            CopyMatT< double, BlockSize >(&panelA[ib], PanelSizeiA, PanelSizeK,
                                                          &A[(ib + ip) * size + kp], size);
                        }

                        __vector double sum[BlockSize];

                        for(size_t idxCol = 0; idxCol < BlockSize; ++idxCol) {
                            sum[idxCol] = vec_splats(0.);
                        }

                        for(size_t idxK = 0; idxK < PanelSizeK; ++idxK) {
                            const __vector double valA = vec_xl(0, &panelA[idxK * PanelSizeiA + ib]);
                            for(size_t idxCol = 0; idxCol < BlockSize; ++idxCol) {
                                sum[idxCol] += valA * vec_xl(0, &panelB[idxCol * PanelSizeK + idxK]);
                            }
                        }

                        double* __restrict__ ptrC = &C[(jp + jb) * size + ip + ib];
                        for(size_t idxCol = 0; idxCol < BlockSize; ++idxCol) {
                            __vector double res = sum[idxCol] + vec_xl(0, &ptrC[idxCol * size]);
                            vec_xst(reinterpret_cast< __vector unsigned int >(res), 0, reinterpret_cast< unsigned int* >(&ptrC[idxCol * size]));
                        }
                    }
                }
            }
        }
    }
}


template < size_t PanelSizeK, size_t PanelSizeiA,
           size_t PanelSizejB >
void InaVecALTIVEC_ScalarGemmIna(const double* __restrict__ A, const double* __restrict__ B,
                                 double* __restrict__ C, const size_t size) {

    const int BlockSize = 2;

    static_assert(PanelSizeK >= BlockSize, "PanelSizeK must be greater than block");
    static_assert(PanelSizeiA >= BlockSize, "PanelSizeiA must be greater than block");
    static_assert(PanelSizejB >= BlockSize, "PanelSizejB must be greater than block");
    static_assert((PanelSizeK / BlockSize) * BlockSize == PanelSizeK, "PanelSizeK must be a multiple of block");
    static_assert((PanelSizeiA / BlockSize) * BlockSize == PanelSizeiA, "PanelSizeiA must be a multiple of block");
    static_assert((PanelSizejB / BlockSize) * BlockSize == PanelSizejB, "PanelSizejB must be a multiple of block");
    // Restrict to a multiple of panelsize for simplcity
    assert((size / PanelSizeK) * PanelSizeK == size);
    assert((size / PanelSizeiA) * PanelSizeiA == size);
    assert((size / PanelSizejB) * PanelSizejB == size);

    for(size_t ip = 0; ip < size; ip += PanelSizeiA) {
        for(size_t jp = 0; jp < size; jp += PanelSizejB) {

            for(size_t kp = 0; kp < size; kp += PanelSizeK) {

                alignas(64) double panelA[PanelSizeiA * PanelSizeK];
                alignas(64) double panelB[PanelSizeK * BlockSize];

                for(size_t jb = 0; jb < PanelSizejB; jb += BlockSize) {

                    CopyMat< double, BlockSize >(panelB, PanelSizeK, &B[jp * size + kp], size);

                    for(size_t ib = 0; ib < PanelSizeiA; ib += BlockSize) {

                        if(jb == 0) {
                            CopyMat< double, BlockSize >(&panelA[PanelSizeK * ib], PanelSizeK, &A[(ib + ip) * size + kp], size);
                        }

                        for(size_t idxRow = 0; idxRow < BlockSize; ++idxRow) {
                            for(size_t idxCol = 0; idxCol < BlockSize; ++idxCol) {
                                __vector double sum = vec_splats(0.);
                                for(size_t idxK = 0; idxK < PanelSizeK; idxK += BlockSize) {
                                    sum += vec_xl(0, &panelA[(idxRow + ib) * PanelSizeK + idxK]) * vec_xl(0, &panelB[idxCol * PanelSizeK + idxK]);
                                }
                                double ptr[2];
                                vec_xst(reinterpret_cast< __vector unsigned int >(sum), 0, reinterpret_cast< unsigned int* >(ptr));
                                C[(jp + jb + idxCol) * size + ip + ib + idxRow] += ptr[0] + ptr[1];
                            }
                        }
                    }
                }
            }
        }
    }
}

template < size_t PanelSizeK, size_t PanelSizeiA,
           size_t PanelSizejB >
void InaVecALTIVEC_ScalarGemmIna(const float* __restrict__ A, const float* __restrict__ B,
                                 float* __restrict__ C, const size_t size) {

    const int BlockSize = 4;

    static_assert(PanelSizeK >= BlockSize, "PanelSizeK must be greater than block");
    static_assert(PanelSizeiA >= BlockSize, "PanelSizeiA must be greater than block");
    static_assert(PanelSizejB >= BlockSize, "PanelSizejB must be greater than block");
    static_assert((PanelSizeK / BlockSize) * BlockSize == PanelSizeK, "PanelSizeK must be a multiple of block");
    static_assert((PanelSizeiA / BlockSize) * BlockSize == PanelSizeiA, "PanelSizeiA must be a multiple of block");
    static_assert((PanelSizejB / BlockSize) * BlockSize == PanelSizejB, "PanelSizejB must be a multiple of block");
    // Restrict to a multiple of panelsize for simplcity
    assert((size / PanelSizeK) * PanelSizeK == size);
    assert((size / PanelSizeiA) * PanelSizeiA == size);
    assert((size / PanelSizejB) * PanelSizejB == size);

    for(size_t ip = 0; ip < size; ip += PanelSizeiA) {
        for(size_t jp = 0; jp < size; jp += PanelSizejB) {

            for(size_t kp = 0; kp < size; kp += PanelSizeK) {

                alignas(64) float panelA[PanelSizeiA * PanelSizeK];
                alignas(64) float panelB[PanelSizeK * BlockSize];

                for(size_t jb = 0; jb < PanelSizejB; jb += BlockSize) {

                    CopyMat< float, BlockSize >(panelB, PanelSizeK, &B[jp * size + kp], size);

                    for(size_t ib = 0; ib < PanelSizeiA; ib += BlockSize) {

                        if(jb == 0) {
                            CopyMat< float, BlockSize >(&panelA[PanelSizeK * ib], PanelSizeK, &A[(ib + ip) * size + kp], size);
                        }

                        for(size_t idxRow = 0; idxRow < BlockSize; ++idxRow) {
                            for(size_t idxCol = 0; idxCol < BlockSize; ++idxCol) {
                                __vector float sum = vec_splats(0.f);
                                for(size_t idxK = 0; idxK < PanelSizeK; idxK += BlockSize) {
                                    sum += vec_xl(0, &panelA[(idxRow + ib) * PanelSizeK + idxK]) * vec_xl(0, &panelB[idxCol * PanelSizeK + idxK]);
                                }
                                float ptr[4];
                                vec_xst(sum, 0, ptr);
                                C[(jp + jb + idxCol) * size + ip + ib + idxRow] += ptr[0] + ptr[1] + ptr[2] + ptr[3];
                            }
                        }
                    }
                }
            }
        }
    }
}

#endif


///////////////////////////////////////////////////////////////////////////////////
///
/// InaVec with a blocking scheme more close Goto-blas
///
///////////////////////////////////////////////////////////////////////////////////

template < class RealType, size_t PanelSizeK, size_t PanelSizeiA,
           size_t PanelSizejB, class VecType >
void ScalarGemmIna(const RealType* __restrict__ A, const RealType* __restrict__ B,
                   RealType* __restrict__ C, const size_t size) {

    const int BlockSize = VecType::VecLength;

    static_assert(PanelSizeK >= BlockSize, "PanelSizeK must be greater than block");
    static_assert(PanelSizeiA >= BlockSize, "PanelSizeiA must be greater than block");
    static_assert(PanelSizejB >= BlockSize, "PanelSizejB must be greater than block");
    static_assert((PanelSizeK / BlockSize) * BlockSize == PanelSizeK, "PanelSizeK must be a multiple of block");
    static_assert((PanelSizeiA / BlockSize) * BlockSize == PanelSizeiA, "PanelSizeiA must be a multiple of block");
    static_assert((PanelSizejB / BlockSize) * BlockSize == PanelSizejB, "PanelSizejB must be a multiple of block");
    // Restrict to a multiple of panelsize for simplcity
    assert((size / PanelSizeK) * PanelSizeK == size);
    assert((size / PanelSizeiA) * PanelSizeiA == size);
    assert((size / PanelSizejB) * PanelSizejB == size);

    for(size_t ip = 0; ip < size; ip += PanelSizeiA) {
        for(size_t jp = 0; jp < size; jp += PanelSizejB) {

            for(size_t kp = 0; kp < size; kp += PanelSizeK) {

                alignas(64) RealType panelA[PanelSizeiA * PanelSizeK];
                alignas(64) RealType panelB[PanelSizeK * BlockSize];

                for(size_t jb = 0; jb < PanelSizejB; jb += BlockSize) {

                    CopyMat< RealType, BlockSize >(panelB, PanelSizeK, &B[jp * size + kp], size);

                    for(size_t ib = 0; ib < PanelSizeiA; ib += BlockSize) {

                        if(jb == 0) {
                            CopyMat< RealType, BlockSize >(&panelA[PanelSizeK * ib], PanelSizeK, &A[(ib + ip) * size + kp], size);
                        }

                        for(size_t idxRow = 0; idxRow < BlockSize; ++idxRow) {
                            for(size_t idxCol = 0; idxCol < BlockSize; ++idxCol) {
                                VecType sum = 0.;
                                for(size_t idxK = 0; idxK < PanelSizeK; idxK += BlockSize) {
                                    sum += VecType(&panelA[(idxRow + ib) * PanelSizeK + idxK]) * VecType(&panelB[idxCol * PanelSizeK + idxK]);
                                }
                                C[(jp + jb + idxCol) * size + ip + ib + idxRow] += sum.horizontalSum();
                            }
                        }
                    }
                }
            }
        }
    }
}


///////////////////////////////////////////////////////////////////////////////////
///
/// InaVec with a blocking scheme more close Goto-blas
/// but with more instructions
///
///////////////////////////////////////////////////////////////////////////////////

template < class RealType, size_t PanelSizeK, size_t PanelSizeiA,
           size_t PanelSizejB, class VecType >
void ScalarGemmInaV2(const RealType* __restrict__ A, const RealType* __restrict__ B,
                     RealType* __restrict__ C, const size_t size) {

    const int BlockSize = VecType::VecLength;

    static_assert(PanelSizeK >= BlockSize, "PanelSizeK must be greater than block");
    static_assert(PanelSizeiA >= BlockSize, "PanelSizeiA must be greater than block");
    static_assert(PanelSizejB >= BlockSize, "PanelSizejB must be greater than block");
    static_assert((PanelSizeK / BlockSize) * BlockSize == PanelSizeK, "PanelSizeK must be a multiple of block");
    static_assert((PanelSizeiA / BlockSize) * BlockSize == PanelSizeiA, "PanelSizeiA must be a multiple of block");
    static_assert((PanelSizejB / BlockSize) * BlockSize == PanelSizejB, "PanelSizejB must be a multiple of block");
    // Restrict to a multiple of panelsize for simplcity
    assert((size / PanelSizeK) * PanelSizeK == size);
    assert((size / PanelSizeiA) * PanelSizeiA == size);
    assert((size / PanelSizejB) * PanelSizejB == size);

    for(size_t ip = 0; ip < size; ip += PanelSizeiA) {
        for(size_t jp = 0; jp < size; jp += PanelSizejB) {

            for(size_t kp = 0; kp < size; kp += PanelSizeK) {

                alignas(64) RealType panelA[PanelSizeiA * PanelSizeK];
                alignas(64) RealType panelB[PanelSizeK * BlockSize];

                for(size_t jb = 0; jb < PanelSizejB; jb += BlockSize) {

                    CopyMat< RealType, BlockSize >(panelB, PanelSizeK, &B[jp * size + kp], size);

                    for(size_t ib = 0; ib < PanelSizeiA; ib += BlockSize) {

                        if(jb == 0) {
                            CopyMatT< RealType, BlockSize >(&panelA[ib], PanelSizeiA, PanelSizeK,
                                                            &A[(ib + ip) * size + kp], size);
                        }

                        VecType sum[BlockSize];

                        for(size_t idxCol = 0; idxCol < BlockSize; ++idxCol) {
                            sum[idxCol] = 0.;
                        }

                        for(size_t idxK = 0; idxK < PanelSizeK; ++idxK) {
                            const VecType valA(&panelA[idxK * PanelSizeiA + ib]);
                            for(size_t idxCol = 0; idxCol < BlockSize; ++idxCol) {
                                sum[idxCol] += valA * VecType(panelB[idxCol * PanelSizeK + idxK]);
                            }
                        }

                        RealType* __restrict__ ptrC = &C[(jp + jb) * size + ip + ib];
                        for(size_t idxCol = 0; idxCol < BlockSize; ++idxCol) {
                            VecType res = sum[idxCol] + VecType(&ptrC[idxCol * size]);
                            res.storeInArray(&ptrC[idxCol * size]);
                        }
                    }
                }
            }
        }
    }
}


template < class VecType, class RealType, const size_t PanelSizeA, const size_t PanelSizeB, const size_t PanelSizeK >
void ComputeGemmIna(const size_t NbOverLoop, const size_t matDim, const size_t nbFlops,
                    const RealType A[], const RealType B[]) {
    {
        std::unique_ptr< RealType[] > CIna(new RealType[matDim * matDim]);
        memset(CIna.get(), 0, sizeof(RealType) * matDim * matDim);

        InaTimer timer;

        for(size_t idxLoop = 0; idxLoop < NbOverLoop; ++idxLoop) {
            ScalarGemmIna< RealType, PanelSizeK, PanelSizeA, PanelSizeB,
                           VecType >(A, B, CIna.get(), matDim);
        }

        timer.stop();
        std::cout << "Vector " << VecType::GetName() << " for size " << matDim
                  << " took " << timer.getElapsed() << "s (" << (double(NbOverLoop * nbFlops) / timer.getElapsed()) / 1E9 << "GFlop/s)\n";
    }

    /////////////////////////////////////////////////////////////

    {
        std::unique_ptr< RealType[] > CIna(new RealType[matDim * matDim]);
        memset(CIna.get(), 0, sizeof(RealType) * matDim * matDim);

        InaTimer timer;

        for(size_t idxLoop = 0; idxLoop < NbOverLoop; ++idxLoop) {
            ScalarGemmInaV2< RealType, PanelSizeK, PanelSizeA, PanelSizeB,
                             VecType >(A, B, CIna.get(), matDim);
        }

        timer.stop();
        std::cout << "Vector V2 " << VecType::GetName() << " for size " << matDim
                  << " took " << timer.getElapsed() << "s (" << (double(NbOverLoop * nbFlops) / timer.getElapsed()) / 1E9 << "GFlop/s)\n";
    }
}


template < class RealType >
void compareGemmTime(const size_t NbOverLoop, const size_t matDim) {
    const size_t nbFlops = matDim * matDim * matDim * 2;

    std::unique_ptr< RealType[] > A(new RealType[matDim * matDim]);
    std::unique_ptr< RealType[] > B(new RealType[matDim * matDim]);

    for(size_t idxVal = 0; idxVal < matDim * matDim; ++idxVal) {
        A[idxVal] = B[idxVal] = RealType((idxVal + 1) % matDim);
    }

    /////////////////////////////////////////////////////////////

    std::unique_ptr< RealType[] > CScalarNoBlock(new RealType[matDim * matDim]);
    memset(CScalarNoBlock.get(), 0, sizeof(RealType) * matDim * matDim);
    if(matDim < 1000) {

        InaTimer timer;

        for(size_t idxLoop = 0; idxLoop < NbOverLoop; ++idxLoop) {
            ScalarGemmNoBlock< RealType >(A.get(), B.get(), CScalarNoBlock.get(), matDim);
        }

        timer.stop();
        std::cout << "Scalar no block for size " << matDim
                  << " took " << timer.getElapsed() << "s (" << (double(NbOverLoop * nbFlops) / timer.getElapsed()) / 1E9 << "GFlop/s)\n";
    } else {
        std::cout << "Scalar no block for size skipped - matrix too large which makes this algorithm too slow!\n";
    }

    /////////////////////////////////////////////////////////////

    {
        std::unique_ptr< RealType[] > CScalar(new RealType[matDim * matDim]);
        memset(CScalar.get(), 0, sizeof(RealType) * matDim * matDim);

        InaTimer timer;

        for(size_t idxLoop = 0; idxLoop < NbOverLoop; ++idxLoop) {
            ScalarGemm< RealType, InaVecBestType< RealType >::VecLength >(A.get(), B.get(), CScalar.get(), matDim);
        }

        timer.stop();
        std::cout << "Scalar for size " << matDim
                  << " took " << timer.getElapsed() << "s (" << (double(NbOverLoop * nbFlops) / timer.getElapsed()) / 1E9 << "GFlop/s)\n";
    }

    /////////////////////////////////////////////////////////////

    const size_t PanelSizeA = 32;
    const size_t PanelSizeB = 32;
    const size_t PanelSizeK = 64;
    {
        std::unique_ptr< RealType[] > CScalar(new RealType[matDim * matDim]);
        memset(CScalar.get(), 0, sizeof(RealType) * matDim * matDim);

        InaTimer timer;

        for(size_t idxLoop = 0; idxLoop < NbOverLoop; ++idxLoop) {
            ScalarGemmV2< RealType, PanelSizeK, PanelSizeA, PanelSizeB,
                          InaVecBestType< RealType >::VecLength >(A.get(), B.get(), CScalar.get(), matDim);
        }

        timer.stop();
        std::cout << "Scalar for size " << matDim
                  << " took " << timer.getElapsed() << "s (" << (double(NbOverLoop * nbFlops) / timer.getElapsed()) / 1E9 << "GFlop/s)\n";
    }

    /////////////////////////////////////////////////////////////

    ComputeGemmIna< InaVecSCALAR< RealType >, RealType, PanelSizeA, PanelSizeB, PanelSizeK >(NbOverLoop, matDim, nbFlops, A.get(), B.get());

/////////////////////////////////////////////////////////////
#ifdef INASTEMP_USE_SSE41
    ComputeGemmIna< InaVecSSE41< RealType >, RealType, PanelSizeA, PanelSizeB, PanelSizeK >(NbOverLoop, matDim, nbFlops, A.get(), B.get());

    {
        std::unique_ptr< RealType[] > CIna(new RealType[matDim * matDim]);
        memset(CIna.get(), 0, sizeof(RealType) * matDim * matDim);

        InaTimer timer;

        for(size_t idxLoop = 0; idxLoop < NbOverLoop; ++idxLoop) {
            InaVecSSE41_ScalarGemmInaV2< PanelSizeK, PanelSizeA, PanelSizeB,
                                         InaVecSSE41< RealType >::VecLength >(A.get(), B.get(), CIna.get(), matDim);
        }

        timer.stop();
        std::cout << "Vector V2 "
                  << "SSE41"
                  << " for size " << matDim
                  << " took " << timer.getElapsed() << "s (" << (double(NbOverLoop * nbFlops) / timer.getElapsed()) / 1E9 << "GFlop/s)\n";
    }
#endif


/////////////////////////////////////////////////////////////
#ifdef INASTEMP_USE_AVX
    ComputeGemmIna< InaVecAVX< RealType >, RealType, PanelSizeA, PanelSizeB, PanelSizeK >(NbOverLoop, matDim, nbFlops, A.get(), B.get());

    {
        std::unique_ptr< RealType[] > CIna(new RealType[matDim * matDim]);
        memset(CIna.get(), 0, sizeof(RealType) * matDim * matDim);

        InaTimer timer;

        for(size_t idxLoop = 0; idxLoop < NbOverLoop; ++idxLoop) {
            InaVecAVX_ScalarGemmInaV2< PanelSizeK, PanelSizeA, PanelSizeB,
                                       InaVecAVX< RealType >::VecLength >(A.get(), B.get(), CIna.get(), matDim);
        }


        timer.stop();
        std::cout << "Vector V2 "
                  << "AVX"
                  << " for size " << matDim
                  << " took " << timer.getElapsed() << "s (" << (double(NbOverLoop * nbFlops) / timer.getElapsed()) / 1E9 << "GFlop/s)\n";
    }
#endif

/////////////////////////////////////////////////////////////
#ifdef INASTEMP_USE_AVX512KNL
    ComputeGemmIna< InaVecAVX512KNL< RealType >, RealType, PanelSizeA, PanelSizeB, PanelSizeK >(NbOverLoop, matDim, nbFlops, A.get(), B.get());

    {
        std::unique_ptr< RealType[] > CIna(new RealType[matDim * matDim]);
        memset(CIna.get(), 0, sizeof(RealType) * matDim * matDim);

        InaTimer timer;

        for(size_t idxLoop = 0; idxLoop < NbOverLoop; ++idxLoop) {
            InaVecAVX512KNL_ScalarGemmInaV2< PanelSizeK, PanelSizeA, PanelSizeB,
                                             InaVecAVX512KNL< RealType >::VecLength >(A.get(), B.get(), CIna.get(), matDim);
        }

        timer.stop();
        std::cout << "Vector V2 "
                  << "AVX512KNL"
                  << " for size " << matDim
                  << " took " << timer.getElapsed() << "s (" << (double(NbOverLoop * nbFlops) / timer.getElapsed()) / 1E9 << "GFlop/s)\n";
    }
#endif

/////////////////////////////////////////////////////////////
#ifdef INASTEMP_USE_ALTIVEC
    ComputeGemmIna< InaVecALTIVEC< RealType >, RealType, PanelSizeA, PanelSizeB, PanelSizeK >(NbOverLoop, matDim, nbFlops, A.get(), B.get());

    {
        std::unique_ptr< RealType[] > CIna(new RealType[matDim * matDim]);
        memset(CIna.get(), 0, sizeof(RealType) * matDim * matDim);

        InaTimer timer;

        for(size_t idxLoop = 0; idxLoop < NbOverLoop; ++idxLoop) {
            InaVecALTIVEC_ScalarGemmInaV2< PanelSizeK, PanelSizeA, PanelSizeB,
                                           InaVecALTIVEC< RealType >::VecLength >(A.get(), B.get(), CIna.get(), matDim);
        }

        timer.stop();
        std::cout << "Vector V2 "
                  << "ALTIVEC"
                  << " for size " << matDim
                  << " took " << timer.getElapsed() << "s (" << (double(NbOverLoop * nbFlops) / timer.getElapsed()) / 1E9 << "GFlop/s)\n";
    }
    {
        std::unique_ptr< RealType[] > CIna(new RealType[matDim * matDim]);
        memset(CIna.get(), 0, sizeof(RealType) * matDim * matDim);

        InaTimer timer;

        for(size_t idxLoop = 0; idxLoop < NbOverLoop; ++idxLoop) {
            InaVecALTIVEC_ScalarGemmIna< PanelSizeK, PanelSizeA, PanelSizeB >(A.get(), B.get(), CIna.get(), matDim);
        }

        timer.stop();
        std::cout << "Vector "
                  << "ALTIVEC"
                  << " for size " << matDim
                  << " took " << timer.getElapsed() << "s (" << (double(NbOverLoop * nbFlops) / timer.getElapsed()) / 1E9 << "GFlop/s)\n";
    }
#endif
}

int main(int /*argc*/, char* /*argv*/ []) {
    std::cout << "[INFO] This program runs the computation of gemm using scalar, intrinsic vectors or inastemp vectors. \n";

    const size_t NbOverLoop   = 3;
    const size_t matDimDouble = 2048;
    const size_t matDimFloat  = 2048 + 512;
    std::cout << "[INFO] It will compute square matrix-matrix product of dim " << matDimFloat << " in float, and " << matDimDouble << " in double. \n";
    std::cout << "[INFO] This process will be done " << NbOverLoop << " times. \n";

    std::cout << "In Float:" << std::endl;
    compareGemmTime< float >(NbOverLoop, matDimFloat);

    std::cout << "In Double:" << std::endl;
    compareGemmTime< double >(NbOverLoop, matDimDouble);

    return 0;
}
