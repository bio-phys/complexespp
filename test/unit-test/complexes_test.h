#ifndef COMPLEXES_TEST_H
#define COMPLEXES_TEST_H

#include "gtest/gtest.h"
#include <unistd.h>

#include "energy/forcefield.h"
#include "energy/pairparameter.h"

const int RAND_SEED = 1194442884;

inline std::string get_current_dir() {
  auto buf = std::array<char, 8192>();
  const auto err = std::make_unique<char>(*getcwd(buf.data(), buf.size()));
  if (*err == ERANGE) {
    throw std::ios_base::failure(
        "Test could'nt get name of current directory.\n");
  }
  return std::string(buf.data());
}

template <typename T>
void EXPECT_EQ_VEC(const util::vec<T>& a, const util::vec<T>& b) {
  for (auto i = 0; i < a.size(); ++i) {
    EXPECT_EQ(a[i], b[i]);
  }
}

template <typename T>
void EXPECT_NEAR_VEC(const util::vec<T>& a, const util::vec<T>& b,
                     const double prec) {
  for (auto i = 0; i < a.size(); ++i) {
    EXPECT_NEAR(a[i], b[i], prec);
  }
}

#define EXPECT_THROW_MESSAGE(statement, exception, message)  \
  try {                                                      \
    statement;                                               \
    EXPECT_FALSE(true); /*hey you didn't throw*/             \
  } catch (exception const& e) {                             \
    EXPECT_STREQ(message, e.what());                         \
  } catch (...) {                                            \
    EXPECT_FALSE(true); /* whopsie caught wrong exception */ \
  }

#define EXPECT_THROW_MESSAGE_CONTAINS(statement, exception, message)       \
  try {                                                                    \
    statement;                                                             \
    EXPECT_FALSE(true); /*hey you didn't throw*/                           \
  } catch (exception const& e) {                                           \
    EXPECT_TRUE(std::string(e.what()).find(message) != std::string::npos); \
  } catch (...) {                                                          \
    EXPECT_FALSE(true); /* whopsie caught wrong exception */               \
  }

const std::string cwd = get_current_dir();
const std::string dataDir = cwd + "/data/";

inline energy::PairParameter<double> constantPairParam(const double fac,
                                                       const int n) {
  auto arr = util::rArray(n, n);
  for (auto i = 0; i < arr.rows(); ++i) {
    for (auto j = 0; j < arr.rows(); ++j) {
      arr(i, j) = fac;
    }
  }
  return energy::PairParameter<double>(arr);
}

inline energy::ForceField dummy_forcefield(const int n,
                                           const double fac_interaction,
                                           const double fac_radius,
                                           const double debyeLength,
                                           const double dielectricContsant) {
  const auto interaction = constantPairParam(fac_interaction, n);
  const auto radius = constantPairParam(fac_radius, n);

  const auto membrane =
      std::vector<std::array<double, 8>>{{{1, 1, 1, 1, 1, 1, 1, 1}}};

  return energy::ForceField({"ALA"}, interaction, radius,
                            std::vector<double>(n, 1), membrane, debyeLength,
                            dielectricContsant, 1);
}

template <typename T, class Func>
std::vector<T> cdf(const std::vector<T>& vals, Func f) {
  auto res = std::vector<T>();
  res.reserve(vals.size());
  for (const auto x : vals) {
    res.push_back(f(x));
  }
  return res;
}

// Copy of scipy implementation
template <typename T>
T kolmogorov(T y) {
  if (y < 1.1e-16) {
    return 1;
  }
  const auto x = -2 * y * y;

  auto sign = T(1);
  auto p = T(0);
  auto r = T(1);
  auto t = T(0);

  do {
    t = std::exp(x * r * r);
    p += sign * t;
    if (t == 0) {
      break;
    }
    r += 1;
    sign = -sign;
  } while ((t / p) > 1.1e-16);

  return p + p;
}

// Need rvs.size() > 2666
// Copy of scipy implementation
// returns p-value if rvs was sampled from CDF
template <typename T, class Func>
double kstest(std::vector<T> rvs, Func f) {
  std::sort(rvs.begin(), rvs.end());
  const auto N = rvs.size();
  const auto cdfvals = cdf<T>(rvs, f);

  double Dplus = 0;
  for (auto i = 0u; i < N; ++i) {
    const auto est_cdf = (i + 1) / static_cast<double>(N);
    const auto diff = est_cdf - cdfvals[i];
    if (diff > Dplus) {
      Dplus = diff;
    }
  }

  double Dmin = 0;
  for (auto i = 0u; i < N; ++i) {
    const auto est_cdf = i / static_cast<double>(N);
    const auto diff = cdfvals[i] - est_cdf;
    if (diff > Dmin) {
      Dmin = diff;
    }
  }

  const auto D = std::max(Dplus, Dmin);
  return kolmogorov(D * std::sqrt(N));
}

#define MY_EXPECT_DOUBLE_EQ(X, Y, EPS)                                        \
  if ((Y) == 0)                                                               \
    EXPECT_EQ((X), (Y));                                                      \
  else                                                                        \
  EXPECT_LT(std::abs(((X) - (Y)) /                                            \
                     (std::abs(Y) + std::numeric_limits<double>::epsilon())), \
            EPS)

#endif  // COMPLEXES_TEST_H
