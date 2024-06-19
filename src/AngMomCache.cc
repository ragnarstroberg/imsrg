#include "AngMomCache.hh"

#include <cstddef>
#include <vector>

#include "AngMom.hh"

static inline std::size_t Dim1BFromJJ1BMax(int jj_1b_max) {
    return JJ1BToIndex(jj_1b_max) + 1;
}

static inline std::size_t Dim2BFromJJ1BMax(int jj_1b_max) {
    return J2BToIndex((jj_1b_max + jj_1b_max) / 2) + 1;
}

static std::vector<double> GenerateSixJCache(int jj_1b_max, std::size_t dim_1, std::size_t dim_2);

SixJCache_112112::SixJCache_112112() {}

SixJCache_112112::SixJCache_112112(const int jj_1b_max)
    : jj_1_max_(jj_1b_max),
      dim_1_(Dim1BFromJJ1BMax(jj_1_max_)),
      dim_2_(Dim2BFromJJ1BMax(jj_1_max_)),
      six_js_(GenerateSixJCache(jj_1_max_, dim_1_, dim_2_)) {}


std::vector<double> GenerateSixJCache(const int jj_1b_max, const std::size_t dim_1, const std::size_t dim_2) {
    std::vector<double> six_js(dim_1 * dim_1 * dim_2 * dim_1 * dim_1 * dim_2, 0.0);

    const int j_2b_max = (jj_1b_max + jj_1b_max) / 2;

    for (int jj1 = 1; jj1 <= jj_1b_max; jj1 += 2) {
        double j1d = jj1 * 0.5;
        for (int jj2 = 1; jj2 <= jj_1b_max; jj2 += 2) {
            double j2d = jj2 * 0.5;
            for (int j_3 = 0; j_3 <= j_2b_max; j_3 += 1) {
                double j3d = j_3 * 1.0;
                for (int jj4 = 1; jj4 <= jj_1b_max; jj4 += 2) {
                    double j4d = jj4 * 0.5;
                    for (int jj5 = 1; jj5 <= jj_1b_max; jj5 += 2) {
                        double j5d = jj5 * 0.5;
                        for (int j_6 = 0; j_6 <= j_2b_max; j_6 += 1) {
                            double j6d = j_6 * 1.0;
                            std::size_t index = JJSToIndex(jj1, jj2, j_3, jj4, jj5, j_6, dim_1, dim_2);
                            six_js[index] = AngMom::SixJ(j1d, j2d, j3d, j4d, j5d, j6d);
                        }
                    }
                }
            }
        }
    }

    return six_js;
}
