#ifndef ANGMOMCACHE_H_
#define ANGMOMCACHE_H_

#include <cstddef>
#include <vector>

inline std::size_t JJ1BToIndex(int jj1) {
    return (jj1 - 1) / 2;
}

inline std::size_t J2BToIndex(int j_2) {
    return j_2;
}

inline std::size_t JJSToIndex(int jj1, int jj2, int j_3, int jj4, int jj5, int j_6, std::size_t dim_1, std::size_t dim_2) {
    return dim_1 * dim_2 * dim_1 * dim_1 * dim_2 * JJ1BToIndex(jj1)
      +            dim_2 * dim_1 * dim_1 * dim_2 * JJ1BToIndex(jj2)
      +                    dim_1 * dim_1 * dim_2 * J2BToIndex(j_3)
      +                            dim_1 * dim_2 * JJ1BToIndex(jj4)
      +                                    dim_2 * JJ1BToIndex(jj5)
      +                                            J2BToIndex(j_6);
}

class SixJCache_112112 {
  public:
    // Default constructor.
    SixJCache_112112();

    // Constructor from maximum 1-body jj value.
    explicit SixJCache_112112(int jj_1b_max);

    // Get SixJ from cache given jj values for 1-body angular momenta
    // and j values for 2-body angular momenta.
    double SixJ(int jj1, int jj2, int j_3, int jj4, int jj5, int j_6) const {
        std::size_t index = JJSToIndex(jj1, jj2, j_3, jj4, jj5, j_6, dim_1_, dim_2_);
        return six_js_[index];
    }

  private:
    int jj_1_max_ = 1;
    std::size_t dim_1_ = 1;
    std::size_t dim_2_ = 2;
    std::vector<double> six_js_ = std::vector<double>(
        dim_1_ * dim_1_ * dim_2_ * dim_1_ * dim_1_ * dim_2_,
        0.0
    );
};


#endif  // ANGMOMCACHE_H_
