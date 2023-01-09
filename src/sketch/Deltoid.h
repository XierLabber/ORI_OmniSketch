/**
 * @file Deltoid.h
 * @author deadlycat <lsmfttb@gmail.com> XierLabber<yangshibo@stu.pku.edu.cn>(modified)
 * @brief Implementation of Deltoid
 *
 * @copyright Copyright (c) 2022
 *
 */
#pragma once

#include <common/hash.h>
#include <common/sketch.h>
#include <vector>
#include <algorithm>
#include <cmath>

namespace OmniSketch::Sketch {
/**
 * @brief Deltoid
 *
 * @tparam key_len  length of flowkey
 * @tparam T        type of the counter
 * @tparam hash_t   hashing class
 */
template <int32_t key_len, typename T, 
          typename hash_t = Hash::AwareHash>
class Deltoid : public SketchBase<key_len, T> {
private:
  T sum_;
  int32_t num_hash_;
  int32_t num_group_;
  int32_t nbits_;
  T ***arr1_; // 3d array:
              //  num_hash*num_group*(nbits_+1). Corresponds to T_{a,b,c} in the
              //  paper Here arr1_[a][b][nbits_] is used to represent T_{a, b,
              //  0} in the paper
  T ***arr0_; // 3d array:
              //  num_hash*num_group*(nbits_). Corresponds to T'_{a,b,c} in the
              //  paper
  hash_t *hash_fns_; // hash funcs

public:
  /**
   * @brief Construct by specifying hash number and group number
   *
   */
  Deltoid(int32_t num_hash, int32_t num_group);
  /**
   * @brief Release the pointer
   *
   */
  ~Deltoid();
  /**
   * @brief Update a flowkey with certain value
   *
   */
  void update(const FlowKey<key_len> &flowkey, T val) override;
  /**
   * @brief Query a flowkey
   *
   */
  T query(const FlowKey<key_len> &flowkey) const override;
  /**
   * @brief Get all the heavy hitters
   *
   */
  Data::Estimation<key_len, T> getHeavyHitter(double threshold) const override;  /**
   * @brief Get the size of the sketch
   *
   */
  size_t size() const override;
  /**
   * @brief Reset the sketch
   *
   */
  void clear();
};

} // namespace OmniSketch::Sketch

//-----------------------------------------------------------------------------
//
///                        Implementation of template methods
//
//-----------------------------------------------------------------------------

namespace OmniSketch::Sketch {

template <int32_t key_len, typename T, typename hash_t>
Deltoid<key_len, T, hash_t>::Deltoid(int32_t num_hash, int32_t num_group)
    : num_hash_(num_hash), num_group_(Util::NextPrime(num_group)),
      nbits_(key_len * 8), sum_(0) {

  // allocate continuous memory
  arr1_ = new T **[num_hash_];
  for (int32_t i = 0; i < num_hash_; ++i) {
    arr1_[i] = new T *[num_group_];
  }
  arr1_[0][0] = new T[num_hash_ * num_group_ * (nbits_ + 1)]();
  for (int32_t j = 1; j < num_group_; ++j) {
    arr1_[0][j] = arr1_[0][j - 1] + (nbits_ + 1);
  }
  for (int32_t i = 1; i < num_hash_; ++i) {
    arr1_[i][0] = arr1_[i - 1][0] + num_group_ * (nbits_ + 1);
    for (int32_t j = 1; j < num_group_; ++j) {
      arr1_[i][j] = arr1_[i][j - 1] + (nbits_ + 1);
    }
  }
  arr0_ = new T **[num_hash_];
  for (int32_t i = 0; i < num_hash_; ++i) {
    arr0_[i] = new T *[num_group_];
  }
  arr0_[0][0] = new T[num_hash_ * num_group_ * nbits_]();
  for (int32_t j = 1; j < num_group_; ++j) {
    arr0_[0][j] = arr0_[0][j - 1] + nbits_;
  }
  for (int32_t i = 1; i < num_hash_; ++i) {
    arr0_[i][0] = arr0_[i - 1][0] + num_group_ * nbits_;
    for (int32_t j = 1; j < num_group_; ++j) {
      arr0_[i][j] = arr0_[i][j - 1] + nbits_;
    }
  }
  // hash functions
  hash_fns_ = new hash_t[num_hash_];

  sum_ = 0;
}

template <int32_t key_len, typename T, typename hash_t>
Deltoid<key_len, T, hash_t>::~Deltoid() {
  if (arr1_ != nullptr) {
    delete[] arr1_[0][0];
    for (int32_t i = 0; i < num_hash_; ++i) {
      delete[] arr1_[i];
    }
    delete[] arr1_;
  }
  if (arr0_ != nullptr) {
    delete[] arr0_[0][0];
    for (int32_t i = 0; i < num_hash_; ++i) {
      delete[] arr0_[i];
    }
    delete[] arr0_;
  }
  delete[] hash_fns_;
}

template <int32_t key_len, typename T, typename hash_t>
void Deltoid<key_len, T, hash_t>::update(const FlowKey<key_len> &flowkey,
                                         T val) {
  sum_ += val;
  for (int32_t i = 0; i < num_hash_; ++i) {
    int32_t idx = hash_fns_[i](flowkey) % num_group_;
    for (int32_t j = 0; j < nbits_; ++j) {
      if (flowkey.getBit(j)) {
        arr1_[i][idx][j] += val;
      } else {
        arr0_[i][idx][j] += val;
      }
    }
    arr1_[i][idx][nbits_] += val;
  }
}

template <int32_t key_len, typename T, typename hash_t>
T Deltoid<key_len, T, hash_t>::query(const FlowKey<key_len> &flowkey) const {
  T min_val = std::numeric_limits<T>::max();
  for (int32_t i = 0; i < num_hash_; ++i) {
    int32_t idx = hash_fns_[i](flowkey) % num_group_;
    for (int32_t j = 0; j < nbits_; ++j) {
      if (flowkey.getBit(j)) {
        min_val = std::min(min_val, arr1_[i][idx][j]);
      } else {
        min_val = std::min(min_val, arr0_[i][idx][j]);
      }
    }
  }
  return min_val;
}

template <int32_t key_len, typename T, typename hash_t>
Data::Estimation<key_len, T>
Deltoid<key_len, T, hash_t>::getHeavyHitter(double threshold) const {
  T thresh = threshold;
  double val1 = 0;
  double val0 = 0;
  Data::Estimation<key_len, T> heavy_hitters;
  for (int32_t i = 0; i < num_hash_; i++) {
    for (int32_t j = 0; j < num_group_; j++) {
      if (arr1_[i][j][nbits_] <= thresh) { // no heavy hitter in this group
        continue;
      }
      FlowKey<key_len> fk{}; // create a flowkey with full 0
      bool reject = false;
      for (int32_t k = 0; k < nbits_; k++) {

        bool t1 = (arr1_[i][j][k] > thresh);
        bool t0 = (arr0_[i][j][k] > thresh);
        if (t1 == t0) {
          reject = true;
          break;
        }
        if (t1) {
          fk.setBit(k, true);
        }
      }
      if (reject) {
        continue;
      } else if (!heavy_hitters.count(fk)) {
        T esti_val = query(fk);
        heavy_hitters[fk] = esti_val;
      }
    }
  }
  return heavy_hitters;
}

template <int32_t key_len, typename T, typename hash_t>
size_t Deltoid<key_len, T, hash_t>::size() const {
  return sizeof(Deltoid<key_len, T, hash_t>) +
         (2 * num_group_ * num_hash_ * nbits_ + 1) * sizeof(T) +
         num_hash_ * sizeof(hash_t);
}

template <int32_t key_len, typename T, typename hash_t>
void Deltoid<key_len, T, hash_t>::clear() {
  sum_ = 0;
  std::fill(arr0_[0][0], arr0_[0][0] + num_hash_ * num_group_ * nbits_, 0);
  std::fill(arr1_[0][0], arr1_[0][0] + num_hash_ * num_group_ * (nbits_ + 1),
            0);
}

} // namespace OmniSketch::Sketch