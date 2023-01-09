/**
 * @file CHDeltoid2Tuple.h
 * @author XierLabber<yangshibo@stu.pku.edu.cn>
 * @brief Implementation of CHDeltoid2Tuple
 *
 * @copyright Copyright (c) 2022
 *
 */
#pragma once

#include <common/hash.h>
#include <common/hierarchy.h>
#include <common/sketch.h>
#include <vector>
#include <algorithm>
#include <cmath>

// #define MY_DEBUG

namespace OmniSketch::Sketch {
/**
 * @brief Deltoid with CH
 *
 * @tparam key_len  length of flowkey
 * @tparam no_layer layer of CH
 * @tparam T        type of the counter
 * @tparam hash_t   hashing class
 */
template <int32_t key_len, int32_t no_layer, typename T, 
          typename hash_t = Hash::AwareHash>
class CHDeltoid2Tuple : public SketchBase<key_len, T> {
private:

#ifdef MY_DEBUG
  int flow_cnt;
#endif

  T sum_;
  int32_t num_hash_;
  int32_t num_group_;
  int32_t nbits_;

  std::vector<size_t> no_cnt1_;
  std::vector<size_t> no_cnt0_;
  std::vector<size_t> width_cnt;
  std::vector<size_t> no_hash;

  hash_t *hash_fns_; // hash funcs
  CounterHierarchy<no_layer,T,hash_t> *ch1_;    // 3d array:
                                                //  num_hash*num_group*(nbits_+1). Corresponds to T_{a,b,c} in the
                                                //  paper Here arr1_[a][b][nbits_] is used to represent T_{a, b,
                                                //  0} in the paper
  CounterHierarchy<no_layer,T,hash_t> *ch0_;    // 3d array:
                                                //  num_hash*num_group*(nbits_). Corresponds to T'_{a,b,c} in the
                                                //  paper

  CHDeltoid2Tuple(const CHDeltoid2Tuple &) = delete;
  CHDeltoid2Tuple(CHDeltoid2Tuple &&) = delete;

public:
  /**
   * @brief Construct by specifying hash number and group number
   *
   */
  CHDeltoid2Tuple(int32_t num_hash, int32_t num_group, double cnt_no_ratio,
            const std::vector<size_t> &width_cnt,
            const std::vector<size_t> &no_hash);
  /**
   * @brief Release the pointer
   *
   */
  ~CHDeltoid2Tuple();
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
  /**
   * @brief calculate the corrosponding index  
   */
  size_t get_idx1_(size_t i, size_t j, size_t k) const;
  size_t get_idx0_(size_t i, size_t j, size_t k) const;
};

} // namespace OmniSketch::Sketch

//-----------------------------------------------------------------------------
//
///                        Implementation of template methods
//
//-----------------------------------------------------------------------------

namespace OmniSketch::Sketch {

template <int32_t key_len, int32_t no_layer, typename T, typename hash_t>
size_t CHDeltoid2Tuple<key_len, no_layer, T, hash_t>::get_idx1_(size_t i, size_t j, size_t k) const{
    return i * num_group_ * (nbits_ + 1) + j * (nbits_ + 1) + k;
}

template <int32_t key_len, int32_t no_layer, typename T, typename hash_t>
size_t CHDeltoid2Tuple<key_len, no_layer, T, hash_t>::get_idx0_(size_t i, size_t j, size_t k) const{
    return i * num_group_ * (nbits_) + j * (nbits_) + k;
}

template <int32_t key_len, int32_t no_layer, typename T, typename hash_t>
CHDeltoid2Tuple<key_len, no_layer, T, hash_t>::CHDeltoid2Tuple(
    int32_t num_hash, int32_t num_group, 
    double cnt_no_ratio,
    const std::vector<size_t> &width_cnt,
    const std::vector<size_t> &no_hash)
    : num_hash_(num_hash), num_group_(Util::NextPrime(num_group)),
      nbits_(key_len * 8), sum_(0), 
      ch1_(nullptr), ch0_(nullptr),
      width_cnt(width_cnt), no_hash(no_hash){
#ifdef MY_DEBUG
  flow_cnt = 0;
#endif

  // check ratio
  if (cnt_no_ratio <= 0.0 || cnt_no_ratio >= 1.0) {
    throw std::out_of_range("Out of Range: Ratio of #counters of adjacent "
                            "layers in CH should be in (0, 1), but got " +
                            std::to_string(cnt_no_ratio) + " instead.");
  }

  // hash functions
  hash_fns_ = new hash_t[num_hash_];

  sum_ = 0;

  // prepare no_cnt
  no_cnt1_.push_back(static_cast<size_t>(num_hash_ * num_group_ * (nbits_ + 1)));
  no_cnt0_.push_back(static_cast<size_t>(num_hash_ * num_group_ * (nbits_)));
  for (int32_t i = 1; i < no_layer; ++i) {
    size_t last_layer = no_cnt1_.back();
    no_cnt1_.push_back(Util::NextPrime(std::ceil(last_layer * cnt_no_ratio)));
  }
  for (int32_t i = 1; i < no_layer; ++i) {
    size_t last_layer = no_cnt0_.back();
    no_cnt0_.push_back(Util::NextPrime(std::ceil(last_layer * cnt_no_ratio)));
  }

  ch1_ = new CounterHierarchy<no_layer, T, hash_t>(no_cnt1_, width_cnt, no_hash);
  ch0_ = new CounterHierarchy<no_layer, T, hash_t>(no_cnt0_, width_cnt, no_hash);
}

template <int32_t key_len, int32_t no_layer, typename T, typename hash_t>
CHDeltoid2Tuple<key_len, no_layer, T, hash_t>::~CHDeltoid2Tuple() {
  delete[] hash_fns_;
  delete[] ch1_;
  delete[] ch0_;
}

template <int32_t key_len, int32_t no_layer, typename T, typename hash_t>
void CHDeltoid2Tuple<key_len, no_layer, T, hash_t>::update(const FlowKey<key_len> &flowkey,
                                         T val) {
#ifdef MY_DEBUG
  if(flow_cnt % 10000 == 0)
  {
    printf("FINISH UPDATING %d FLOWS!\n", flow_cnt);
  }
  flow_cnt++;
#endif
  sum_ += val;
  for (int32_t i = 0; i < num_hash_; ++i) {
    int32_t idx = hash_fns_[i](flowkey) % num_group_;
    for (int32_t j = 0; j < nbits_; ++j) {
      if (flowkey.getBit(j)) {
        ch1_->updateCnt(get_idx1_(i, idx, j), val);
      } else {
        ch0_->updateCnt(get_idx0_(i, idx, j), val);
      }
    }
    ch1_->updateCnt(get_idx1_(i, idx, nbits_), val);
  }
}

template <int32_t key_len, int32_t no_layer, typename T, typename hash_t>
T CHDeltoid2Tuple<key_len, no_layer, T, hash_t>::query(const FlowKey<key_len> &flowkey) const {  
  
  T min_val = std::numeric_limits<T>::max();
  for (int32_t i = 0; i < num_hash_; ++i) {
    int32_t idx = hash_fns_[i](flowkey) % num_group_;
    for (int32_t j = 0; j < nbits_; ++j) {
      if (flowkey.getBit(j)) {
        min_val = std::min(min_val, ch1_->getCnt(get_idx1_(i, idx, j)));
      } else {
        min_val = std::min(min_val, ch0_->getCnt(get_idx0_(i, idx, j)));
      }
    }
  }
  return min_val;
}

template <int32_t key_len, int32_t no_layer, typename T, typename hash_t>
Data::Estimation<key_len, T>
CHDeltoid2Tuple<key_len, no_layer, T, hash_t>::getHeavyHitter(double threshold) const {
  T thresh = threshold;
  double val1 = 0;
  double val0 = 0;
  Data::Estimation<key_len, T> heavy_hitters;
  for (int32_t i = 0; i < num_hash_; i++) {
    for (int32_t j = 0; j < num_group_; j++) {
      if (ch1_->getCnt(get_idx1_(i, j, nbits_)) <= thresh) { // no heavy hitter in this group
        continue;
      }
      FlowKey<key_len> fk{}; // create a flowkey with full 0
      bool reject = false;
      for (int32_t k = 0; k < nbits_; k++) {

        bool t1 = (ch1_->getCnt(get_idx1_(i, j, k)) > thresh);
        bool t0 = (ch0_->getCnt(get_idx0_(i, j, k)) > thresh);
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

template <int32_t key_len, int32_t no_layer, typename T, typename hash_t>
size_t CHDeltoid2Tuple<key_len, no_layer, T, hash_t>::size() const {
  ch1_->print_rate("CH1");
  ch0_->print_rate("CH0");
  return sizeof(CHDeltoid2Tuple<key_len, no_layer, T, hash_t>) +
         num_hash_ * sizeof(hash_t) +
         ch1_->size() + ch0_->size();
}

template <int32_t key_len, int32_t no_layer, typename T, typename hash_t>
void CHDeltoid2Tuple<key_len, no_layer, T, hash_t>::clear() {
  sum_ = 0;
  ch1_->clear();
  ch0_->clear();
}

} // namespace OmniSketch::Sketch