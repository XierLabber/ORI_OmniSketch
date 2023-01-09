/**
 * @file CHNitroSketch
 * @author XierLabber
 * @brief Implementation of CHNitroSketch
 *
 * @copyright Copyright (c) 2022
 *
 */
#pragma once

#include <common/hash.h>
#include <common/hierarchy.h>
#include <common/sketch.h>

#include <algorithm>
#include <cmath>
#include <iostream>
#include <memory>
#include <random>

namespace OmniSketch::Sketch {
/**
 * @brief CHNitroSketch
 *
 * @tparam key_len  length of flowkey
 * @tparam hash_t   hashing class
 * @tparam T        type of the counter
 */

template <int32_t key_len, int32_t no_layer, typename T,
          typename hash_t = Hash::AwareHash>
class CHNitroSketch : public SketchBase<key_len, T> {
public:
  /**
   * @brief Construct by specifying depth and width
   *
   */
  CHNitroSketch(int depth, int width, 
                double cnt_no_ratio,
                const std::vector<size_t> &width_cnt, 
                const std::vector<size_t> &no_hash);
  /**
   * @brief Release the pointer
   *
   */
  ~CHNitroSketch();
  /**
   * @brief Update a flowkey with certain value
   *
   */
  void update(const FlowKey<key_len> &flowkey, T value);

  void alwaysLineRateUpdate(const FlowKey<key_len> &flowkey, T value);

  void alwaysCorrectUpdate(const FlowKey<key_len> &flowkey, T value);

  T query(const FlowKey<key_len> &flowkey) const;

  void adjustUpdateProb(double traffic_rate);

  void clear();
  
  size_t size() const;

private:
  int depth_;
  int width_;

  std::vector<size_t> no_cnt;
  std::vector<size_t> width_cnt;
  std::vector<size_t> no_hash;
  CounterHierarchy<no_layer, T, hash_t> *ch;

  hash_t *hash_fns_;

  int32_t next_bucket_; // index of the next update packet
  int32_t next_packet_; // number of skipped packets

  double update_prob_; // sampling probability
  std::default_random_engine generator;

  bool line_rate_enable_; // enable always line rate
  double switch_thresh_;  // switch to always line rate update

  static double update_probs[8];

  void getNextUpdate(double prob); // update next_bbucket_ and next_packet_
  bool isLineRateUpdate();         // judge whether enable line rate update

  void __do_update(const FlowKey<key_len> &flowkey, T value, double prob);

  size_t get_idx(size_t i, size_t j) const;

  T CHNS_getCnt(size_t idx) const;
};

} // namespace OmniSketch::Sketch

//-----------------------------------------------------------------------------
//
///                        Implementation of templated methods
//
//-----------------------------------------------------------------------------

namespace OmniSketch::Sketch {

template <int32_t key_len, int32_t no_layer, typename T, typename hash_t>
double CHNitroSketch<key_len, no_layer, T, hash_t>::update_probs[8] = {
    1.0, 1.0 / 2, 1.0 / 4, 1.0 / 8, 1.0 / 16, 1.0 / 32, 1.0 / 64, 1.0 / 128};

template <int32_t key_len, int32_t no_layer, typename T, typename hash_t>
size_t CHNitroSketch<key_len, no_layer, T, hash_t>::get_idx(size_t i, size_t j) const{
    return i * width_ + j;
}

template <int32_t key_len, int32_t no_layer, typename T, typename hash_t>
T CHNitroSketch<key_len, no_layer, T, hash_t>::CHNS_getCnt(size_t idx) const{
    return ch->getCnt(idx);
}

template <int32_t key_len, int32_t no_layer, typename T, typename hash_t>
CHNitroSketch<key_len, no_layer, T, hash_t>::CHNitroSketch(
    int depth, int width, 
    double cnt_no_ratio,
    const std::vector<size_t> &width_cnt, 
    const std::vector<size_t> &no_hash)
    : depth_(depth), width_(Util::NextPrime(width)), 
      width_cnt(width_cnt), no_hash(no_hash) {

  // check ratio
  if (cnt_no_ratio <= 0.0 || cnt_no_ratio >= 1.0) {
    throw std::out_of_range("Out of Range: Ratio of #counters of adjacent "
                            "layers in CH should be in (0, 1), but got " +
                            std::to_string(cnt_no_ratio) + " instead.");
  }

  switch_thresh_ = (1.0 + std::sqrt(11.0 / width_)) * width_ * width_;

  line_rate_enable_ = false; // always correct update by default
  update_prob_ = 1.0;

  next_packet_ = 1;
  next_bucket_ = 0;

  hash_fns_ = new hash_t[depth_ * 2];

  // prepare no_cnt
  no_cnt.push_back(static_cast<size_t>(depth_) * width_);
  for (int32_t i = 1; i < no_layer; ++i) {
    size_t last_layer = no_cnt.back();
    no_cnt.push_back(Util::NextPrime(std::ceil(last_layer * cnt_no_ratio)));
  }
  // CH
  ch = new CounterHierarchy<no_layer, T, hash_t>(no_cnt, width_cnt,
                                                 no_hash, true);
}


template <int32_t key_len, int32_t no_layer, typename T, typename hash_t>
CHNitroSketch<key_len, no_layer, T, hash_t>::~CHNitroSketch() {
  delete[] hash_fns_;
  delete[] ch;
}

template <int32_t key_len, int32_t no_layer, typename T, typename hash_t>
void CHNitroSketch<key_len, no_layer, T, hash_t>::alwaysLineRateUpdate(
    const FlowKey<key_len> &flowkey, T value) {
  __do_update(flowkey, value, update_prob_);
}

template <int32_t key_len, int32_t no_layer, typename T, typename hash_t>
void CHNitroSketch<key_len, no_layer, T, hash_t>::alwaysCorrectUpdate(
    const FlowKey<key_len> &flowkey, T value) {
  if (isLineRateUpdate()) {
    __do_update(flowkey, value, update_prob_);
  } else {
    __do_update(flowkey, value, 1.0);
  }
}

template <int32_t key_len, int32_t no_layer, typename T, typename hash_t>
std::size_t CHNitroSketch<key_len, no_layer, T, hash_t>::size() const{
  ch->print_rate("");
  printf("\n");
  return sizeof(CHNitroSketch<key_len, no_layer, T, hash_t>) 
        + depth_ * 2 * sizeof(hash_t) 
        + ch->size();
}

template <int32_t key_len, int32_t no_layer, typename T, typename hash_t>
void CHNitroSketch<key_len, no_layer, T, hash_t>::clear() {
  ch->clear();
}

template <int32_t key_len, int32_t no_layer, typename T, typename hash_t>
void CHNitroSketch<key_len, no_layer, T, hash_t>::__do_update(const FlowKey<key_len> &flowkey,
                                         T value, double prob) {
  next_packet_--; // skip packets
  if (next_packet_ == 0) {
    int i;
    for (;;) {
      i = next_bucket_;
      int index = hash_fns_[i](flowkey) % width_;

      double delta =
          1.0 * value / prob *
          (2 * static_cast<int>(hash_fns_[depth_ + i](flowkey) & 1) - 1);
      ch->updateCnt(get_idx(i, index), -static_cast<T>(delta));

      getNextUpdate(prob);
      if (next_packet_ > 0) {
        break;
      }
    }
  }
}

template <int32_t key_len, int32_t no_layer, typename T, typename hash_t>
T CHNitroSketch<key_len, no_layer, T, hash_t>::query(const FlowKey<key_len> &flowkey) const{

  static bool my_cnt_distrib = true;
  if (my_cnt_distrib) {
      std::vector<double> distrib(32);
      for (int32_t i = 0; i < no_cnt[0]; ++i) {
      T val = ch->getOriginalCnt(i);
      for (int32_t k = 0; k < 32; ++k) {
          if (std::abs(val) >= (1 << k))
          distrib[k] += 1.0;
          else
          break;
      }
      }
      for (int32_t k = 0; k < 32; ++k) {
      std::cout << k << ": " << distrib[k] / no_cnt[0] << " \n";
      if (distrib[k] == 0.0) {
          std::cout << std::endl;
          break;
      }
      }
      my_cnt_distrib = false;
  }
  
  T median;
  T values[depth_];
  for (int i = 0; i < depth_; i++) {
    int index = hash_fns_[i](flowkey) % width_;
    values[i] = CHNS_getCnt(get_idx(i, index)) *
                (2 * static_cast<int>(hash_fns_[depth_ + i](flowkey) & 1) - 1);
  }
  std::sort(values, values + depth_);
  if (depth_ & 1) {
    median = std::abs(values[depth_ / 2]);
  } else {
    median = std::abs((values[depth_ / 2 - 1] + values[depth_ / 2]) / 2);
  }
  return median;
}

template <int32_t key_len, int32_t no_layer, typename T, typename hash_t>
void CHNitroSketch<key_len, no_layer, T, hash_t>::getNextUpdate(double prob) {
  int sample = 1;
  if (prob < 1.0) {
    std::geometric_distribution<int> dist(prob);
    sample = 1 + dist(generator);
  }
  next_bucket_ = next_bucket_ + sample;
  next_packet_ = ((int)(next_bucket_ / depth_));
  next_bucket_ %= depth_;
}

template <int32_t key_len, int32_t no_layer, typename T, typename hash_t>
bool CHNitroSketch<key_len, no_layer, T, hash_t>::isLineRateUpdate() {
  if (line_rate_enable_) {
    return true;
  } else {
    double median;
    double values[depth_];
    for (int i = 0; i < depth_; i++) {
      values[i] =0.0;
      for(int j = 0; j < width_; j++)
      {
        double tmp = CHNS_getCnt(get_idx(i, j));
        values[i] += tmp * tmp;
      }
    }
    std::sort(values, values + depth_);
    if (depth_ & 1) {
      median = values[depth_ / 2];
    } else {
      median = (values[depth_ / 2 - 1] + values[depth_ / 2]) / 2;
    }
    if (median >= switch_thresh_) {
      std::cout << "line rate update enable\n";
      line_rate_enable_ = true;
    }
    return line_rate_enable_;
  }
}

template <int32_t key_len, int32_t no_layer, typename T, typename hash_t>
void CHNitroSketch<key_len, no_layer, T, hash_t>::adjustUpdateProb(double traffic_rate) {
  int log_rate = static_cast<int>(std::log2(traffic_rate));
  int update_index = std::max(0, std::min(log_rate, 7));
  update_prob_ = update_probs[update_index];
}

template <int32_t key_len, int32_t no_layer, typename T, typename hash_t>
void CHNitroSketch<key_len, no_layer, T, hash_t>::update(
    const FlowKey<key_len> &flowkey, T value) {
  __do_update(flowkey, value, update_prob_);
}

} // namespace OmniSketch::Sketch
