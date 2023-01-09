/**
 * @file CHCUSketch.h
 * @author XierLabber (you@domain.com)
 * @brief
 *
 *
 */
#pragma once

#include <common/hash.h>
#include <common/hierarchy.h>
#include <common/sketch.h>

namespace OmniSketch::Sketch {
/**
 * @brief CHCU Sketch
 *
 * @tparam key_len  length of flowkey
 * @tparam T        type of the counter
 * @tparam hash_t   hashing class
 */
template <int32_t key_len, int32_t no_layer, typename T,
          typename hash_t>
class CHCUSketch : public SketchBase<key_len, T> {
private:
  int32_t depth;
  int32_t width;
  hash_t *hash_fns;

  std::vector<size_t> no_cnt;
  std::vector<size_t> width_cnt;
  std::vector<size_t> no_hash;

  CounterHierarchy<no_layer, T, hash_t> *ch;

  CHCUSketch(const CHCUSketch &) = delete;
  CHCUSketch(CHCUSketch &&) = delete;

public:
  /**
   * @brief Construct by specifying depth and width
   *
   */
  CHCUSketch(int32_t depth_, int32_t width_, double cnt_no_ratio_,
             const std::vector<size_t> &width_cnt_,
             const std::vector<size_t> &no_hash_, 
             const int32_t ch_cm_r_, 
             const int32_t ch_cm_w_);
  /**
   * @brief Release the pointer
   *
   */
  ~CHCUSketch();
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
   * @brief Get the size of the sketch
   *
   */
  size_t size() const override;
  /**
   * @brief Reset the sketch
   *
   */
  void clear();

  int32_t getCHIdx(int32_t i, int32_t j) const{
    return i * width + j;
  }
};

} // namespace OmniSketch::Sketch

//-----------------------------------------------------------------------------
//
///                        Implementation of templated methods
//
//-----------------------------------------------------------------------------

namespace OmniSketch::Sketch {
template <int32_t key_len, int32_t no_layer, typename T,
          typename hash_t>
CHCUSketch<key_len, no_layer, T, hash_t>::CHCUSketch(int32_t depth_, int32_t width_,
             double cnt_no_ratio_,
             const std::vector<size_t> &width_cnt_,
             const std::vector<size_t> &no_hash_, 
             const int32_t ch_cm_r_, 
             const int32_t ch_cm_w_)
    : depth(depth_), width(Util::NextPrime(width_)),
    width_cnt(width_cnt_), no_hash(no_hash_) {
  hash_fns = new hash_t[depth];
  // check ratio
  if (cnt_no_ratio_ <= 0.0 || cnt_no_ratio_ >= 1.0) {
    throw std::out_of_range("Out of Range: Ratio of #counters of adjacent "
                            "layers in CH should be in (0, 1), but got " +
                            std::to_string(cnt_no_ratio_) + " instead.");
  }
  // prepare no_cnt
  no_cnt.push_back(static_cast<size_t>(this->depth) * this->width);
  for (int32_t i = 1; i < no_layer; ++i) {
    size_t last_layer = no_cnt.back();
    no_cnt.push_back(Util::NextPrime(std::ceil(last_layer * cnt_no_ratio_)));
  }
  // CH
  ch = new CounterHierarchy<no_layer, T, hash_t>(no_cnt, this->width_cnt,
                                                 this->no_hash, false, true, 
                                                 ch_cm_r_, ch_cm_w_);
}

template <int32_t key_len, int32_t no_layer, typename T,
          typename hash_t>
CHCUSketch<key_len, no_layer, T, hash_t>::~CHCUSketch() {
  delete[] hash_fns;
  delete[] ch;
}

template <int32_t key_len, int32_t no_layer, typename T,
          typename hash_t>
void CHCUSketch<key_len, no_layer, T, hash_t>::update(const FlowKey<key_len> &flowkey,
                                          T val) {
  int32_t indices[depth];
  T EstVal[depth];
  T min_val = std::numeric_limits<T>::max();
  for (int32_t i = 0; i < depth; ++i) {
    int32_t idx = hash_fns[i](flowkey) % width;
    int32_t CHIdx = getCHIdx(i, idx);
    indices[i] = CHIdx;

    int32_t tmp = ch->getEstCnt(CHIdx);
    EstVal[i] = tmp;
    min_val = std::min(min_val, tmp);
  }
  min_val += val;
  for (int32_t i = 0; i < depth; ++i) {
    if(EstVal[i] < min_val){
      ch->updateCnt(indices[i], min_val - EstVal[i]);
    }
  }
}

template <int32_t key_len, int32_t no_layer, typename T,
          typename hash_t>
T CHCUSketch<key_len, no_layer, T, hash_t>::query(const FlowKey<key_len> &flowkey) const {
  T min_val = std::numeric_limits<T>::max();
  for (int32_t i = 0; i < depth; ++i) {
    int32_t index = hash_fns[i](flowkey) % width;
    int32_t chIdx = getCHIdx(i, index);
    min_val = std::min(min_val, ch->getCnt(chIdx));
  }
  return min_val;
}

template <int32_t key_len, int32_t no_layer, typename T,
          typename hash_t>
size_t CHCUSketch<key_len, no_layer, T, hash_t>::size() const {
  ch->print_rate("CH for CU");
  return sizeof(*this)                // instance
         + sizeof(hash_t) * depth     // hashing class
         + ch->size();
}

template <int32_t key_len, int32_t no_layer, typename T,
          typename hash_t>
void CHCUSketch<key_len, no_layer, T, hash_t>::clear() {
  ch->clear();
}

} // namespace OmniSketch::Sketch