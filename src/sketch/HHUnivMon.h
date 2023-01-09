/**
 * @file HHUnivMon.h
 * @author dromniscience XierLabber (you@domain.com)
 * @brief Implementation of UnivMon, which is used to
 *        get heavy hitters 
 *
 * @copyright Copyright (c) 2022
 *
 */
#pragma once
#include <iostream>

#include <common/hash.h>
#include <common/StreamSummary.h>
#include <common/sketch.h>
#include <sketch/CountSketch.h>

// #define NO_HEAP_SIZE

namespace OmniSketch::Sketch {
/**
 * @brief HHUnivMon
 *
 * @tparam key_len  length of flowkey
 * @tparam T        type of the counter
 * @tparam hash_t   hashing class
 */
template <int32_t key_len, typename T, typename hash_t = Hash::AwareHash>
class HHUnivMon : public SketchBase<key_len, T> {
private:
  int32_t depth;
  int32_t width;
  int32_t logn;
  hash_t *hash_fns;
  CountSketch<key_len, T, hash_t> **sketch;
  Data::Estimation<key_len> *flows;

  HHUnivMon(const HHUnivMon &) = delete;
  HHUnivMon(HHUnivMon &&) = delete;

  T* sum;
  StreamSummary<key_len, T, hash_t>* HHHeaps;

public:
  /**
   * @brief Construct by specifying depth, width and $\log n$, where $n$ is the
   * number of flows to insert.
   *
   */
  HHUnivMon(int32_t depth_, int32_t width_, int32_t log_n, 
            int32_t heap_size, double SSalpha);
  /**
   * @brief Release the pointer
   *
   */
  ~HHUnivMon();
  /**
   * @brief Update a flowkey with certain value
   *
   */
  void update(const FlowKey<key_len> &flowkey, T val) override;
  /**
   * @brief Get Heavy Hitter
   *
   */
  Data::Estimation<key_len, T> getHeavyHitter(double val_threshold) const;
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
  int32_t getHash(const FlowKey<key_len> &flowkey, int32_t layer) const;
};

} // namespace OmniSketch::Sketch

//-----------------------------------------------------------------------------
//
///                        Implementation of templated methods
//
//-----------------------------------------------------------------------------

namespace OmniSketch::Sketch {

template <int32_t key_len, typename T, typename hash_t>
int32_t HHUnivMon<key_len, T, hash_t>::getHash(const FlowKey<key_len> &flowkey, int32_t layer) const{
  if(layer == 0){
    return 1;
  }
  return (hash_fns[layer - 1](flowkey) & 1);
}

template <int32_t key_len, typename T, typename hash_t>
HHUnivMon<key_len, T, hash_t>::HHUnivMon(int32_t depth_, int32_t width_,
                                     int32_t log_n, int32_t heap_size, double SSalpha)
    : depth(depth_), width(Util::NextPrime(width_)), logn(log_n) {

  hash_fns = new hash_t[logn - 1];
  // Allocate continuous memory
  sketch = new CountSketch<key_len, T, hash_t> *[logn];
  for (int32_t i = 0; i < logn; ++i){
    sketch[i] = new CountSketch<key_len, T, hash_t>(depth, width_);
    width_ = std::max(1, width_ / 2);
  }
  flows = new Data::Estimation<key_len>[logn];

  HHHeaps = new StreamSummary<key_len, T, hash_t>[logn];
  for(int i = 0; i < logn; i++){
    HHHeaps[i].init(heap_size, SSalpha);
  }
  sum = new T[logn];
}

template <int32_t key_len, typename T, typename hash_t>
HHUnivMon<key_len, T, hash_t>::~HHUnivMon() {
  if (hash_fns)
    delete[] hash_fns;
  if (sketch) {
    for (int32_t i = 0; i < logn; ++i)
      delete sketch[i];
    delete[] sketch;
  }
  if (flows)
    delete[] flows;
  for(int32_t i = 0; i < logn; i++){
    HHHeaps[i].destroy();
  }
  delete[] HHHeaps;
  delete[] sum;
}

template <int32_t key_len, typename T, typename hash_t>
Data::Estimation<key_len, T> HHUnivMon<key_len, T, hash_t>::getHeavyHitter(double val_threshold) const{
  Data::Estimation<key_len, T> ret;
  for(int i = 0; i < logn; i++){
    Data::Estimation<key_len, T> tmp = HHHeaps[i].getHeavyHitter(val_threshold);
    for(auto kv : tmp){
      if(!ret.count(kv.get_left())){
        ret[kv.get_left()] = kv.get_right();
      }
    }
  }
  return ret;
}

template <int32_t key_len, typename T, typename hash_t>
void HHUnivMon<key_len, T, hash_t>::update(const FlowKey<key_len> &flowkey,
                                         T val) {
  for (int32_t i = 0; i < logn; ++i) {
    if (getHash(flowkey, i)) {
      sum[i] += val;
      sketch[i]->update(flowkey, val);
      T est = sketch[i]->query(flowkey);
      HHHeaps[i].insert(flowkey, est);
    } else
      break;
  }
}

template <int32_t key_len, typename T, typename hash_t>
T HHUnivMon<key_len, T, hash_t>::query(const FlowKey<key_len> &flowkey) const {
  int level;
  for(level = 0; level < logn; level++){
    if(!getHash(flowkey, level)){
      break;
    }
  }
  level--;
  T ret = sketch[level]->query(flowkey);
  for(int i = level - 1; i >= 0; i--){
    ret = 2 * ret - sketch[i]->query(flowkey);
  }
  return ret;
}

template <int32_t key_len, typename T, typename hash_t>
size_t HHUnivMon<key_len, T, hash_t>::size() const {
  size_t total = sizeof(*this); // instance
  for (int32_t i = 0; i < logn; ++i)
    total += sketch[i]->size(); // L2 HH
  total += sizeof(hash_t) * (logn - 1);
  #ifndef NO_HEAP_SIZE
  for(int32_t i = 0; i < logn; i++){
    total += HHHeaps[i].memory_size();
  }
  #endif
  total += sizeof(T) * logn;
  return total;
}

template <int32_t key_len, typename T, typename hash_t>
void HHUnivMon<key_len, T, hash_t>::clear() {
  // std::fill(counter[0], counter[0] + depth * width, 0);
}

} // namespace OmniSketch::Sketch


