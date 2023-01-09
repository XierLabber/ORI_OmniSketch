/**
 * @file LDSketch.h
 * @author Ferric XierLabber (you@domain.com)
 * @brief Implementation of LD Sketch
 *
 * @copyright Copyright (c) 2022
 *
 */

#define RECORD_PACKET_NUM
#define RECORD_MAP_MEM

#pragma once

#include <common/hash.h>
#include <common/sketch.h>

namespace OmniSketch::Sketch {
/**
 * @brief LD Sketch
 *
 * @tparam key_len  length of flowkey
 * @tparam T        type of the counter
 * @tparam hash_t   hashing class
 */
template <int32_t key_len, typename T, typename hash_t = Hash::AwareHash>
class LDSketch : public SketchBase<key_len, T> {
private:
  int32_t depth_;
  int32_t width_;
  double expansion_;
  double thre;

  hash_t *hash_fns_;

  struct Bounds {
    T lower, upper;
  };

  struct Bucket {
    T V, e;
    int32_t l;
    std::map<FlowKey<key_len>, T> A;

    Bucket();
    void update(const FlowKey<key_len> &flowkey, T val, double expansion);
    Bounds query(const FlowKey<key_len> &flowkey) const;
    void clear();
    size_t size() const;
  };
  Bucket **counter_;

  LDSketch(const LDSketch &) = delete;
  LDSketch(LDSketch &&) = delete;

public:
  /**
   * @brief Construct by specifying depth and width
   *
   */
  LDSketch(int32_t depth, int32_t width, 
           double eps, double threshold);
  /**
   * @brief Release the pointer
   *
   */
  ~LDSketch();
  /**
   * @brief Update a flowkey with certain value
   *
   */
  void update(const FlowKey<key_len> &flowkey, T val) override;
  /**
   * @brief Query a flowkey
   *
   */
  Data::Estimation<key_len, T> getHeavyHitter(double val_threshold) const;
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
};

} // namespace OmniSketch::Sketch

//-----------------------------------------------------------------------------
//
///                        Implementation of templated methods
//
//-----------------------------------------------------------------------------

namespace OmniSketch::Sketch {


template <int32_t key_len, typename T, typename hash_t>
LDSketch<key_len, T, hash_t>::Bucket::Bucket() : V(0), e(0), l(0) {}

template <int32_t key_len, typename T, typename hash_t>
void LDSketch<key_len, T, hash_t>::Bucket::update(
    const FlowKey<key_len> &flowkey, T val, double expansion) {

#ifdef RECORD_PACKET_NUM
  static int tmp = 0;
  tmp++;
  if(tmp % 100000 == 0){
    printf("FINISH UPDATING %d PACKETS\n", tmp);
  }
#endif
  
  V += val;

  auto it = A.find(flowkey);
  if (it != A.end())
    it->second += val;
  else if (A.size() < l)
    A.emplace(flowkey, val);
  else {
    T k = V / expansion;
    if ((k + 1) * (k + 2) - 1 <= l) {
      T e_cap = val;
      for (const auto &kv : A)
        e_cap = std::min(e_cap, kv.second);
      e += e_cap;
      for (auto it = A.begin(); it != A.end();) {
        it->second -= e_cap;
        if (it->second <= 0)
          it = A.erase(it);
        else
          ++it;
      }
      if (val > e_cap)
        A.emplace(flowkey, val - e_cap);
    } else {
      l = (k + 1) * (k + 2) - 1;
      A.emplace(flowkey, val);
    }
  }
}

template <int32_t key_len, typename T, typename hash_t>
typename LDSketch<key_len, T, hash_t>::Bounds
LDSketch<key_len, T, hash_t>::Bucket::query(
    const FlowKey<key_len> &flowkey) const {
  auto it = A.find(flowkey);
  if (it == A.end())
    return {0, e};
  else
    return {it->second, it->second + e};
}

template <int32_t key_len, typename T, typename hash_t>
void LDSketch<key_len, T, hash_t>::Bucket::clear() {
  V = e = 0;
  l = 0;
  A.clear();
}

template <int32_t key_len, typename T, typename hash_t>
size_t LDSketch<key_len, T, hash_t>::Bucket::size() const{
  // Not accurate
  return sizeof(LDSketch<key_len, T, hash_t>::Bucket) +
         A.size() * (key_len + sizeof(T));
}

template <int32_t key_len, typename T, typename hash_t>
LDSketch<key_len, T, hash_t>::LDSketch(int32_t depth, int32_t width, 
                                       double eps, double threshold)
    : depth_(depth), width_(Util::NextPrime(width)), 
      expansion_(eps * threshold), thre(threshold) {

  hash_fns_ = new hash_t[depth_];

  // Allocate continuous memory;
  counter_ = new Bucket *[depth_];
  counter_[0] = new Bucket[depth_ * width_]();
  for (int i = 1; i < depth_; ++i)
    counter_[i] = counter_[i - 1] + width_;
}

template <int32_t key_len, typename T, typename hash_t>
LDSketch<key_len, T, hash_t>::~LDSketch() {
  delete[] hash_fns_;

  delete[] counter_[0];
  delete[] counter_;
}

template <int32_t key_len, typename T, typename hash_t>
void LDSketch<key_len, T, hash_t>::update(const FlowKey<key_len> &flowkey,
                                          T val) {
  for (int i = 0; i < depth_; ++i) {
    int32_t index = hash_fns_[i](flowkey) % width_;
    counter_[i][index].update(flowkey, val, expansion_);
  }
}

template <int32_t key_len, typename T, typename hash_t>
Data::Estimation<key_len, T> LDSketch<key_len, T, hash_t>::getHeavyHitter(
  double val_threshold) const{

  assert((thre > val_threshold - 1e-6) && 
         (thre < val_threshold + 1e-6));

  Data::Estimation<key_len, T> heavy_hitters;

  for (int i = 0; i < depth_; ++i)
    for (int j = 0; j < width_; ++j) {
      const Bucket &bucket = counter_[i][j];
      if (bucket.V < val_threshold)
        continue;

      for (const auto &kv : bucket.A) {
        if (heavy_hitters.count(kv.first))
          continue;

        std::vector<T> upper_bounds(depth_);
        for (int k = 0; k < depth_; ++k) {
          int32_t index = hash_fns_[k](kv.first) % width_;
          upper_bounds[k] = counter_[k][index].query(kv.first).upper;
        }

        T u = *std::min_element(upper_bounds.begin(), upper_bounds.end());
        if (u >= val_threshold)
          heavy_hitters[kv.first] = u;
      }
    }

  return heavy_hitters;
  
}

template <int32_t key_len, typename T, typename hash_t>
size_t LDSketch<key_len, T, hash_t>::size() const {
  #ifdef RECORD_MAP_MEM
  size_t map_mem = 0;
  #endif
  size_t s = sizeof(LDSketch<key_len, T, hash_t>) + // Instance
             sizeof(hash_t) * depth_ +              // hash_fns
             sizeof(Bucket *) * depth_;
  for (int i = 0; i < depth_; ++i){
    for (int j = 0; j < width_; ++j){
      #ifndef RECORD_MAP_MEM
      s += counter_[i][j].size(); // Respective sizes
      #else
      size_t tmp = counter_[i][j].size();
      s += tmp;
      map_mem += tmp;
      #endif
    }
  }
  #ifdef RECORD_MAP_MEM
  printf("MAP MEM: %ld BYTES, %lf KB, %lf MB\n", 
          map_mem, (double)map_mem / 1024, (double)map_mem / (1024 * 1024));
  #endif
  return s;
}

template <int32_t key_len, typename T, typename hash_t>
void LDSketch<key_len, T, hash_t>::clear() {
  for (int i = 0; i < depth_; ++i)
    for (int j = 0; j < width_; ++j)
      counter_[i][j].clear();
}

} // namespace OmniSketch::Sketch
