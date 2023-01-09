/**
 * @file MVSketch.h
 * @author 
 * @brief Implementation of Count Min Sketch
 *
 * @copyright Copyright (c) 2022
 *
 */

#include <algorithm>
#include <map>
#include <vector>

#include <common/hash.h>
#include <common/sketch.h>
#include <common/utils.h>

namespace OmniSketch::Sketch {
template <int32_t key_len, typename T, typename hash_t> 
class MVSketch : public SketchBase<key_len, T> {
  int32_t depth_;
  int32_t width_;

  hash_t *hash_fns_;

  struct Bounds {
    T lower;
    T upper;
  };

  struct Bucket {
    T V;
    FlowKey<key_len> K;
    T C;
  };

  Bucket **counter_;

public:
  MVSketch(int32_t depth, int32_t width);
  MVSketch(MVSketch &&) = delete;
  ~MVSketch();
  MVSketch &operator=(const MVSketch &) = delete;
  MVSketch &operator=(MVSketch &&) = delete;

  void update(const FlowKey<key_len> &flow_key, T val);
  T query(const FlowKey<key_len> &flow_key) const;

  void clear();
  size_t size() const;

  Bounds queryBounds(const FlowKey<key_len> &flow_key) const;
  Data::Estimation<key_len, T>  getHeavyHitter(double threshold) const;
};

} // namespace OmniSketch::Sketch

//-----------------------------------------------------------------------------
//
///                        Implementation of templated methods
//
//-----------------------------------------------------------------------------

namespace OmniSketch::Sketch {

template <int32_t key_len, typename T, typename hash_t>
MVSketch<key_len, T, hash_t>::MVSketch(int32_t depth, int32_t width)
    : depth_(depth), width_(Util::NextPrime(width)) {
  hash_fns_ = new hash_t[depth_];

  // Allocate continuous memory
  counter_ = new Bucket *[depth_];
  counter_[0] = new Bucket[depth_ * width_]();
  for (int i = 1; i < depth_; ++i)
    counter_[i] = counter_[i - 1] + width_;
}

template <int32_t key_len, typename T, typename hash_t>
MVSketch<key_len, T, hash_t>::~MVSketch() {
  delete[] hash_fns_;

  delete[] counter_[0];
  delete[] counter_;
}

template <int32_t key_len, typename T, typename hash_t>
void MVSketch<key_len, T, hash_t>::update(const FlowKey<key_len> &flow_key,
                                          T val) {
  for (int i = 0; i < depth_; ++i) {
    int index = hash_fns_[i](flow_key) % width_;
    counter_[i][index].V += val;
    if (counter_[i][index].K == flow_key)
      counter_[i][index].C += val;
    else {
      counter_[i][index].C -= val;
      if (counter_[i][index].C < 0) {
        counter_[i][index].K = flow_key;
        counter_[i][index].C *= -1;
      }
    }
  }
}

template <int32_t key_len, typename T, typename hash_t>
T MVSketch<key_len, T, hash_t>::query(const FlowKey<key_len> &flow_key) const {
  std::vector<T> S_cap(depth_);

  for (int i = 0; i < depth_; ++i) {
    int index = hash_fns_[i](flow_key) % width_;
    if (counter_[i][index].K == flow_key)
      S_cap[i] = (counter_[i][index].V + counter_[i][index].C) / 2;
    else
      S_cap[i] = (counter_[i][index].V - counter_[i][index].C) / 2;
  }

  return *min_element(S_cap.begin(), S_cap.end());
}

template <int32_t key_len, typename T, typename hash_t>
void MVSketch<key_len, T, hash_t>::clear() {
  std::fill(counter_[0], counter_[0] + depth_ * width_, {0, 0, 0});
}

template <int32_t key_len, typename T, typename hash_t>
size_t MVSketch<key_len, T, hash_t>::size() const {
  return sizeof(MVSketch<key_len, T, hash_t>) + // Instance
         depth_ * sizeof(hash_t) +              // hash_fns
         sizeof(Bucket *) * depth_ +            // counter
         (2 * sizeof(T) + sizeof(FlowKey<key_len>)) * depth_ * width_;
}

template <int32_t key_len, typename T, typename hash_t>
typename MVSketch<key_len, T, hash_t>::Bounds
MVSketch<key_len, T, hash_t>::queryBounds(
    const FlowKey<key_len> &flow_key) const {
  std::vector<T> L(depth_);

  for (int i = 0; i < depth_; ++i) {
    int index = hash_fns_[i](flow_key) % width_;
    L[i] = counter_[i][index].K == flow_key ? counter_[i][index].C : 0;
  }

  return {*max_element(L.begin(), L.end()), query(flow_key)};
}

template <int32_t key_len, typename T, typename hash_t>
Data::Estimation<key_len, T> 
MVSketch<key_len, T, hash_t>::getHeavyHitter(double threshold) const {
  Data::Estimation<key_len, T> heavy_hitters;
  std::set<FlowKey<key_len>> heavy_set;

  for (int i = 0; i < depth_; ++i)
    for (int j = 0; j < width_; ++j) {
      if (counter_[i][j].V < threshold)
        continue;

      const FlowKey<key_len> &flow_key = counter_[i][j].K;
      if (heavy_set.count(flow_key))
        continue;

      T S_cap = query(flow_key);
      if (S_cap >= threshold){
        heavy_set.insert(flow_key);
        heavy_hitters[flow_key] = S_cap;
      }
    }

  return heavy_hitters;
}

} // namespace OmniSketch::Sketch
