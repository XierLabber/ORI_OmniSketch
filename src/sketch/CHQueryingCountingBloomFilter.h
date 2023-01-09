/**
 * @file CHQueryingCountingBloomFilter.h
 * @author XierLabber (you@domain.com)
 * @brief CHCounting Bloom Filter
 *
 * @copyright Copyright (c) 2022
 *
 */
#pragma once

#include <common/hash.h>
#include <common/hierarchy.h>

namespace OmniSketch::Sketch {
/**
 * @brief CHCounting Bloom Filter
 *
 * @tparam key_len  length of flowkey
 * @tparam hash_t   hashing class
 */
template <int32_t key_len, int32_t no_layer, typename T, typename hash_t>
class CHQueryingCountingBloomFilter : public SketchBase<key_len, T> {
  using CH = CounterHierarchy<no_layer, T, hash_t>;
  using CH_ = CounterHierarchy<1, T, hash_t>;

private:
  int32_t ncnt;
  int32_t nhash;
  hash_t *hash_fns;

  std::vector<size_t> no_cnt;
  std::vector<size_t> width_cnt;
  std::vector<size_t> no_hash;
  CH *counter;

  CHQueryingCountingBloomFilter(const CHQueryingCountingBloomFilter &) = delete;
  CHQueryingCountingBloomFilter(CHQueryingCountingBloomFilter &&) = delete;
  CHQueryingCountingBloomFilter &operator=(CHQueryingCountingBloomFilter) = delete;

public:
  /**
   * @brief Construct by specifying #counters, #hash and length of counters
   *
   * @param num_cnt    #counter
   * @param num_hash    #hash
   */
  CHQueryingCountingBloomFilter(int32_t num_cnt, int32_t num_hash,
             double cnt_no_ratio,
             const std::vector<size_t> &width_cnt,
             const std::vector<size_t> &no_hash, 
             const int32_t cm_r,
             const int32_t cm_w);
  /**
   * @brief Destructor
   *
   */
  ~CHQueryingCountingBloomFilter();
  
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
   * @brief Size of the sketch
   * @details An overriding method
   */
  size_t size() const override;
  /**
   * @brief Reset the Bloom Filter
   * @details A non-overriding method
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

template <int32_t key_len, int32_t no_layer, typename T, typename hash_t>
CHQueryingCountingBloomFilter<key_len, no_layer, T, hash_t>::CHQueryingCountingBloomFilter(int32_t num_cnt,
                                                          int32_t num_hash,
                                                          double cnt_no_ratio,
                                                          const std::vector<size_t> &width_cnt,
                                                          const std::vector<size_t> &no_hash, 
                                                          const int32_t cm_r,
                                                          const int32_t cm_w)
    : ncnt(Util::NextPrime(num_cnt)), nhash(num_hash),
      width_cnt(width_cnt), no_hash(no_hash)  {
  // hash functions
  hash_fns = new hash_t[num_hash];
  // counter array
  size_t cnt_length = 0;
  for(int i = 0; i < width_cnt.size(); i++){
    cnt_length += width_cnt[i];
  }
  // check ratio
  if (cnt_no_ratio <= 0.0 || cnt_no_ratio >= 1.0) {
    throw std::out_of_range("Out of Range: Ratio of #counters of adjacent "
                            "layers in CH should be in (0, 1), but got " +
                            std::to_string(cnt_no_ratio) + " instead.");
  }
  // prepare no_cnt
  no_cnt.push_back(static_cast<size_t>(this->ncnt));
  for (int32_t i = 1; i < no_layer; ++i) {
    size_t last_layer = no_cnt.back();
    no_cnt.push_back(Util::NextPrime(std::ceil(last_layer * cnt_no_ratio)));
  }
  // CH
  counter = new CounterHierarchy<no_layer, T, hash_t>(no_cnt, this->width_cnt,
                                                      this->no_hash, false, true, 
                                                      cm_r, cm_w);
}

template <int32_t key_len, int32_t no_layer, typename T, typename hash_t>
CHQueryingCountingBloomFilter<key_len, no_layer, T, hash_t>::~CHQueryingCountingBloomFilter() {
  delete[] hash_fns;
  delete[] counter;
}

template <int32_t key_len, int32_t no_layer, typename T, typename hash_t>
void CHQueryingCountingBloomFilter<key_len, no_layer, T, hash_t>::update(const FlowKey<key_len> 
    &flowkey, T val) {
    for (int32_t j = 0; j < nhash; ++j) {
      counter->updateCnt(hash_fns[j](flowkey) % ncnt, val);
    }
}

template <int32_t key_len, int32_t no_layer, typename T, typename hash_t>
T CHQueryingCountingBloomFilter<key_len, no_layer, T, hash_t>::query(const FlowKey<key_len> &flowkey) const{
  T minn = counter->getCnt(hash_fns[0](flowkey) % ncnt);
  // if every counter is non-zero, return true
  for (int32_t i = 1; i < nhash; ++i) {
    int32_t idx = hash_fns[i](flowkey) % ncnt;
    T cur = counter->getCnt(idx);
    minn = std::min(minn, cur);
  }
  return minn;
}

template <int32_t key_len, int32_t no_layer, typename T, typename hash_t>
size_t CHQueryingCountingBloomFilter<key_len, no_layer, T, hash_t>::size() const {
  counter->print_rate("");
  return sizeof(*this)            // instance
         + sizeof(hash_t) * nhash // hash functions
         + counter->size();       // counter size
}

template <int32_t key_len, int32_t no_layer, typename T, typename hash_t>
void CHQueryingCountingBloomFilter<key_len, no_layer, T, hash_t>::clear() {
  counter->clear();
}

} // namespace OmniSketch::Sketch