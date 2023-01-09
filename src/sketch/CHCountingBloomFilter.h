/**
 * @file CHCountingBloomFilter.h
 * @author XierLabber (you@domain.com)
 * @brief CHCounting Bloom Filter
 *
 * @copyright Copyright (c) 2022
 *
 */
#pragma once

#include <common/hash.h>
#include <common/hierarchy.h>

// #define DEBUG

namespace OmniSketch::Sketch {
/**
 * @brief CHCounting Bloom Filter
 *
 * @tparam key_len  length of flowkey
 * @tparam hash_t   hashing class
 */
template <int32_t key_len, int32_t no_layer, typename hash_t>
class CHCountingBloomFilter : public SketchBase<key_len> {
  // for convenience
  using T = int64_t;
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
#ifdef DEBUG
  CH_ *counter_;
#endif

  CHCountingBloomFilter(const CHCountingBloomFilter &) = delete;
  CHCountingBloomFilter(CHCountingBloomFilter &&) = delete;
  CHCountingBloomFilter &operator=(CHCountingBloomFilter) = delete;

public:
  /**
   * @brief Construct by specifying #counters, #hash and length of counters
   *
   * @param num_cnt    #counter
   * @param num_hash    #hash
   */
  CHCountingBloomFilter(int32_t num_cnt, int32_t num_hash,
             double cnt_no_ratio,
             const std::vector<size_t> &width_cnt,
             const std::vector<size_t> &no_hash, 
             const int32_t cm_r,
             const int32_t cm_w);
  /**
   * @brief Destructor
   *
   */
  ~CHCountingBloomFilter();

  /**
   * @brief Insert a flowkey into the bloom filter
   * @details An overriding method
   *
   */
  void insert(const FlowKey<key_len> &flowkey) override;
  /**
   * @brief Look up a flowkey to see whether it exists
   * @details An overriding method
   *
   */
  bool lookup(const FlowKey<key_len> &flowkey) const override;
  /**
   * @brief Remove a flowkey from the bloom filter
   * @details A non-overriding method
   *
   */
  void remove(const FlowKey<key_len> &flowkey);
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

template <int32_t key_len, int32_t no_layer, typename hash_t>
CHCountingBloomFilter<key_len, no_layer, hash_t>::CHCountingBloomFilter(int32_t num_cnt,
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
#ifdef DEBUG
  counter_ = new CH_({static_cast<size_t>(ncnt)},
                   {static_cast<size_t>(cnt_length)}, {});
#endif
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

template <int32_t key_len, int32_t no_layer, typename hash_t>
CHCountingBloomFilter<key_len, no_layer, hash_t>::~CHCountingBloomFilter() {
  delete[] hash_fns;
  delete[] counter;
}

template <int32_t key_len, int32_t no_layer, typename hash_t>
void CHCountingBloomFilter<key_len, no_layer, hash_t>::insert(
    const FlowKey<key_len> &flowkey) {
  // if there is a 0
  int32_t i = 0;
  while (i < nhash) {
    int32_t idx = hash_fns[i](flowkey) % ncnt;
    if (counter->getEstCnt(idx) == 0){
    #ifdef DEBUG
      if(counter_->getCnt(idx) != 0){
        printf("INSERT JUDGE WRONG! COUNTER_: %ld, COUNTER: %ld\n", counter_->getCnt(idx), counter->getEstCnt(idx));
      }
    #endif
      break;
    }
    i++;
  }
  // increment the buckets
  if (i < nhash) {
    for (int32_t j = 0; j < nhash; ++j) {
    #ifdef DEBUG
      counter_->updateCnt(hash_fns[j](flowkey) % ncnt, 1);
    #endif
      counter->updateCnt(hash_fns[j](flowkey) % ncnt, 1);
    }
  }
}

template <int32_t key_len, int32_t no_layer, typename hash_t>
bool CHCountingBloomFilter<key_len, no_layer, hash_t>::lookup(
    const FlowKey<key_len> &flowkey) const {
  // if every counter is non-zero, return true
  for (int32_t i = 0; i < nhash; ++i) {
    int32_t idx = hash_fns[i](flowkey) % ncnt;
    #ifdef DEBUG
    if (counter_->getCnt(idx) == 0) {
      if(counter->getCnt(idx) != 0){
        printf("LOOK UP WRONG!COUNTER_: %ld, COUNTER: %ld, COUNT_ ORI: %ld, COUNT ORI: %ld\n", 
                counter_->getCnt(idx), counter->getCnt(idx), counter_->getOriginalCnt(idx), counter->getOriginalCnt(idx));
      }
      return false;
    }
    #else
    if (counter->getCnt(idx) == 0) {
      return false;
    }
    #endif
  }
  return true;
}

template <int32_t key_len, int32_t no_layer, typename hash_t>
void CHCountingBloomFilter<key_len, no_layer, hash_t>::remove(
    const FlowKey<key_len> &flowkey) {
  // if there is a 0
  int32_t i = 0;
  while (i < nhash) {
    int32_t idx = hash_fns[i](flowkey) % ncnt;
    #ifdef DEBUG
    if (counter_->getCnt(idx) == 0)
      break;
    #else
    if (counter->getCnt(idx) == 0)
      break;
    #endif
    i++;
  }
  // decrement the buckets
  if (i == nhash) {
    for (int32_t j = 0; j < nhash; ++j) {
    #ifdef DEBUG
      counter_->updateCnt(hash_fns[j](flowkey) % ncnt, -1);
    #endif
      counter->updateCnt(hash_fns[j](flowkey) % ncnt, -1);
    }
  }
}

template <int32_t key_len, int32_t no_layer, typename hash_t>
size_t CHCountingBloomFilter<key_len, no_layer, hash_t>::size() const {
  counter->print_rate("");
  return sizeof(*this)            // instance
         + sizeof(hash_t) * nhash // hash functions
         + counter->size();       // counter size
}

template <int32_t key_len, int32_t no_layer, typename hash_t>
void CHCountingBloomFilter<key_len, no_layer, hash_t>::clear() {
  counter->clear();
}

} // namespace OmniSketch::Sketch