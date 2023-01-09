/**
 * @file QueryingCountingBloomFilter.h
 * @author XierLabber (you@domain.com)
 * @brief Querying Counting Bloom Filter
 *
 * @copyright Copyright (c) 2022
 *
 */
#pragma once

#include <common/hash.h>
#include <common/hierarchy.h>

namespace OmniSketch::Sketch {
/**
 * @brief Querying Counting Bloom Filter
 *
 * @tparam key_len  length of flowkey
 * @tparam hash_t   hashing class
 */
template <int32_t key_len, typename T, typename hash_t>
class QueryingCountingBloomFilter : public SketchBase<key_len, T> {
  using CH = CounterHierarchy<1, T, hash_t>;

private:
  int32_t ncnt;
  int32_t nhash;
  hash_t *hash_fns;
  CH *counter;

  QueryingCountingBloomFilter(const QueryingCountingBloomFilter &) = delete;
  QueryingCountingBloomFilter(QueryingCountingBloomFilter &&) = delete;
  QueryingCountingBloomFilter &operator=(QueryingCountingBloomFilter) = delete;

public:
  /**
   * @brief Construct by specifying #counters, #hash and length of counters
   *
   * @param num_cnt    #counter
   * @param num_hash    #hash
   * @param cnt_length  length of each counter
   */
  QueryingCountingBloomFilter(int32_t num_cnt, int32_t num_hash, int32_t cnt_length);
  /**
   * @brief Destructor
   *
   */
  ~QueryingCountingBloomFilter();
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

template <int32_t key_len, typename T, typename hash_t>
QueryingCountingBloomFilter<key_len, T, hash_t>::QueryingCountingBloomFilter(int32_t num_cnt,
                                                          int32_t num_hash,
                                                          int32_t cnt_length)
    : ncnt(Util::NextPrime(num_cnt)), nhash(num_hash) {
  // hash functions
  hash_fns = new hash_t[num_hash];
  // counter array
  counter = new CH({static_cast<size_t>(ncnt)},
                   {static_cast<size_t>(cnt_length)}, {});
}

template <int32_t key_len, typename T, typename hash_t>
QueryingCountingBloomFilter<key_len, T, hash_t>::~QueryingCountingBloomFilter() {
  delete[] hash_fns;
  delete counter;
}

template <int32_t key_len, typename T, typename hash_t>
void QueryingCountingBloomFilter<key_len, T, hash_t>::update(const FlowKey<key_len> 
    &flowkey, T val) {
    for (int32_t j = 0; j < nhash; ++j) {
      counter->updateCnt(hash_fns[j](flowkey) % ncnt, val);
    }

}

template <int32_t key_len, typename T, typename hash_t>
T QueryingCountingBloomFilter<key_len, T, hash_t>::query(const FlowKey<key_len> &flowkey) const{
  T minn = counter->getCnt(hash_fns[0](flowkey) % ncnt);
  // if every counter is non-zero, return true
  for (int32_t i = 1; i < nhash; ++i) {
    int32_t idx = hash_fns[i](flowkey) % ncnt;
    T cur = counter->getCnt(idx);
    minn = std::min(minn, cur);
  }
  return minn;
}

template <int32_t key_len, typename T, typename hash_t>
size_t QueryingCountingBloomFilter<key_len, T, hash_t>::size() const {
  counter->print_rate("");
  return sizeof(*this)            // instance
         + sizeof(hash_t) * nhash // hash functions
         + counter->size()        // counter size
         - ncnt / 8;
}

template <int32_t key_len, typename T, typename hash_t>
void QueryingCountingBloomFilter<key_len, T, hash_t>::clear() {
  counter->clear();
}

} // namespace OmniSketch::Sketch