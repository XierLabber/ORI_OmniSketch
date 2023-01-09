/**
 * @file CounterBraids.h
 * @author dromniscience (you@domain.com)
 * @brief Implementation of Counter Braids
 *
 * @copyright Copyright (c) 2022
 *
 */
#pragma once

#include <boost/dynamic_bitset.hpp>
#include <common/hash.h>
#include <common/sketch.h>
#include <iostream>
#include <unordered_map>

namespace OmniSketch::Sketch {
/**
 * @brief Counter Braids
 *
 * @tparam key_len  length of flowkey
 * @tparam T        type of the counter
 * @tparam hash_t   hashing class
 */
template <int32_t key_len, int32_t no_layer, typename T,
          typename hash_t = Hash::AwareHash>
class CounterBraids : public SketchBase<key_len, T> {
private:
  using CarryOver = std::map<std::size_t, T>;
  /**
   * @brief Number of counters on each layer, from low to high.
   *
   */
  std::vector<size_t> no_cnt;
  /**
   * @brief Width of counters on each layer, from low to high.
   *
   */
  const std::vector<size_t> width_cnt;
  /**
   * @brief Number of hash function used on each layer, from low to
   * high (except for the last layer)
   *
   */
  const std::vector<size_t> no_hash;
  /**
   * @brief Counters
   *
   */
  std::vector<Util::DynamicIntX<T>> *cnt_array;
  /**
   * @brief Status bits
   *
   */
  boost::dynamic_bitset<uint8_t> *status_bits;
  /**
   * @brief array of hashing classes
   *
   */
  std::vector<hash_t> *hash_fns;
  /**
   * @brief For lazy update policy
   *
   */
  CarryOver lazy_update;
  /**
   * @brief All flowkeys
   * @details Typically this should be provided by other entities.
   * For convenience, CB collects keys for itself.
   * One may substitute this with GndTruth + BloomFilter.
   */
  Data::Estimation<key_len, T> key_set;
  /**
   * @brief get decoded counter
   */
  std::vector<T> decoded_cnt;

  CounterBraids(const CounterBraids &) = delete;
  CounterBraids(CounterBraids &&) = delete;

  /**
   * @brief Update a layer (aggregation)
   *
   * @details An overflow error would be thrown if there is an overflow at the
   * last layer. For the other layers, since a counter may first oveflow and
   * then be substracted to withdraw any carry over, the number of overflows of
   * counter whose status bit is set is not assumed to be 1 at least.
   *
   * @param layer   the current layer
   * @param updates updates to be propagated to the current layer
   * @return updates to be propagated to the next layer
   */
  [[nodiscard]] CarryOver updateLayer(const int32_t layer, CarryOver &&updates);
  /**
   * @brief Decode a layer
   *
   * @param layer   the higher layer
   * @param higher  decoded results of the higher layer
   * @return results of the current layer
   */
  [[nodiscard]] std::vector<T> decodeLayer(const int32_t layer,
                                           std::vector<T> &&higher) const;

public:
  /**
   * @brief Constructor
   * @param no_cnt      number of counters on each layer, from low to high
   * @param width_cnt   width of counters on each layer, from low to high
   * @param no_hash     number of hash functions used on each layer, from low
   * to high
   *
   * @details The meaning of the three parameters stipulates the following
   * requirements:
   * - size of `no_cnt` should equal `no_layer`.
   * - size of `width_cnt` should equal `no_layer`.
   * - size of `no_hash` should equal `no_layer`.
   *
   * If any of these is violated, an exception would be thrown. Other
   * conditions that trigger an exception:
   * - Items in these vectors contains a 0.
   * - `no_layer <= 0`
   * - Sum of `width_cnt` exceeds `sizeof(T) * 8`. This constraint is imposed to
   * guarantee proper shifting of counters when decoding.
   *
   */
  CounterBraids(const std::vector<size_t> &no_cnt,
                const std::vector<size_t> &width_cnt,
                const std::vector<size_t> &no_hash);
  /**
   * @brief Release the pointer
   *
   */
  ~CounterBraids();
  /**
   * @brief Update a flowkey with certain value
   *
   */
  void update(const FlowKey<key_len> &flowkey, T val) override;
  /**
   * @brief Query a flowkey
   *
   */
  Data::Estimation<key_len, T> decode() override;
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

template <int32_t key_len, int32_t no_layer, typename T, typename hash_t>
typename CounterBraids<key_len, no_layer, T, hash_t>::CarryOver
CounterBraids<key_len, no_layer, T, hash_t>::updateLayer(const int32_t layer,
                                                         CarryOver &&updates) {
  CarryOver ret; // aggregate all updates on the current layer
  for (const auto &kv : updates) {
    T overflow = cnt_array[layer][kv.first] + kv.second;
    if (overflow) {
      // mark status bits
      status_bits[layer][kv.first] = true;
      if (layer == no_layer - 1) { // last layer
        throw std::overflow_error(
            "Counter overflow at the last layer in CB, overflow by " +
            std::to_string(overflow) + ".");
      } else { // hash to upper-layer counters
        for (size_t i = 0; i < no_hash[layer + 1]; i++) {
          std::size_t index =
              hash_fns[layer + 1][i](kv.first) % no_cnt[layer + 1];
          ret[index] += overflow;
        }
      }
    }
  }
  return ret;
}

#define ITER 10
#define LBOUND 1
#define UBOUND std::numeric_limits<T>::max()
template <int32_t key_len, int32_t no_layer, typename T, typename hash_t>
std::vector<T> CounterBraids<key_len, no_layer, T, hash_t>::decodeLayer(
    const int32_t layer, std::vector<T> &&higher) const {

  using Graph = std::vector<std::unordered_map<int32_t, T>>;
  // sparse random bipartite graph
  int32_t left, right = no_cnt[layer];
  // if layer == 0, `left` equals #flows.
  left = layer ? no_cnt[layer - 1] : key_set.size();
  // from right to left
  Graph graph(right);
  // estimate
  std::vector<T> estimate(left);
  // Optimization
  std::vector<T> second_to_last(left);

  // build the graph
  if (layer) {
    for (int32_t i = 0; i < left; ++i) {
      if (status_bits[layer - 1][i]) {
        estimate[i] = LBOUND;
        for (int32_t j = 0; j < no_hash[layer]; ++j) {
          int32_t k = hash_fns[layer][j](i) % right;
          graph[k].emplace(i, 0);
        }
      }
    }
  } else {
    int32_t i = 0;
    for (const auto &kv : key_set) {
      estimate[i] = LBOUND;
      for (int32_t j = 0; j < no_hash[0]; ++j) {
        int32_t k = hash_fns[0][j](kv.get_left()) % right;
        graph[k].emplace(i, 0);
      }
      i++;
    }
  }

  // iteration
  for (int32_t iter = 0; iter < ITER; ++iter) {
    // forward message
    for (int32_t j = 0; j < right; ++j) {
      T acc = 0;
      for (auto &k : graph[j]) {
        acc += estimate[k.first];
      }
      acc = cnt_array[layer][j].getVal() - acc;
      for (auto &k : graph[j]) {
        k.second = std::max(acc + estimate[k.first], LBOUND);
      }
    }
    // backward message
    for (int32_t j = 0; j < left; ++j) {
      if (!estimate[j])
        continue;
      if (iter & 1) {
        estimate[j] = LBOUND;
      } else {
        estimate[j] = UBOUND;
      }
    }
    for (int32_t j = 0; j < right; ++j) {
      for (auto &k : graph[j]) {
        if (iter & 1) {
          estimate[k.first] = std::max(estimate[k.first], k.second);
        } else {
          estimate[k.first] = std::min(estimate[k.first], k.second);
        }
      }
    }
    // Optimization
    // second to last
    if (iter == ITER - 2) {
      std::copy(estimate.begin(), estimate.end(), second_to_last.begin());
    } else if (iter == ITER - 1) {
      for (int32_t i = 0; i < left; ++i) {
        estimate[i] = (estimate[i] + second_to_last[i]) >> 1;
      }
    }
  } // End iteration

  // shift and add
  if (layer) {
    for (int32_t i = 0; i < left; ++i) {
      estimate[i] = (estimate[i] << width_cnt[layer - 1]) +
                    cnt_array[layer - 1][i].getVal();
    }
  }
  return estimate;
}
#undef UBOUND
#undef LBOUND
#undef ITER

template <int32_t key_len, int32_t no_layer, typename T, typename hash_t>
CounterBraids<key_len, no_layer, T, hash_t>::CounterBraids(
    const std::vector<size_t> &no_cnt0, const std::vector<size_t> &width_cnt,
    const std::vector<size_t> &no_hash)
    : no_cnt(no_cnt0), width_cnt(width_cnt), no_hash(no_hash) {
  for(int i = 0; i < no_cnt.size(); i++)
  {
    no_cnt[i] = Util::NextPrime(no_cnt[i]);
  }
  // validity check
  if (no_layer < 1) {
    throw std::invalid_argument(
        "Invalid Template Argument: `no_layer` must > 1, got " +
        std::to_string(no_layer) + ".");
  }
  if (no_cnt.size() != no_layer) {
    throw std::invalid_argument(
        "Invalid Argument: `no_cnt` should be of size " +
        std::to_string(no_layer) + ", but got size " +
        std::to_string(no_cnt.size()) + ".");
  }
  if (width_cnt.size() != no_layer) {
    throw std::invalid_argument(
        "Invalid Argument: `width_cnt` should be of size " +
        std::to_string(no_layer) + ", but got size " +
        std::to_string(width_cnt.size()) + ".");
  }
  if (no_hash.size() != no_layer) {
    throw std::invalid_argument(
        "Invalid Argument: `no_hash` should be of size " +
        std::to_string(no_layer) + ", but got size " +
        std::to_string(no_hash.size()) + ".");
  }
  for (auto i : no_cnt) {
    if (i == 0) {
      throw std::invalid_argument(
          "Invalid Argument: There is a zero in `no_cnt`.");
    }
  }
  for (auto i : width_cnt) {
    if (i == 0) {
      throw std::invalid_argument(
          "Invalid Argument: There is a zero in `width_cnt`.");
    }
  }
  for (auto i : no_hash) {
    if (i == 0) {
      throw std::invalid_argument(
          "Invalid Argument: There is a zero in `no_hash`.");
    }
  }
  size_t length = 0;
  for (auto i : width_cnt) {
    size_t tmp = length + i;
    if (tmp < length || tmp > sizeof(T) * 8) {
      throw std::invalid_argument(
          "Invalid Argument: Aggregate length of `width_cnt` is too large.");
    }
    length = tmp;
  }

  // allocate in heap
  hash_fns = new std::vector<hash_t>[no_layer];
  for (int32_t i = 0; i < no_layer; ++i) {
    hash_fns[i] = std::vector<hash_t>(no_hash[i]);
  }
  cnt_array = new std::vector<Util::DynamicIntX<T>>[no_layer];
  for (int32_t i = 0; i < no_layer; ++i) {
    cnt_array[i] = std::vector<Util::DynamicIntX<T>>(no_cnt[i], width_cnt[i]);
  }
  status_bits = new boost::dynamic_bitset<uint8_t>[no_layer];
  for (int32_t i = 0; i < no_layer; ++i) {
    status_bits[i].resize(no_cnt[i], false);
  }
  // decoded counters, value initialized
  // decoded_cnt.resize(no_cnt[0]);
}

template <int32_t key_len, int32_t no_layer, typename T, typename hash_t>
CounterBraids<key_len, no_layer, T, hash_t>::~CounterBraids() {
  if (hash_fns)
    delete[] hash_fns;
  if (cnt_array)
    delete[] cnt_array;
  if (status_bits)
    delete[] status_bits;
}

template <int32_t key_len, int32_t no_layer, typename T, typename hash_t>
void CounterBraids<key_len, no_layer, T, hash_t>::update(
    const FlowKey<key_len> &flowkey, T val) {
  // update to layer 0
  for (int32_t i = 0; i < no_hash[0]; ++i) {
    lazy_update[hash_fns[0][i](flowkey) % no_cnt[0]] += val;
  }
  // record the flow key
  key_set.insert(flowkey);
}

template <int32_t key_len, int32_t no_layer, typename T, typename hash_t>
Data::Estimation<key_len, T>
CounterBraids<key_len, no_layer, T, hash_t>::decode() {
  // lazy update
  for (int32_t i = 0; i < no_layer; i++) {
    lazy_update = updateLayer(i, std::move(lazy_update)); // throw exception
  }
  lazy_update.clear();
  // decode
  decoded_cnt.resize(no_cnt[no_layer - 1]);
  for (size_t i = 0; i < no_cnt[no_layer - 1]; ++i) {
    decoded_cnt[i] = cnt_array[no_layer - 1][i].getVal();
  }
  for (int32_t i = no_layer - 1; i >= 0; i--) {
    decoded_cnt = decodeLayer(i, std::move(decoded_cnt));
  }
  int32_t index = 0;
  for (auto &kv : key_set) {
    key_set[kv.get_left()] = decoded_cnt[index];
    index++;
  }
  return key_set;
}

template <int32_t key_len, int32_t no_layer, typename T, typename hash_t>
size_t CounterBraids<key_len, no_layer, T, hash_t>::size() const {
  // counters + status bits
  size_t tot = 0; // first in bits
  for (int32_t i = 0; i < no_layer; ++i) {
    tot += no_cnt[i] * width_cnt[i];
    tot += no_cnt[i];
  }
  tot >>= 3; // to bytes
  // hash functions
  for (int32_t i = 0; i < no_layer; ++i) {
    tot += sizeof(hash_t) * no_hash[i];
  }
  return tot;
}

template <int32_t key_len, int32_t no_layer, typename T, typename hash_t>
void CounterBraids<key_len, no_layer, T, hash_t>::clear() {
  // reset counters
  for (int32_t i = 0; i < no_layer; ++i) {
    cnt_array[i] = std::vector<Util::DynamicIntX<T>>(no_cnt[i], width_cnt[i]);
  }
  // reset status bits
  for (int32_t i = 0; i < no_layer; ++i) {
    status_bits[i].reset();
  }
  // reset lazy update
  lazy_update.clear();
  // reset key set
  key_set = Data::GndTruth<key_len, T>();
}

} // namespace OmniSketch::Sketch
