/**
 * @file HeavyKeeper.h
 * @author XierLabber (you@domain.com)
 * @brief Heavy Keeper
 *
 * @copyright Copyright (c) 2022
 *
 */
#pragma once

#include <common/hash.h>
#include <common/sketch.h>
#include <common/StreamSummary.h>

// #define USE_MAP
// #define DEBUG

namespace OmniSketch::Sketch {
/**
 * @brief Heavy Keeper
 *
 * @tparam key_len  length of flowkey
 * @tparam T        type of the counter
 * @tparam hash_t   hashing class
 */
template <int32_t key_len, typename T, typename hash_t = Hash::AwareHash>
class HeavyKeeper : public SketchBase<key_len, T> {
private:

  static const int32_t MINBUCKET = -1;
  static const int32_t MAXBUCKET = 0xfffffff;

  static const int32_t TOMB = 0xfffffff;

  int32_t depth_;
  int32_t width_;

  double b_;

  hash_t *sketch_hash_fun_;
  hash_t fingerprint_hash_fun_;

  class counter_t {
  public:
    T C;
    int16_t FP;
  };

  /*
  class hash_table_elem;
  class bucket_list_elem;
  class StreamSummary;
  */

  T n_min_;

  counter_t **counter_;


  double hash_table_alpha;


#ifndef USE_MAP
  StreamSummary<key_len, T, hash_t> StreamSummary_;
#else
  std::map<FlowKey<key_len>, int32_t> StreamSummary_;
#endif
  int32_t num_threshold_;

  HeavyKeeper(const HeavyKeeper &) = delete;
  HeavyKeeper(HeavyKeeper &&) = delete;
  HeavyKeeper &operator=(HeavyKeeper) = delete;

public:
  /**
   * @brief Construct by specifying depth, width, threshold size
   *        and the base used to calculate the probability of reduction
   *
   */
  HeavyKeeper(int32_t depth, int32_t width, int32_t num_threshold, double b, double hash_table_alpha);
  /**
   * @brief Release the pointer
   *
   */
  ~HeavyKeeper();
  /**
   * @brief Update a flowkey with certain value
   *        val is a useless parameter here, as we only count
   *        the number of times flowkey occurs
   *
   */
  void update(const FlowKey<key_len> &flowkey, T val);
  /**
   * @brief Get the K_th largest elem of the Array who has
   *        length elems.
   */
  T findKth(T *Array, int32_t k, int32_t length);
  /**
   * @brief Get the size of the sketch
   *
   */
  size_t size() const;
  /**
   * @brief Get the k flows with the largest occurance time
   *
   */
  Data::Estimation<key_len, T> getTopK(int32_t k) const;
  /**
   * @brief Get Heavy Hitter
   *
   */
  Data::Estimation<key_len, T> getHeavyHitter(double val_threshold) const;
};

} // namespace OmniSketch::Sketch

//-----------------------------------------------------------------------------
//
///                        Implementation of template methods
//
//-----------------------------------------------------------------------------

namespace OmniSketch::Sketch {

template <int32_t key_len, typename T, typename hash_t>
HeavyKeeper<key_len, T, hash_t>::HeavyKeeper(
    int32_t depth, int32_t width, int32_t num_threshold, double b, double hash_table_alpha)
    : depth_(depth), width_(Util::NextPrime(width)), hash_table_alpha(hash_table_alpha), 
      num_threshold_(num_threshold), b_(b), n_min_(0) {

  sketch_hash_fun_ = new hash_t[depth_];

  // Allocate continuous memory
  counter_ = new counter_t *[depth_];
  counter_[0] = new counter_t[depth_ * this->width_]();
  for (int32_t i = 1; i < depth_; ++i) {
    counter_[i] = counter_[i - 1] + this->width_;
  }

#ifndef USE_MAP
  StreamSummary_.init(num_threshold, hash_table_alpha);
#endif

}

template <int32_t key_len, typename T, typename hash_t>
HeavyKeeper<key_len, T, hash_t>::~HeavyKeeper() {
  delete[] sketch_hash_fun_;

  delete[] counter_[0];
  delete[] counter_;
  StreamSummary_.destroy();
}

template <int32_t key_len, typename T, typename hash_t>
size_t HeavyKeeper<key_len, T, hash_t>::size() 
  const{
  return depth_ * width_ * (sizeof(T) + sizeof(uint16_t))
         + depth_ * sizeof(hash_t)
         + sizeof(HeavyKeeper<key_len, hash_t, T>)
         + StreamSummary_.memory_size();
}

template <int32_t key_len, typename T, typename hash_t>
void HeavyKeeper<key_len, T, hash_t>::update(
    const FlowKey<key_len> &flowkey, T val) {
  auto iter = StreamSummary_.find(flowkey);
  #ifndef USE_MAP
  bool flag = (iter != NULL);
  #else
  bool flag = (iter != StreamSummary_.end());
  #endif

  int32_t EstimatedCount = MINBUCKET;
  int16_t FlowFP = fingerprint_hash_fun_(flowkey);

  bool done = false;

  for (int i = 0; i < depth_; i++) {
      int32_t index = sketch_hash_fun_[i](flowkey) % width_;
      int32_t counterC = counter_[i][index].C;
      if (counterC > 0 && (flag || counterC < n_min_) && counter_[i][index].FP == FlowFP) {
      counterC += val;
      EstimatedCount = (counterC < EstimatedCount) ? EstimatedCount : counterC;
      counter_[i][index].C = counterC;
      done = true;
      }
  }

  if (!done) {
    for (int i = 0; i < depth_; i++) {
      int32_t index = sketch_hash_fun_[i](flowkey) % width_;
      if (counter_[i][index].C == 0) {
        counter_[i][index].C = val;
        counter_[i][index].FP = FlowFP;
        EstimatedCount = 1;
        done = true;
        break;
      }
    }
  }

  if (!done) {
    int32_t minC = counter_[0][sketch_hash_fun_[0](flowkey) % width_].C;
    int32_t minCounterID = 0;
    for (int i = 1; i < depth_; i++) {
      int32_t index = sketch_hash_fun_[i](flowkey) % width_;
      int32_t counterC = counter_[i][index].C;
      if (counterC < minC) {
        minC = counterC;
        minCounterID = i;
      }
    }
    int32_t minIndex = sketch_hash_fun_[minCounterID](flowkey) % width_;
    if (((double)rand() / RAND_MAX) < (pow(b_, -minC))) {
      counter_[minCounterID][minIndex].C -= val;
      if (counter_[minCounterID][minIndex].C <= 0) {
        counter_[minCounterID][minIndex].C = 1;
        counter_[minCounterID][minIndex].FP = FlowFP;
        EstimatedCount = 1;
      }
    }
  }

  if (EstimatedCount > 0) {
    if (flag) {
      if(EstimatedCount > iter->second){
        StreamSummary_.increment(iter, EstimatedCount - iter->second);
      }
      if (n_min_ == iter->second) {
        #ifndef USE_MAP
        auto min_iter = StreamSummary_.get_least_elem();
        n_min_ = (min_iter == NULL)? 0 : min_iter->second;
        #else
        auto min_iter =
          std::min_element(StreamSummary_.begin(), StreamSummary_.end(),
                          [](const std::pair<FlowKey<key_len>, int32_t> &left,
                              const std::pair<FlowKey<key_len>, int32_t> &right) {
                          return left.second < right.second;
                          });
          n_min_=min_iter->second;
        #endif
      }
    } else if (StreamSummary_.size() < num_threshold_) {
      StreamSummary_.emplace(flowkey, EstimatedCount);
      n_min_ = std::min(EstimatedCount, n_min_);
    } else if (EstimatedCount > n_min_) {
      #ifndef USE_MAP
      auto min_iter = StreamSummary_.get_least_elem();
      StreamSummary_.emplace(flowkey, EstimatedCount);
      StreamSummary_.erase(min_iter);
      min_iter = StreamSummary_.get_least_elem();
      n_min_ = (min_iter == NULL)? 0 : min_iter->second;
      #else
      auto min_iter =
      std::min_element(StreamSummary_.begin(), StreamSummary_.end(),
                      [](const std::pair<FlowKey<key_len>, int32_t> &left,
                          const std::pair<FlowKey<key_len>, int32_t> &right) {
                      return left.second < right.second;
                      });
      StreamSummary_.emplace(flowkey, EstimatedCount);
      StreamSummary_.erase(min_iter);
      min_iter =
          std::min_element(StreamSummary_.begin(), StreamSummary_.end(),
                          [](const std::pair<FlowKey<key_len>, int32_t> &left,
                              const std::pair<FlowKey<key_len>, int32_t> &right) {
                          return left.second < right.second;
                          });
      n_min_=min_iter->second;
      #endif
    }
  }
}

template <int32_t key_len, typename T, typename hash_t>
T HeavyKeeper<key_len, T, hash_t>::findKth(
    T *Array, int32_t k, int32_t length) {
  int32_t left = 0, right = length - 1;
  while (left < right) {
    while (Array[right] < Array[left])
      right--;
    int32_t tmp = Array[right];
    Array[right] = Array[left];
    Array[left] = tmp;
    if (right <= left)
      break;
    while (Array[left] > Array[right])
      left++;
    tmp = Array[right];
    Array[right] = Array[left];
    Array[left] = tmp;
  }
  if (left == k - 1) {
    return Array[left];
  } else if (left < k - 1) {
    return findKth(Array + left + 1, k - left - 1, length - left - 1);
  } else {
    return findKth(Array, k, left);
  }
}

template <int32_t key_len, typename T, typename hash_t>
Data::Estimation<key_len, T>
HeavyKeeper<key_len, T, hash_t>::getTopK(int32_t k) const{
  T ValueArray[num_threshold_];
  int32_t ArrayLength = 0;
  #ifndef USE_MAP
  for(int i=0; i < StreamSummary_.get_hashtable_length(); i++)
  {
    T val = StreamSummary_.get_hashtable_val(i);
    if(val != MINBUCKET && val != TOMB)
    {
        ValueArray[ArrayLength++] = val;
    }
  }
  #else
  auto iter = StreamSummary_.begin();
  while (iter != StreamSummary_.end()) {
    ValueArray[ArrayLength++] = iter->second;
    iter++;
  }
  #endif
  int32_t KthValue = findKth(ValueArray, k, ArrayLength);
  Data::Estimation<key_len, T> TopK;
  #ifndef USE_MAP
  for(int i=0; i < StreamSummary_.get_hashtable_length(); i++)
  {
    T val = StreamSummary_.get_hashtable_val(i);
    if(val != MINBUCKET && val != TOMB && val >= KthValue)
    {
        TopK[StreamSummary_.get_hashtable_key(i)] = StreamSummary_.get_hashtable_val(i);
    }
  }
  #else
  iter = StreamSummary_.begin();
  while (iter != StreamSummary_.end()) {
    if (iter->second >= KthValue)
      TopK[iter->first] = iter->second;
    iter++;
  }
  #endif
  return TopK;
}

template <int32_t key_len, typename T, typename hash_t>
Data::Estimation<key_len, T>
HeavyKeeper<key_len, T, hash_t>::getHeavyHitter(
    double val_threshold) const {
  #ifndef USE_MAP
  return StreamSummary_.getHeavyHitter(val_threshold);
  #else
  Data::Estimation<key_len, T> heavy_Hitter;
  auto iter = StreamSummary_.begin();
  int found = 0;
  iter = StreamSummary_.begin();
  while (iter != StreamSummary_.end()) {
    if (iter->second >= val_threshold)
    {
      heavy_Hitter[iter->first] = iter->second;
      found++;
    }
    iter++;
  }
  return heavy_Hitter;
  #endif
}

} // namespace OmniSketch