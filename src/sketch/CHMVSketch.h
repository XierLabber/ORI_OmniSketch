/**
 * @file CHMVSketch.h
 * @author XierLabber
 * @brief Implementation of Count Min Sketch
 *
 * @copyright Copyright (c) 2022
 *
 */

#define VAL_WILL_ALWAYS_BE_ONE
#define RECORD_GUESS_WRONG_TIME

#include <algorithm>
#include <map>
#include <vector>

#include <common/hash.h>
#include <common/hierarchy.h>
#include <common/sketch.h>
#include <common/utils.h>

namespace OmniSketch::Sketch {
template <int32_t key_len, int32_t no_layer, typename T,
          typename hash_t = Hash::AwareHash>
class CHMVSketch : public SketchBase<key_len, T> {
  int32_t depth_;
  int32_t width_;

  hash_t *hash_fns_;

  struct Bounds {
    T lower;
    T upper;
  };

  FlowKey<key_len>** CounterK;

  std::vector<size_t> V_no_cnt;
  std::vector<size_t> V_width_cnt;
  std::vector<size_t> V_no_hash;

  CounterHierarchy<no_layer, T, hash_t>* CounterV;

  std::vector<size_t> C_no_cnt;
  std::vector<size_t> C_width_cnt;
  std::vector<size_t> C_no_hash;

  int32_t guess_negative_weight;
  int32_t guess_positive_weight;

  CounterHierarchy<no_layer, T, hash_t>* CounterC;
#ifdef RECORD_GUESS_WRONG_TIME
  int guess_wrong_time = 0;
#endif

public:

  size_t getIdx(size_t i, size_t j) const;

  CHMVSketch(
    int32_t depth, int32_t width, double V_cnt_no_ratio,
    const std::vector<size_t> &V_width_cnt, const std::vector<size_t> &V_no_hash, 
    double C_cnt_no_ratio,
    const std::vector<size_t> &C_width_cnt, const std::vector<size_t> &C_no_hash, 
    int32_t guess_negative_weight = 3, int32_t guess_positive_weight = 1);
  
  CHMVSketch(CHMVSketch &&) = delete;
  ~CHMVSketch();
  CHMVSketch &operator=(const CHMVSketch &) = delete;
  CHMVSketch &operator=(CHMVSketch &&) = delete;

  void update(const FlowKey<key_len> &flow_key, T val);
  T query(const FlowKey<key_len> &flow_key) const;

  void clear();
  size_t size() const;

  Data::Estimation<key_len, T>  getHeavyHitter(double threshold) const;
};

} // namespace OmniSketch::Sketch

//-----------------------------------------------------------------------------
//
///                        Implementation of templated methods
//
//-----------------------------------------------------------------------------

namespace OmniSketch::Sketch {

template <int32_t key_len, int32_t no_layer, typename T, typename hash_t>
size_t CHMVSketch<key_len, no_layer, T, hash_t>::getIdx(size_t i, size_t j) const{
  return i * width_ + j;
}

template <int32_t key_len, int32_t no_layer, typename T, typename hash_t>
CHMVSketch<key_len, no_layer, T, hash_t>::CHMVSketch(
    int32_t depth, int32_t width, double V_cnt_no_ratio,
    const std::vector<size_t> &V_width_cnt, const std::vector<size_t> &V_no_hash, 
    double C_cnt_no_ratio,
    const std::vector<size_t> &C_width_cnt, const std::vector<size_t> &C_no_hash,
    int32_t guess_negative_weight, int32_t guess_positive_weight)
    : depth_(depth), width_(Util::NextPrime(width)), V_width_cnt(V_width_cnt),
    V_no_hash(V_no_hash), C_width_cnt(C_width_cnt), C_no_hash(C_no_hash),
    guess_negative_weight(guess_negative_weight), guess_positive_weight(guess_positive_weight) {
  // check ratio
  if (V_cnt_no_ratio <= 0.0 || V_cnt_no_ratio >= 1.0) {
    throw std::out_of_range("Out of Range: Ratio of #counters of adjacent "
                            "layers in CH should be in (0, 1), but got " +
                            std::to_string(V_cnt_no_ratio) + " instead.");
  }

  if (C_cnt_no_ratio <= 0.0 || C_cnt_no_ratio >= 1.0) {
    throw std::out_of_range("Out of Range: Ratio of #counters of adjacent "
                            "layers in CH should be in (0, 1), but got " +
                            std::to_string(C_cnt_no_ratio) + " instead.");
  }
  
  hash_fns_ = new hash_t[depth_];

  CounterK = new FlowKey<key_len> *[depth_];
  CounterK[0] = new FlowKey<key_len>[depth_ * width_]();
  for(int i = 1; i < depth_; i++)
    CounterK[i] = CounterK[i - 1] + width_;

  // prepare no_cnt
  V_no_cnt.push_back(static_cast<size_t>(depth_) * width_);
  for (int32_t i = 1; i < no_layer; ++i) {
    size_t last_layer = V_no_cnt.back();
    V_no_cnt.push_back(Util::NextPrime(std::ceil(last_layer * V_cnt_no_ratio)));
  }
  C_no_cnt.push_back(static_cast<size_t>(depth_) * width_);
  for (int32_t i = 1; i < no_layer; ++i) {
    size_t last_layer = C_no_cnt.back();
    C_no_cnt.push_back(Util::NextPrime(std::ceil(last_layer * C_cnt_no_ratio)));
  }
  
  // CH
  CounterV = new CounterHierarchy<no_layer, T, hash_t>(V_no_cnt, this->V_width_cnt,
                                                      this->V_no_hash);
  CounterC = new CounterHierarchy<no_layer, T, hash_t>(C_no_cnt, this->C_width_cnt,
                                                      this->C_no_hash, true);
#ifdef RECORD_GUESS_WRONG_TIME
  guess_wrong_time = 0;
#endif
}

template <int32_t key_len, int32_t no_layer, typename T, typename hash_t>
CHMVSketch<key_len, no_layer, T, hash_t>::~CHMVSketch() {
  delete[] hash_fns_;

  delete[] CounterK[0];
  delete[] CounterK;

  delete[] CounterV;
  delete[] CounterC;
}

template <int32_t key_len, int32_t no_layer, typename T, typename hash_t>
void CHMVSketch<key_len, no_layer, T, hash_t>::update(const FlowKey<key_len> &flow_key,
                                          T val) {
  if(val != 1){
    printf("NOT 1! VAL = %d\n",val);
  }
  for (int i = 0; i < depth_; ++i) {
    int index = hash_fns_[i](flow_key) % width_;
    CounterV->updateCnt(getIdx(i, index), val);
    if(CounterK[i][index] == flow_key)
      CounterC->updateCnt(getIdx(i, index), val);
    else{
      CounterC->updateCnt(getIdx(i, index), -val);
      #ifndef VAL_WILL_ALWAYS_BE_ONE
      size_t chIdx = getIdx(i, index);
      bool my_guess = (CounterC->guess_sign(guess_negative_weight, guess_positive_weight, chIdx) < 0);
      if(my_guess){
        CounterK[i][index] = flow_key;
        CounterC->updateCnt(chIdx, 2 * val);
      }
      #else
      size_t chIdx = getIdx(i, index);
      bool my_guess = (CounterC->get_current_cnt(chIdx) + val == (1 << (C_width_cnt[0] - 1)) &&
                       CounterC->guess_sign(guess_negative_weight, guess_positive_weight, chIdx) < 0);
      if(my_guess){
        CounterK[i][index] = flow_key;
        CounterC->updateCnt(chIdx, 2 * val);
      }
      #endif

      #ifdef RECORD_GUESS_WRONG_TIME
      if((CounterC->getOriginalCnt(chIdx) >= 2 * val && my_guess) || 
         (CounterC->getOriginalCnt(chIdx) < 2 * val && (!my_guess))){
            guess_wrong_time++; 
         }
      #endif
    }
  }
}

template <int32_t key_len, int32_t no_layer, typename T, typename hash_t>
T CHMVSketch<key_len, no_layer, T, hash_t>::query(const FlowKey<key_len> &flow_key) const {
  std::vector<T> S_cap(depth_);

  for(int i = 0; i < depth_; i++){
    int index = hash_fns_[i](flow_key) % width_;
    size_t chIdx = getIdx(i, index);
    if(CounterK[i][index] == flow_key)
      S_cap[i] = (CounterV->getCnt(chIdx) + CounterC->getCnt(chIdx)) / 2;
    else
      S_cap[i] = (CounterV->getCnt(chIdx) - CounterC->getCnt(chIdx)) / 2;
  }

  return *min_element(S_cap.begin(), S_cap.end());
}

template <int32_t key_len, int32_t no_layer, typename T, typename hash_t>
void CHMVSketch<key_len, no_layer, T, hash_t>::clear() {
  std::fill(CounterK[0], 0);
  CounterC->clear();
  CounterV->clear();
}

template <int32_t key_len, int32_t no_layer, typename T, typename hash_t>
size_t CHMVSketch<key_len, no_layer, T, hash_t>::size() const {
#ifdef RECORD_GUESS_WRONG_TIME
  printf("GUESS WRONG %d TIMES!\n",guess_wrong_time);
#endif
  CounterV->print_rate("CounterV");
  CounterC->print_rate("CounterC");
  return sizeof(CHMVSketch<key_len, no_layer, T, hash_t>) + // Instance
         depth_ * sizeof(hash_t) +                          // hash_fns
         sizeof(FlowKey<key_len> *) * depth_ +              // counter
         sizeof(FlowKey<key_len>) * depth_ * width_ + 
         CounterC->size() + CounterV->size();
}

template <int32_t key_len, int32_t no_layer, typename T, typename hash_t>
Data::Estimation<key_len, T> 
CHMVSketch<key_len, no_layer, T, hash_t>::getHeavyHitter(double threshold) const {
  Data::Estimation<key_len, T> heavy_hitters;
  std::set<FlowKey<key_len>> heavy_set;

  for (int i = 0; i < depth_; ++i)
    for (int j = 0; j < width_; ++j) {
      size_t chIdx = getIdx(i, j);
      // if (counter_[i][j].V < threshold)
      if (CounterV->getCnt(chIdx) < threshold)
        continue;

      const FlowKey<key_len> &flow_key = CounterK[i][j];
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
