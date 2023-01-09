/**
 * @file THD_CHHHUnivMon.h
 * @author XierLabber (you@domain.com)
 * @brief Implementation of THD_CHHHUnivMon
 *
 * @copyright Copyright (c) 2022
 *
 */
#pragma once
#include <iostream>

#include <common/hash.h>
#include <common/hierarchy.h>
#include <common/StreamSummary.h>
#include <common/sketch.h>
#include <sketch/CountSketch.h>

namespace OmniSketch::Sketch {
/**
 * @brief THD_CHHHUnivMon
 *
 * @tparam key_len  length of flowkey
 * @tparam T        type of the counter
 * @tparam hash_t   hashing class
 */
template <int32_t key_len, int32_t no_layer, typename T, typename hash_t>
class THD_CHHHUnivMon : public SketchBase<key_len, T> {
private:
  int32_t depth;
  int32_t total_length;
  int32_t* width;
  int32_t* width_idx;
  int32_t logn;
  hash_t *global_hash_fns;
  hash_t **CS_hash_fns;
  hash_t **CS_update_hash_fns;
  
  std::vector<size_t> no_cnt;
  std::vector<size_t> width_cnt;
  std::vector<size_t> no_hash;

  CounterHierarchy<no_layer, T, hash_t> *ch;

  Data::Estimation<key_len> *flows;

  THD_CHHHUnivMon(const THD_CHHHUnivMon &) = delete;
  THD_CHHHUnivMon(THD_CHHHUnivMon &&) = delete;

  T* sum;
  StreamSummary<key_len, T, hash_t>* HHHeaps;

public:
  /**
   * @brief Construct by specifying depth, width and $\log n$, where $n$ is the
   * number of flows to insert.
   *
   */
  THD_CHHHUnivMon(int32_t depth_, int32_t width_, int32_t log_n, 
            int32_t heap_size, double SSalpha,
            double cnt_no_ratio,
            const std::vector<size_t> &width_cnt,
            const std::vector<size_t> &no_hash,
            int32_t ch_cm_r, int32_t ch_cm_w);
  /**
   * @brief Release the pointer
   *
   */
  ~THD_CHHHUnivMon();
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
  int32_t getGlobalHash(const FlowKey<key_len> &flowkey, int32_t layer) const;
  int32_t getHashIdx(int32_t sketch_idx, int32_t hash_idx) const;
  int32_t getCHIdx(int32_t sketch_id, int32_t r, int32_t c) const;
  void updateSketch(int32_t sketch_idx, const FlowKey<key_len> &flowkey, T val);
  T querySketch(int32_t sketch_idx, const FlowKey<key_len> &flowkey) const;
  T estSketch(int32_t sketch_idx, const FlowKey<key_len> &flowkey) const;
};

} // namespace OmniSketch::Sketch

//-----------------------------------------------------------------------------
//
///                        Implementation of templated methods
//
//-----------------------------------------------------------------------------

namespace OmniSketch::Sketch {

template <int32_t key_len, int32_t no_layer, typename T, typename hash_t>
T THD_CHHHUnivMon<key_len, no_layer, T, hash_t>::estSketch(int32_t sketch_idx, const FlowKey<key_len> &flowkey) const{
  T values[depth];
  for (int i = 0; i < depth; ++i) {
    int32_t idx = CS_hash_fns[sketch_idx][i](flowkey) % width[sketch_idx];
    int32_t chIdx = getCHIdx(sketch_idx, i, idx);
    values[i] = ch->getEstCnt(chIdx) *
                (static_cast<int>(CS_update_hash_fns[sketch_idx][i](flowkey) & 1) * 2 - 1);
  }
  std::sort(values, values + depth);
  if (!(depth & 1)) { // even
    return std::abs((values[depth / 2 - 1] + values[depth / 2]) / 2);
  } else { // odd
    return std::abs(values[depth / 2]);
  }
}

template <int32_t key_len, int32_t no_layer, typename T, typename hash_t>
T THD_CHHHUnivMon<key_len, no_layer, T, hash_t>::querySketch(int32_t sketch_idx, const FlowKey<key_len> &flowkey) const{
  T values[depth];
  for (int i = 0; i < depth; ++i) {
    int32_t idx = CS_hash_fns[sketch_idx][i](flowkey) % width[sketch_idx];
    int32_t chIdx = getCHIdx(sketch_idx, i, idx);
    values[i] = ch->getCnt(chIdx) *
                (static_cast<int>(CS_update_hash_fns[sketch_idx][i](flowkey) & 1) * 2 - 1);
  }
  std::sort(values, values + depth);
  if (!(depth & 1)) { // even
    return std::abs((values[depth / 2 - 1] + values[depth / 2]) / 2);
  } else { // odd
    return std::abs(values[depth / 2]);
  }
}

template <int32_t key_len, int32_t no_layer, typename T, typename hash_t>
void THD_CHHHUnivMon<key_len, no_layer, T, hash_t>::updateSketch(int32_t sketch_idx, const FlowKey<key_len> &flowkey, T val){
  for(int i = 0; i < depth; i++){
    int32_t idx = CS_hash_fns[sketch_idx][i](flowkey) % width[sketch_idx];
    T update_val = val * (static_cast<int>(CS_update_hash_fns[sketch_idx][i](flowkey) & 1) * 2 - 1);
    ch->updateCnt(getCHIdx(sketch_idx, i, idx), update_val);
  }
}

template <int32_t key_len, int32_t no_layer, typename T, typename hash_t>
int32_t THD_CHHHUnivMon<key_len, no_layer, T, hash_t>::getCHIdx(int32_t sketch_id, int32_t r, int32_t c) const{
  return (width_idx[sketch_id] + c) + r * total_length;
}

template <int32_t key_len, int32_t no_layer, typename T, typename hash_t>
int32_t THD_CHHHUnivMon<key_len, no_layer, T, hash_t>::getGlobalHash(const FlowKey<key_len> &flowkey, int32_t layer) const{
  if(layer == 0){
    return 1;
  }
  return (global_hash_fns[layer - 1](flowkey) & 1);
}

template <int32_t key_len, int32_t no_layer, typename T, typename hash_t>
THD_CHHHUnivMon<key_len, no_layer, T, hash_t>::THD_CHHHUnivMon(int32_t depth_, int32_t width_,
                                     int32_t log_n, int32_t heap_size,
                                     double SSalpha, 
                                     double cnt_no_ratio,
                                     const std::vector<size_t> &width_cnt,
                                     const std::vector<size_t> &no_hash,
                                     int32_t ch_cm_r, int32_t ch_cm_w)
    : depth(depth_), logn(log_n),
      width_cnt(width_cnt), no_hash(no_hash)  {

  global_hash_fns = new hash_t[logn - 1];
  width_ = (Util::NextPrime(width_));

  CS_hash_fns = new hash_t*[logn];
  CS_update_hash_fns = new hash_t*[logn];

  int32_t running_idx = 0;
  width = new int32_t[logn];
  width_idx = new int32_t[logn];

  for(int i = 0; i < logn; i++){
    CS_hash_fns[i] = new hash_t[depth];
    CS_update_hash_fns[i] = new hash_t[depth];
  }

  for(int i = 0; i < logn; i++){
    width[i] = Util::NextPrime(width_);
    width_idx[i] = running_idx;
    running_idx += width[i];
    width_ = std::max(1, width_ / 2);
  }

  total_length = running_idx;

  // check ratio
  if (cnt_no_ratio <= 0.0 || cnt_no_ratio >= 1.0) {
    throw std::out_of_range("Out of Range: Ratio of #counters of adjacent "
                            "layers in CH should be in (0, 1), but got " +
                            std::to_string(cnt_no_ratio) + " instead.");
  }
  // prepare no_cnt
  no_cnt.push_back(static_cast<size_t>(this->depth) * running_idx);
  for (int32_t i = 1; i < no_layer; ++i) {
    size_t last_layer = no_cnt.back();
    no_cnt.push_back(Util::NextPrime(std::ceil(last_layer * cnt_no_ratio)));
  }
  // CH
  ch = new CounterHierarchy<no_layer, T, hash_t>(no_cnt, this->width_cnt,
                                                 this->no_hash, true, true, 
                                                 ch_cm_r, ch_cm_w);

  flows = new Data::Estimation<key_len>[logn];

  HHHeaps = new StreamSummary<key_len, T, hash_t>[logn];
  for(int i = 0; i < logn; i++){
    HHHeaps[i].init(heap_size, SSalpha);
  }
  sum = new T[logn];
}

template <int32_t key_len, int32_t no_layer, typename T, typename hash_t>
THD_CHHHUnivMon<key_len, no_layer, T, hash_t>::~THD_CHHHUnivMon() {
  if (global_hash_fns)
    delete[] global_hash_fns;
  for(int i = 0; i < depth; i++){
    delete[] CS_hash_fns[i];
    delete[] CS_update_hash_fns[i];
  }
  delete[] CS_hash_fns;
  delete[] CS_update_hash_fns;
  delete[] width;
  delete[] width_idx;
  delete[] ch;
  if (flows)
    delete[] flows;
  for(int32_t i = 0; i < logn; i++){
    HHHeaps[i].destroy();
  }
  delete[] HHHeaps;
  delete[] sum;
}

template <int32_t key_len, int32_t no_layer, typename T, typename hash_t>
Data::Estimation<key_len, T> THD_CHHHUnivMon<key_len, no_layer, T, hash_t>::getHeavyHitter(double val_threshold) const{
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

template <int32_t key_len, int32_t no_layer, typename T, typename hash_t>
void THD_CHHHUnivMon<key_len, no_layer, T, hash_t>::update(const FlowKey<key_len> &flowkey,
                                         T val) {
  for (int32_t i = 0; i < logn; ++i) {
    if (getGlobalHash(flowkey, i)) {
      sum[i] += val;
      updateSketch(i, flowkey, val);
      T est = estSketch(i, flowkey);
      HHHeaps[i].insert(flowkey, est);
    } else
      break;
  }
}

template <int32_t key_len, int32_t no_layer, typename T, typename hash_t>
T THD_CHHHUnivMon<key_len, no_layer, T, hash_t>::query(const FlowKey<key_len> &flowkey) const {
  int level;
  for(level = 0; level < logn; level++){
    if(!getGlobalHash(flowkey, level)){
      break;
    }
  }
  level--;
  T ret = querySketch(level, flowkey);
  for(int i = level - 1; i >= 0; i--){
    ret = 2 * ret - querySketch(i, flowkey);
  }
  return ret;
}

template <int32_t key_len, int32_t no_layer, typename T, typename hash_t>
size_t THD_CHHHUnivMon<key_len, no_layer, T, hash_t>::size() const {
  ch->print_rate("HHUnivMon CH");
  size_t heap_size = 0;
  for(int i = 0; i < logn; i++){
    heap_size += HHHeaps[i].memory_size();
  }
  return sizeof(*this) + 
         sizeof(hash_t) * (logn + 2 * logn * depth - 1) + 
         sizeof(int32_t) * 2 * logn + 
         heap_size + 
         ch->size();
}

template <int32_t key_len, int32_t no_layer, typename T, typename hash_t>
void THD_CHHHUnivMon<key_len, no_layer, T, hash_t>::clear() {
  // std::fill(counter[0], counter[0] + depth * width, 0);
}

} // namespace OmniSketch::Sketch


