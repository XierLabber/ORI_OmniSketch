/**
 * @file CHHashPipe.h
 * @author XierLabber (you@domain.com)
 * @brief CH Hash Pipe
 *
 * @copyright Copyright (c) 2022
 *
 */
#pragma once

#include <common/hash.h>
#include <common/hierarchy.h>
#include <common/sketch.h>

// #define DO_NOT_CONSIDER_FLOWKEY_SIZE
#define USE_ONE_CH

namespace OmniSketch::Sketch {
/**
 * @brief CH Hash Pipe
 *
 * @tparam key_len  length of flowkey
 * @tparam T        type of the counter
 * @tparam hash_t   hashing class
 */
template <int32_t key_len, int32_t no_layer, typename T,
          typename hash_t = Hash::AwareHash>
class CHHashPipe : public SketchBase<key_len, T> {
private:
  class Entry {
  public:
    FlowKey<key_len> flowkey;
    T val;
  };
  int32_t depth;
  int32_t normal_depth;
  int32_t ch_depth;
  int32_t width;
  hash_t *hash_fns;
  Entry **slots;

  std::vector<size_t> no_cnt;
  std::vector<size_t> width_cnt;
  std::vector<size_t> no_hash;

  FlowKey<key_len>** Keys;
#ifndef USE_ONE_CH
  CounterHierarchy<no_layer, T, hash_t> **ch;
#else
  CounterHierarchy<no_layer, T, hash_t> *ch;
#endif

  CHHashPipe(const CHHashPipe &) = delete;
  CHHashPipe(CHHashPipe &&) = delete;
  CHHashPipe &operator=(CHHashPipe) = delete;

public:
  /**
   * @brief Construct by specifying depth and width
   *
   */
  CHHashPipe(int32_t depth_, int32_t width_, double cnt_no_ratio,
             const std::vector<size_t> &width_cnt,
             const std::vector<size_t> &no_hash, 
             int32_t chcm_r,
             int32_t chcm_c, 
             int32_t chdepth_ = -1);
  /**
   * @brief Release the pointer
   *
   */
  ~CHHashPipe();
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
   * @brief Get Heavy Hitter
   * @param threshold A flowkey is a HH iff its counter `>= threshold`
   *
   */
  Data::Estimation<key_len, T> getHeavyHitter(double threshold) const override;
  /**
   * @brief Get sketch size
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
///                        Implementation of template methods
//
//-----------------------------------------------------------------------------

namespace OmniSketch::Sketch {

template <int32_t key_len, int32_t no_layer, typename T, typename hash_t>
CHHashPipe<key_len, no_layer, T, hash_t>::CHHashPipe(int32_t depth_, int32_t width_, double cnt_no_ratio,
                                                    const std::vector<size_t> &width_cnt_,
                                                    const std::vector<size_t> &no_hash_, 
                                                    int32_t chcm_r,
                                                    int32_t chcm_c, 
                                                    int32_t chdepth_)
    : depth(depth_), ch_depth((chdepth_ == -1)? depth_ : chdepth_), width(Util::NextPrime(width_)),
      width_cnt(width_cnt_), no_hash(no_hash_) {
  
  normal_depth = (depth - ch_depth); 
  hash_fns = new hash_t[depth];
  // Allocate continuous memory

  // check ratio
  if (cnt_no_ratio <= 0.0 || cnt_no_ratio >= 1.0) {
    throw std::out_of_range("Out of Range: Ratio of #counters of adjacent "
                            "layers in CH should be in (0, 1), but got " +
                            std::to_string(cnt_no_ratio) + " instead.");
  }
  // check depth
  if (ch_depth > depth || ch_depth < 0) {
    throw std::out_of_range("Out of Range: ch_depth of #counters of adjacent "
                            "layers should be in [0, depth), but got " +
                            std::to_string(ch_depth) + " instead.");
  }

  Keys = new FlowKey<key_len> *[ch_depth];
  Keys[0] = new FlowKey<key_len>[ch_depth * width]();
  for(int i = 1; i < ch_depth; i++){
    Keys[i] = Keys[i - 1] + width;
  }

#ifndef USE_ONE_CH
  // prepare no_cnt
  no_cnt.push_back(static_cast<size_t>(this->width));
  for (int32_t i = 1; i < no_layer; ++i) {
    size_t last_layer = no_cnt.back();
    no_cnt.push_back(Util::NextPrime(std::ceil(last_layer * cnt_no_ratio)));
  }

  ch = new CounterHierarchy<no_layer, T, hash_t> *[ch_depth];
  for(int i = 0; i < ch_depth; i++){
    ch[i] = new CounterHierarchy<no_layer, T, hash_t>(no_cnt, this->width_cnt,
                                                      this->no_hash);
  }
#else
  // prepare no_cnt
  no_cnt.push_back(static_cast<size_t>(this->width * this->ch_depth));
  for (int32_t i = 1; i < no_layer; ++i) {
    size_t last_layer = no_cnt.back();
    no_cnt.push_back(Util::NextPrime(std::ceil(last_layer * cnt_no_ratio)));
  }
  
  ch = new CounterHierarchy<no_layer, T, hash_t>(no_cnt, this->width_cnt, this->no_hash, false, 
                                                 true, chcm_r, chcm_c);
#endif

  // Allocate continuous memory
  if(normal_depth >= 1){
    slots = new Entry *[normal_depth];
    slots[0] = new Entry[normal_depth * width](); // Init with zero
    for (int32_t i = 1; i < normal_depth; ++i) {
      slots[i] = slots[i - 1] + width;
    }
  }
}

template <int32_t key_len, int32_t no_layer, typename T, typename hash_t>
CHHashPipe<key_len, no_layer, T, hash_t>::~CHHashPipe() {
  delete[] hash_fns;
  delete[] Keys[0];
  delete[] Keys;
  delete[] ch[0];
  delete[] ch;
}

template <int32_t key_len, int32_t no_layer, typename T, typename hash_t>
void CHHashPipe<key_len, no_layer, T, hash_t>::update(const FlowKey<key_len> &flowkey,
                                          T val) {
  if(normal_depth >= 1){
    // The first stage
    int idx = hash_fns[0](flowkey) % width;
    FlowKey<key_len> empty_key;
    FlowKey<key_len> c_key;
    T c_val;
    if (slots[0][idx].flowkey == flowkey) {
      // flowkey hit
      slots[0][idx].val += val;
      return;
    } else if (slots[0][idx].flowkey == empty_key) {
      // empty
      slots[0][idx].val = val;
      slots[0][idx].flowkey = flowkey;
      return;
    } else {
      // swap
      c_key = slots[0][idx].flowkey;
      c_val = slots[0][idx].val;
      slots[0][idx].flowkey = flowkey;
      slots[0][idx].val = val;
    }
    // Later stages
    for (int i = 1; i < normal_depth; ++i) {
      idx = hash_fns[i](c_key) % width;
      if (slots[i][idx].flowkey == c_key) {
        slots[i][idx].val += c_val;
        return;
      } else if (slots[i][idx].flowkey == empty_key) {
        slots[i][idx].flowkey = c_key;
        slots[i][idx].val = c_val;
        return;
      } else if (slots[i][idx].val < c_val) {
        // swap
        // slots[i][idx].swap(c_key);
        auto tmpkey = c_key;
        c_key = slots[i][idx].flowkey;
        slots[i][idx].flowkey = tmpkey;
        std::swap(c_val, slots[i][idx].val);
      }
    }
    for(int j = normal_depth; j < depth; j++){
      idx = hash_fns[j](c_key) % width;
      int i = j - normal_depth;
      #ifndef USE_ONE_CH
      if(Keys[i][idx] == c_key){
        ch[i]->updateCnt(idx, c_val);
        return;
      } else if(Keys[i][idx] == empty_key){
        Keys[i][idx] = c_key;
        ch[i]->updateCnt(idx, c_val);
        return;
      } else{
        T EstVal = ch[i]->getEstCnt(idx, 0);
        if(EstVal < c_val){
          auto tmpkey = c_key;
          c_key = Keys[i][idx];
          Keys[i][idx] = tmpkey;
          ch[i]->updateCnt(idx, c_val - EstVal);
          c_val = EstVal;
        }
      }
      #else
      int offset = i * width;
      if(Keys[i][idx] == c_key){
        ch->updateCnt(offset + idx, c_val);
        return;
      } else if(Keys[i][idx] == empty_key){
        Keys[i][idx] = c_key;
        ch->updateCnt(offset + idx, c_val);
        return;
      } else{
        T EstVal = ch->getEstCnt(offset + idx);
        if(EstVal < c_val){
          auto tmpkey = c_key;
          c_key = Keys[i][idx];
          Keys[i][idx] = tmpkey;
          ch->updateCnt(offset + idx, c_val - EstVal);
          c_val = EstVal;
        }
      }
      #endif
    }
  }
  else{
    // The first stage
    int idx = hash_fns[0](flowkey) % width;
    FlowKey<key_len> empty_key;
    FlowKey<key_len> c_key;
    T c_val;
    if (Keys[0][idx] == flowkey) {
      // flowkey hit
      ch->updateCnt(idx, val);
      return;
    } else if (Keys[0][idx] == empty_key) {
      // empty
      ch->updateCnt(idx, val);
      Keys[0][idx] = flowkey;
      return;
    } else {
      // swap
      c_key = Keys[0][idx];
      c_val = ch->getEstCnt(idx);
      Keys[0][idx] = flowkey;
      ch->updateCnt(idx, val - c_val);
    }
    for(int j = 1; j < depth; j++){
      idx = hash_fns[j](c_key) % width;
      int i = j;
      #ifndef USE_ONE_CH
      if(Keys[i][idx] == c_key){
        ch[i]->updateCnt(idx, c_val);
        return;
      } else if(Keys[i][idx] == empty_key){
        Keys[i][idx] = c_key;
        ch[i]->updateCnt(idx, c_val);
        return;
      } else{
        T EstVal = ch[i]->getEstCnt(idx, 0);
        if(EstVal < c_val){
          auto tmpkey = c_key;
          c_key = Keys[i][idx];
          Keys[i][idx] = tmpkey;
          ch[i]->updateCnt(idx, c_val - EstVal);
          c_val = EstVal;
        }
      }
      #else
      int offset = i * width;
      if(Keys[i][idx] == c_key){
        ch->updateCnt(offset + idx, c_val);
        return;
      } else if(Keys[i][idx] == empty_key){
        Keys[i][idx] = c_key;
        ch->updateCnt(offset + idx, c_val);
        return;
      } else{
        T EstVal = ch->getEstCnt(offset + idx);
        if(EstVal < c_val){
          auto tmpkey = c_key;
          c_key = Keys[i][idx];
          Keys[i][idx] = tmpkey;
          ch->updateCnt(offset + idx, c_val - EstVal);
          c_val = EstVal;
        }
      }
      #endif
    }
  }
}

template <int32_t key_len, int32_t no_layer, typename T, typename hash_t>
T CHHashPipe<key_len, no_layer, T, hash_t>::query(const FlowKey<key_len> &flowkey) const {
  T ret = 0;
  for (int i = 0; i < normal_depth; ++i) {
    int idx = hash_fns[i](flowkey) % width;
    if (slots[i][idx].flowkey == flowkey) {
      ret += slots[i][idx].val;
    }
  }
  #ifndef USE_ONE_CH
  for (int j = normal_depth; j < depth; ++j) {
    int i = j - normal_depth;
    int idx = hash_fns[j](flowkey) % width;
    if (Keys[i][idx] == flowkey) {
      ret += ch[i]->getCnt(idx);
    }
  }
  #else
  for (int j = normal_depth; j < depth; ++j) {
    int i = j - normal_depth;
    int idx = hash_fns[j](flowkey) % width;
    if (Keys[i][idx] == flowkey) {
      ret += ch->getCnt(i * width + idx);
    }
  }
  #endif
  return ret;
}

template <int32_t key_len, int32_t no_layer, typename T, typename hash_t>
Data::Estimation<key_len, T>
CHHashPipe<key_len, no_layer, T, hash_t>::getHeavyHitter(double threshold) const {
  Data::Estimation<key_len, T> heavy_hitters;
  std::set<FlowKey<key_len>> checked;
  for (int i = 0; i < normal_depth; ++i) {
    for (int j = 0; j < width; ++j) {
      const auto &flowkey = slots[i][j].flowkey;
      if (checked.find(flowkey) != checked.end()) {
        continue;
      }
      checked.insert(flowkey);
      auto estimate_val = query(flowkey);
      if (estimate_val >= threshold) {
        heavy_hitters[flowkey] = estimate_val;
      }
    }
  }
  for (int t = normal_depth; t < depth; ++t) {
    int i = t -  normal_depth;
    for (int j = 0; j < width; ++j) {
      const auto &flowkey = Keys[i][j];
      if (checked.find(flowkey) != checked.end()) {
        continue;
      }
      checked.insert(flowkey);
      auto estimate_val = query(flowkey);
      if (estimate_val >= threshold) {
        heavy_hitters[flowkey] = estimate_val;
      }
    }
  }
  return heavy_hitters;
}

template <int32_t key_len, int32_t no_layer, typename T, typename hash_t>
size_t CHHashPipe<key_len, no_layer, T, hash_t>::size() const {
  #ifndef USE_ONE_CH
  for(int i = normal_depth; i < depth; i++){
    printf("CH %d: ",i);
    ch[i - normal_depth]->print_rate("");
  }
  size_t ans = sizeof(*this)                           // instance
               + sizeof(hash_t) * depth                // hashing class
               + (sizeof(FlowKey<key_len>) + sizeof(T)) * normal_depth * width  // slots
               + sizeof(CounterHierarchy<no_layer, T, hash_t>*) * ch_depth
               + sizeof(FlowKey<key_len>*) * ch_depth
               + sizeof(FlowKey<key_len>) * ch_depth * width;
  for(int i = normal_depth; i < depth; i++){
    ans += ch[i - normal_depth]->size();
  }
  #else
  ch->print_rate("");
  
  size_t ans = sizeof(*this)                           // instance
               + sizeof(hash_t) * depth                // hashing class
               + (sizeof(FlowKey<key_len>) + sizeof(T)) * normal_depth * width  // slots
               + sizeof(CounterHierarchy<no_layer, T, hash_t>*) * ch_depth
               + sizeof(FlowKey<key_len>*) * ch_depth
               + (sizeof(FlowKey<key_len>)) * ch_depth * width
               + ch->size();
  #endif
  return ans;
}

template <int32_t key_len, int32_t no_layer, typename T, typename hash_t>
void CHHashPipe<key_len, no_layer, T, hash_t>::clear() {
  FlowKey<key_len> empty_key;
  for (int i = 0; i < normal_depth; ++i) {
    for (int j = 0; j < width; ++j) {
      slots[i][j].flowkey = empty_key;
      slots[i][j].val = 0;
    }
  }
  #ifndef USE_ONE_CH
  for (int t = normal_depth; t < depth; ++t) {
    int i = t - normal_depth;
    ch[i]->clear();
    for (int j = 0; j < width; ++j) {
      Keys[i][j] = empty_key;
    }
  }
  #else
  ch->clear();
  for (int t = normal_depth; t < depth; ++t) {
    int i = t - normal_depth;
    for (int j = 0; j < width; ++j) {
      Keys[i][j] = empty_key;
    }
  }
  #endif
}

} // namespace OmniSketch::Sketch