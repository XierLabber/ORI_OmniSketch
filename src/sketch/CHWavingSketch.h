/**
 * @file CHWavingSketch.h
 * @author XierLabber (you@domain.com)
 * @brief Implementation of CHWaving Sketch
 *
 * @copyright Copyright (c) 2022
 *
 */
#pragma once

#include <common/hash.h>
#include <common/hierarchy.h>
#include <common/sketch.h>

namespace OmniSketch::Sketch {
/**
 * @brief CHWaving Sketch
 *
 * @tparam key_len  length of flowkey
 * @tparam T        type of the counter
 * @tparam hash_t   hashing class
 */

template <int32_t key_len, int32_t no_layer, typename T,
          typename hash_t>
class CHWavingSketch : public SketchBase<key_len, T> {
private:

  int32_t counter_num;
  int32_t heavy_part_length;
  hash_t h;
  hash_t s_;

  class heavy_part;
  class counter;

  counter* CHWavingCounter;

  std::vector<size_t> counter_no_cnt;
  std::vector<size_t> counter_width_cnt;
  std::vector<size_t> counter_no_hash;
  CounterHierarchy<no_layer, T, hash_t>* counterCH;

  std::vector<size_t> heavy_no_cnt;
  std::vector<size_t> heavy_width_cnt;
  std::vector<size_t> heavy_no_hash;
  CounterHierarchy<no_layer, T, hash_t>* heavyCH;
  
public:
  /**
   * @brief Construct by specifying counter_num and heavy_part_length
   *
   */
  CHWavingSketch(int32_t bucket_num_, int32_t heavy_part_length_, 
             double counter_cnt_no_ratio_,
             const std::vector<size_t> &counter_width_cnt_,
             const std::vector<size_t> &counter_no_hash_, 
             double heavy_cnt_no_ratio_,
             const std::vector<size_t> &heavy_width_cnt_,
             const std::vector<size_t> &heavy_no_hash_, 
             const size_t counter_cm_r, 
             const size_t counter_cm_w,
             const size_t heavy_cm_r,
             const size_t heavy_cm_w);
  /**
   * @brief Release the pointer
   *
   */
  ~CHWavingSketch();
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
  /**
   * @brief get s(key)
   *
   */
  int32_t s(const FlowKey<key_len> &key) const;
  /**
   * @brief A struct to simulate real structures
   *
   */
  int32_t getHeavyIdx(int32_t counterIdx, int32_t heavyIdx) const;
};

} // namespace OmniSketch::Sketch

//-----------------------------------------------------------------------------
//
///                        Implementation of templated methods
//
//-----------------------------------------------------------------------------

namespace OmniSketch::Sketch {

template <int32_t key_len, int32_t no_layer, typename T,
          typename hash_t>
class CHWavingSketch<key_len, no_layer, T, hash_t>::heavy_part{
public:
    FlowKey<key_len> key;
    uint8_t flag;
    heavy_part(): flag(1){}
};

template <int32_t key_len, int32_t no_layer, typename T,
          typename hash_t>
class CHWavingSketch<key_len, no_layer, T, hash_t>::counter{
public:
    int32_t size_;
    heavy_part* heavy;
    void init(int32_t heavy_part_length_)
    {
        heavy = new heavy_part[heavy_part_length_]();
        size_ = 0;
    }
    void clear()
    {
        delete[] heavy;
    }
    int32_t smallestID(CounterHierarchy<no_layer, T, hash_t>* heavyCH_, int32_t heavy_part_length_, int32_t counterIdx_)
    {
        int32_t ans = 0;
        int32_t the_id = counterIdx_ * heavy_part_length_;
        T minn = heavyCH_->getEstCnt(the_id);
        // T minn = heavyCH_->getOriginalCnt(the_id);
        for(int i = 1; i < size_; i++)
        {
            T fre = heavyCH_->getEstCnt(the_id + i);
            // T fre = heavyCH_->getOriginalCnt(the_id + i);
            if(fre < minn)
            {
                minn = fre;
                ans = i;
            }
        }
        return ans;
    }
    int32_t findID(const FlowKey<key_len>& key)
    {
        for(int i = 0; i < size_; i++)
        {
            if(heavy[i].key == key)
            {
                return i;
            }
        }
        return -1;
    }
    void insert(FlowKey<key_len> e, T val, uint8_t flag, CounterHierarchy<no_layer, T, hash_t>* heavyCH_, int32_t heavy_part_length_, int32_t counterIdx_)
    {
        heavy[size_].key = e;
        heavyCH_->updateCnt(counterIdx_ * heavy_part_length_ + size_, val);
        heavy[size_++].flag = flag;
    }
};

template <int32_t key_len, int32_t no_layer, typename T,
          typename hash_t>
CHWavingSketch<key_len, no_layer, T, hash_t>::CHWavingSketch(
    int32_t bucket_num_, int32_t heavy_part_length_, 
    double counter_cnt_no_ratio_,
    const std::vector<size_t> &counter_width_cnt_,
    const std::vector<size_t> &counter_no_hash_, 
    double heavy_cnt_no_ratio_,
    const std::vector<size_t> &heavy_width_cnt_,
    const std::vector<size_t> &heavy_no_hash_, 
    const size_t counter_cm_r, 
    const size_t counter_cm_w,
    const size_t heavy_cm_r,
    const size_t heavy_cm_w):
    counter_num(Util::NextPrime(bucket_num_)), 
    heavy_part_length(heavy_part_length_),
    counter_width_cnt(counter_width_cnt_),
    counter_no_hash(counter_no_hash_),
    heavy_width_cnt(heavy_width_cnt_),
    heavy_no_hash(heavy_no_hash_)
{
  // check ratio
    if (counter_cnt_no_ratio_ <= 0.0 || counter_cnt_no_ratio_ >= 1.0) {
      throw std::out_of_range("Out of Range: Ratio of #counters of adjacent "
                              "layers in CH should be in (0, 1), but got " +
                              std::to_string(counter_cnt_no_ratio_) + " instead.");
    }
    if (heavy_cnt_no_ratio_ <= 0.0 || heavy_cnt_no_ratio_ >= 1.0) {
      throw std::out_of_range("Out of Range: Ratio of #counters of adjacent "
                              "layers in CH should be in (0, 1), but got " +
                              std::to_string(heavy_cnt_no_ratio_) + " instead.");
    }
    CHWavingCounter = new counter[counter_num];
    counter_no_cnt.push_back(static_cast<size_t>(this->counter_num));
    for (int32_t i = 1; i < no_layer; ++i) {
      size_t last_layer = counter_no_cnt.back();
      counter_no_cnt.push_back(Util::NextPrime(std::ceil(last_layer * counter_cnt_no_ratio_)));
    }
    counterCH = new CounterHierarchy<no_layer, T, hash_t>(counter_no_cnt, this->counter_width_cnt,
                                                          this->counter_no_hash, true, true, 
                                                          counter_cm_r, counter_cm_w);
    for(int i = 0; i < counter_num; i++)
    {
        CHWavingCounter[i].init(heavy_part_length);
    }
    heavy_no_cnt.push_back(static_cast<size_t>(this->heavy_part_length) * this->counter_num);
    for (int32_t i = 1; i < no_layer; ++i) {
      size_t last_layer = heavy_no_cnt.back();
      heavy_no_cnt.push_back(Util::NextPrime(std::ceil(last_layer * heavy_cnt_no_ratio_)));
    }
    heavyCH = new CounterHierarchy<no_layer, T, hash_t>(heavy_no_cnt, this->heavy_width_cnt,
                                                        this->heavy_no_hash, false, true,
                                                        heavy_cm_r, heavy_cm_w);
}

template <int32_t key_len, int32_t no_layer, typename T,
          typename hash_t>
CHWavingSketch<key_len, no_layer, T, hash_t>::~CHWavingSketch()
{
    for(int i = 0; i < heavy_part_length; i++)
    {
        CHWavingCounter[i].clear();
    }
    delete[] CHWavingCounter;
    delete[] counterCH;
    delete[] heavyCH;
}

template <int32_t key_len, int32_t no_layer, typename T,
          typename hash_t>
int32_t CHWavingSketch<key_len, no_layer, T, hash_t>::getHeavyIdx(int32_t counterIdx, int32_t heavyIdx) const{
    return counterIdx * heavy_part_length + heavyIdx;
}

template <int32_t key_len, int32_t no_layer, typename T,
          typename hash_t>
size_t CHWavingSketch<key_len, no_layer, T, hash_t>::size() 
  const{
  counterCH->print_rate("COUNTER CH");
  heavyCH->print_rate("HEAVY CH");
  return counter_num * heavy_part_length * (sizeof(FlowKey<key_len>) + 0.125)
         + counter_num * (sizeof(counter))
         + counterCH->size()
         + heavyCH->size()
         + sizeof(CHWavingSketch<key_len, no_layer, hash_t, T>);
}

template <int32_t key_len, int32_t no_layer, typename T,
          typename hash_t>
int32_t CHWavingSketch<key_len, no_layer, T, hash_t>::s (const FlowKey<key_len> &key) const
{
    return 2 * (s_(key) % 2) - 1;
}

template <int32_t key_len, int32_t no_layer, typename T,
          typename hash_t>
void CHWavingSketch<key_len, no_layer, T, hash_t>::update(const FlowKey<key_len> &flowkey, T val)
{
    int32_t counter_idx = h(flowkey) % counter_num;
    counter* the_counter = CHWavingCounter + counter_idx;
    int32_t heavy_idx = the_counter->findID(flowkey);
    if(heavy_idx != -1)
    {
        heavy_part* the_heavy = the_counter->heavy + heavy_idx;
        int32_t heavy_ch_idx = getHeavyIdx(counter_idx, heavy_idx);
        heavyCH->updateCnt(heavy_ch_idx, val);
        if(the_heavy->flag == false)
        {
            T res = s(flowkey);
            counterCH->updateCnt(counter_idx, res);
        }
    }
    else if(the_counter->size_ < heavy_part_length)
    {
        the_counter->insert(flowkey, val, true, heavyCH, heavy_part_length, counter_idx);
    }
    else
    {
        int32_t sei = s(flowkey);
        T fi = sei * counterCH->getEstCnt(counter_idx);
        // T fi = sei * counterCH->getOriginalCnt(counter_idx);
        counterCH->updateCnt(counter_idx, sei);
        int32_t smallest_idx = the_counter->smallestID(heavyCH, heavy_part_length, counter_idx);
        heavy_part* smallest_heavy = the_counter->heavy + smallest_idx;
        T fr = heavyCH->getEstCnt(getHeavyIdx(counter_idx, smallest_idx));
        // T fr = heavyCH->getOriginalCnt(getHeavyIdx(counter_idx, smallest_idx));
        if(fr <= fi)
        {
            if(smallest_heavy->flag == true)
            {
                T res = fr * s(smallest_heavy->key);
                counterCH->updateCnt(counter_idx, res);
            }
            smallest_heavy->flag = false;
            heavyCH->resetCnt(getHeavyIdx(counter_idx, smallest_idx), fi + 1);
            smallest_heavy->key = flowkey;
        }
    }
}

template <int32_t key_len, int32_t no_layer, typename T,
          typename hash_t>
T CHWavingSketch<key_len, no_layer, T, hash_t>::query(const FlowKey<key_len> &flowkey) const {
    int32_t counter_idx = (h(flowkey) % counter_num);
    counter* the_counter = CHWavingCounter + counter_idx;
    int32_t flowID = the_counter->findID(flowkey);
    if(flowID == -1 || the_counter->heavy[flowID].flag == false)
    {
        return s(flowkey) * counterCH->getCnt(counter_idx);
    }
    else
    {
        return heavyCH->getCnt(getHeavyIdx(counter_idx, flowID));
    }
}

template <int32_t key_len, int32_t no_layer, typename T,
          typename hash_t>
Data::Estimation<key_len, T>
CHWavingSketch<key_len, no_layer, T, hash_t>::getHeavyHitter(
    double val_threshold) const {
    Data::Estimation<key_len, T> heavy_Hitter;
    for(int i = 0; i < counter_num; i++)
    {
        int32_t the_id = getHeavyIdx(i, 0);
        counter* the_counter = CHWavingCounter + i;
        for(int j = 0; j < the_counter->size_; j++)
        {
            heavy_part* the_heavy = the_counter->heavy + j;
            if(heavyCH->getCnt(the_id + j) >= val_threshold)
            {
                heavy_Hitter[the_heavy->key] = heavyCH->getCnt(the_id + j);
            }
        }
    }
    return heavy_Hitter;
}

} // namespace OmniSketch::Sketch