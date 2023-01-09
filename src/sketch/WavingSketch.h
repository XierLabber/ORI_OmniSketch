/**
 * @file WavingSketch.h
 * @author XierLabber (you@domain.com)
 * @brief Implementation of Waving Sketch
 *
 * @copyright Copyright (c) 2022
 *
 */
#pragma once

#include <common/hash.h>
#include <common/sketch.h>

// #define DO_NOT_CONSIDER_FLOWKEY_SIZE

namespace OmniSketch::Sketch {
/**
 * @brief Waving Sketch
 *
 * @tparam key_len  length of flowkey
 * @tparam T        type of the counter
 * @tparam hash_t   hashing class
 */

template <int32_t key_len, typename T, typename hash_t = Hash::AwareHash>
class WavingSketch : public SketchBase<key_len, T> {
private:

  int32_t counter_num;
  int32_t heavy_part_length;
  hash_t h;
  hash_t s_;

  class heavy_part;
  class counter;

  counter* WavingCounter;
public:
  /**
   * @brief Construct by specifying counter_num and heavy_part_length
   *
   */
  WavingSketch(int32_t bucket_num_, int32_t heavy_part_length_);
  /**
   * @brief Release the pointer
   *
   */
  ~WavingSketch();
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
  struct real_heavy_part{ T frequency; FlowKey<key_len> key;};
};

} // namespace OmniSketch::Sketch

//-----------------------------------------------------------------------------
//
///                        Implementation of templated methods
//
//-----------------------------------------------------------------------------

namespace OmniSketch::Sketch {

template <int32_t key_len, typename T, typename hash_t>
class WavingSketch<key_len, T, hash_t>::heavy_part{
public:
    T frequency;
    FlowKey<key_len> key;
    uint8_t flag;
    heavy_part(): frequency(0), flag(1){}
    bool operator<(const heavy_part& h) const
    {
        return frequency < h.frequency;
    }
};

template <int32_t key_len, typename T, typename hash_t>
class WavingSketch<key_len, T, hash_t>::counter{
public:
    T counter;
    int32_t size_;
    heavy_part* heavy;
    void init(int32_t heavy_part_length_)
    {
        heavy = new heavy_part[heavy_part_length_]();
        counter = 0;
        size_ = 0;
    }
    void clear()
    {
        delete[] heavy;
    }
    int32_t smallestID()
    {
        int32_t ans = 0;
        T minn = heavy[0].frequency;
        for(int i = 1; i < size_; i++)
        {
            T fre = heavy[i].frequency;
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
    void insert(FlowKey<key_len> e, T val, uint8_t flag)
    {
        heavy[size_].key = e;
        heavy[size_].frequency = val;
        heavy[size_++].flag = flag;
    }
};

template <int32_t key_len, typename T, typename hash_t>
WavingSketch<key_len, T, hash_t>::WavingSketch(
    int32_t bucket_num_, int32_t heavy_part_length_):
    counter_num(Util::NextPrime(bucket_num_)), heavy_part_length(heavy_part_length_)
{
    WavingCounter = new counter[counter_num];
    for(int i = 0; i < counter_num; i++)
    {
        WavingCounter[i].init(heavy_part_length);
    }
}

template <int32_t key_len, typename T, typename hash_t>
WavingSketch<key_len, T, hash_t>::~WavingSketch()
{
    for(int i = 0; i < heavy_part_length; i++)
    {
        WavingCounter[i].clear();
    }
    delete[] WavingCounter;
}

template <int32_t key_len, typename T, typename hash_t>
size_t WavingSketch<key_len, T, hash_t>::size() 
  const{
#ifndef DO_NOT_CONSIDER_FLOWKEY_SIZE
  return counter_num * heavy_part_length * (sizeof(T) + sizeof(FlowKey<key_len>) + 0.125)
         + counter_num * sizeof(counter)
         + sizeof(WavingSketch<key_len, hash_t, T>);
#else
  return counter_num * sizeof(counter)
         + sizeof(WavingSketch<key_len, hash_t, T>)
         + counter_num * heavy_part_length * (8 * sizeof(T) + 1) / 8;
#endif
}

template <int32_t key_len, typename T, typename hash_t>
int32_t WavingSketch<key_len, T, hash_t>::s (const FlowKey<key_len> &key) const
{
    return 2 * (s_(key) % 2) - 1;
}

template <int32_t key_len, typename T, typename hash_t>
void WavingSketch<key_len, T, hash_t>::update(const FlowKey<key_len> &flowkey, T val)
{
    int32_t counter_idx = h(flowkey) % counter_num;
    counter* the_counter = WavingCounter + counter_idx;
    int32_t heavy_idx = the_counter->findID(flowkey);
    if(heavy_idx != -1)
    {
        heavy_part* the_heavy = the_counter->heavy + heavy_idx;
        the_heavy->frequency += val;
        if(the_heavy->flag == false)
        {
            the_counter->counter += s(flowkey);
        }
    }
    else if(the_counter->size_ < heavy_part_length)
    {
        the_counter->insert(flowkey, val, true);
    }
    else
    {
        int32_t sei = s(flowkey);
        T fi = sei * the_counter->counter;
        the_counter->counter += sei;
        int32_t smallest_idx = the_counter->smallestID();
        heavy_part* smallest_heavy = the_counter->heavy + smallest_idx;
        T fr = smallest_heavy->frequency;
        if(fr <= fi)
        {
            if(smallest_heavy->flag == true)
            {
                the_counter->counter = the_counter->counter + fr * s(smallest_heavy->key);
            }
            smallest_heavy->flag = false;
            smallest_heavy->frequency = fi + 1;
            smallest_heavy->key = flowkey;
        }
    }
}

template <int32_t key_len, typename T, typename hash_t>
T WavingSketch<key_len, T, hash_t>::query(const FlowKey<key_len> &flowkey) const {
    counter* the_counter = WavingCounter + (h(flowkey) % counter_num);
    int32_t flowID = the_counter->findID(flowkey);
    if(flowID == -1 || the_counter->heavy[flowID].flag == false)
    {
        return s(flowkey) * the_counter->counter;
    }
    else
    {
        return the_counter->heavy[flowID].frequency;
    }
}

template <int32_t key_len, typename T, typename hash_t>
Data::Estimation<key_len, T>
WavingSketch<key_len, T, hash_t>::getHeavyHitter(
    double val_threshold) const {
    Data::Estimation<key_len, T> heavy_Hitter;
    for(int i = 0; i < counter_num; i++)
    {
        counter* the_counter = WavingCounter + i;
        for(int j = 0; j < the_counter->size_; j++)
        {
            heavy_part* the_heavy = the_counter->heavy + j;
            if(the_heavy->frequency >= val_threshold)
            {
                heavy_Hitter[the_heavy->key] = the_heavy->frequency;
            }
        }
    }
    return heavy_Hitter;
}

} // namespace OmniSketch::Sketch