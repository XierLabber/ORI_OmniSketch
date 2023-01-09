/**
 * @file PCMSketch.h
 * @author XierLabber (you@domain.com)
 * @brief Implementation of Pyramid Count Min Sketch
 *
 * @copyright Copyright (c) 2022
 *
 */
#pragma once

#include <common/hash.h>
#include <common/sketch.h>

namespace OmniSketch::Sketch {
/**
 * @brief Pyramid Count Min Sketch
 *
 * @tparam key_len  length of flowkey
 * @tparam T        type of the counter
 * @tparam hash_t   hashing class
 */
template <int32_t key_len, typename T, typename hash_t = Hash::AwareHash>
class PCMSketch : public SketchBase<key_len, T> {
private:
  int32_t pyramid_depth;
  int32_t d;                            // number of hash functions
  uint64_t** counter;                   // T is not used here for better performance

  int32_t word_index_size;              // The last word_index_size bits of a hash value
                                        // will be used to calculate the index of words

  int32_t lg_used_bits;                 // Pyramid structure will use (1 << lg_used_bits) bits
                                        // each layer
  int32_t counter_index_size;           // lg(word_length / counter_length)
  int32_t word_num;
  int32_t counter_num;

  int32_t total_counter_num;
  int32_t MY_MASK;
  int32_t HIGHER_MASK;

  int32_t MAX_HASH_NUM;

  hash_t hash_fns;

  T* original_counter;
  double ARE;
  int32_t query_time;

  PCMSketch(const PCMSketch &) = delete;
  PCMSketch(PCMSketch &&) = delete;

public:
  /**
   * @brief Construct by specifying depth and width
   *
   */
  PCMSketch(int32_t word_num_, int32_t hash_num, int32_t lg_used_bits_, int32_t pyramid_depth);
  /**
   * @brief Release the pointer
   *
   */
  ~PCMSketch();
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
   * @brief carry from the lower layer to the higher layer
   *
   */
  void carry(int32_t idx, T val);
  T get_value(int32_t idx) const;
};

} // namespace OmniSketch::Sketch

//-----------------------------------------------------------------------------
//
///                        Implementation of templated methods
//
//-----------------------------------------------------------------------------

namespace OmniSketch::Sketch {

template <int32_t key_len, typename T, typename hash_t>
PCMSketch<key_len, T, hash_t>::PCMSketch(int32_t word_num_, int32_t hash_num, int32_t lg_used_bits_, int32_t pyramid_depth_)
    : word_num(Util::NextPrime(word_num_)), d(hash_num),
    lg_used_bits(lg_used_bits_), total_counter_num(0),
    word_index_size(18), MAX_HASH_NUM(30), pyramid_depth(pyramid_depth_) {
  
  if (lg_used_bits >= 6 || lg_used_bits <= 0) {
  throw std::out_of_range("lg_used_bits: Should be in [1, 5], but got " +
                          std::to_string(lg_used_bits) + " instead.");
  }

  MY_MASK = ((1 << (1 << lg_used_bits)) - 1);
  HIGHER_MASK = ((1 << ((1 << lg_used_bits) - 2)) - 1);

  counter_index_size = 6 - lg_used_bits;
  counter_num = (word_num << counter_index_size);

  counter = new uint64_t*[pyramid_depth];
  for(int i = 0; i < pyramid_depth; i++){
    int cur_length = (word_num + (1 << i) - 1) >> i;
    counter[i] = new uint64_t[cur_length];
	memset(counter[i], 0, sizeof(uint64_t) * cur_length);
    total_counter_num += cur_length;
  }

  original_counter = new T[counter_num];
  ARE = 0;

}

template <int32_t key_len, typename T, typename hash_t>
PCMSketch<key_len, T, hash_t>::~PCMSketch() {
  delete hash_fns;
  for(int i = 0; i < d; i++){
    delete[] counter[i];
  }
  delete[] counter;
}

template <int32_t key_len, typename T, typename hash_t>
void PCMSketch<key_len, T, hash_t>::update(const FlowKey<key_len> &flowkey,
                                          T val) {
  T min_value = 1 << 30;

  T value[MAX_HASH_NUM];
  int index[MAX_HASH_NUM];
  int counter_offset[MAX_HASH_NUM];
	
  uint64_t hash_value = hash_fns(flowkey);
  int my_word_index = (hash_value & ((1 << word_index_size) - 1)) % word_num;
  hash_value >>= word_index_size;

  int flag = 0xFFFFFFF;

  for(int i = 0; i < d; i++)
  {
  	counter_offset[i] = (hash_value & 0xFFF) % (1 << counter_index_size);
  	index[i] = ((my_word_index << counter_index_size) + counter_offset[i]) % counter_num;
    original_counter[index[i]] += val;
  	hash_value >>= counter_index_size;
  
  	value[i] = (counter[0][my_word_index] >> (counter_offset[i] << lg_used_bits)) & MY_MASK;

	if(((flag >> counter_offset[i]) & 1) == 0)
		continue;

	flag &= (~(1 << counter_offset[i]));

    T new_val = value[i] + val;
	if (new_val > MY_MASK)
	{
        T tmp = (new_val & MY_MASK);
		counter[0][my_word_index] &= (~((uint64_t)MY_MASK << (counter_offset[i] << lg_used_bits)));
        counter[0][my_word_index] += ((uint64_t)tmp << (counter_offset[i] << lg_used_bits));
		carry(index[i], new_val >> (1 << lg_used_bits));
	}
	else
	{
		counter[0][my_word_index] += ((uint64_t)val << (counter_offset[i] << lg_used_bits));
	}
  }
  return;
}

template <int32_t key_len, typename T, typename hash_t>
T PCMSketch<key_len, T, hash_t>::query(const FlowKey<key_len> &flowkey) const {
  T min_value = 1 << 30;
  T original_min_value = 1 << 30;

  T value[MAX_HASH_NUM];
  int index[MAX_HASH_NUM];
  int counter_offset[MAX_HASH_NUM];
	
  uint64_t hash_value = hash_fns(flowkey);
  int my_word_index = (hash_value & ((1 << word_index_size) - 1)) % word_num;
  hash_value >>= word_index_size;

  for(int i = 0; i < d; i++)
  {
  	counter_offset[i] = (hash_value & 0xFFF) % (1 << counter_index_size);
  	index[i] = ((my_word_index << counter_index_size) + counter_offset[i]) % counter_num;
    
    T ori_val = original_counter[index[i]];
    original_min_value = ori_val < original_min_value ? ori_val : original_min_value;

  	hash_value >>= counter_index_size;

  	value[i] = (counter[0][my_word_index] >> (counter_offset[i] << lg_used_bits)) & MY_MASK;
  	value[i] += get_value(index[i]);
  	min_value = value[i] < min_value ? value[i] : min_value;
  }

  const_cast<PCMSketch<key_len, T, hash_t>*>(this)->query_time++;

  if(original_min_value != 0)
    const_cast<PCMSketch<key_len, T, hash_t>*>(this)->ARE += ((double)std::abs(min_value - original_min_value)) / original_min_value;

  return min_value;
}

template <int32_t key_len, typename T, typename hash_t>
size_t PCMSketch<key_len, T, hash_t>::size() const {
  printf("ARE: %lf\n", ARE / query_time);
  return sizeof(*this)                // instance
         + sizeof(uint64_t) * total_counter_num; // counter
}

template <int32_t key_len, typename T, typename hash_t>
void PCMSketch<key_len, T, hash_t>::clear() {
  for(int i = 0; i < pyramid_depth; i++){
    int cur_length = (word_num + (1 << i) - 1) >> i;
    std::fill(counter[i], counter[i] + cur_length, 0);
  }
}

template <int32_t key_len, typename T, typename hash_t>
void PCMSketch<key_len, T, hash_t>::carry(int32_t index, T val){
  int left_or_right;	
	
  T value;
  int word_index = index >> counter_index_size;
  int offset = index % (1 << counter_index_size);

  for(int i = 1; i < pyramid_depth; i++)
  {

	left_or_right = word_index & 1;
	word_index >>= 1;

	counter[i][word_index] |= ((uint64_t)0x1 << ((1 << lg_used_bits) - 2 + left_or_right + (offset << lg_used_bits)));
	value = (counter[i][word_index] >> (offset << lg_used_bits)) & MY_MASK;

	if(((value & HIGHER_MASK) + val) <= HIGHER_MASK)
	{
		counter[i][word_index] += ((uint64_t)val << (offset << lg_used_bits));
		return;
	}

    T new_val = (value & HIGHER_MASK) + val;
	counter[i][word_index] &= (~((uint64_t)HIGHER_MASK << (offset << lg_used_bits)));
    counter[i][word_index] += ((uint64_t)(new_val & HIGHER_MASK) << (offset << lg_used_bits));
    val = ((new_val) >> ((1 << lg_used_bits) - 2));
  }

  throw std::overflow_error(
            "Counter overflow at the last layer in PCM, overflow by " +
            std::to_string(val) + ".");
}

template <int32_t key_len, typename T, typename hash_t>
T PCMSketch<key_len, T, hash_t>::get_value(int32_t index) const
{
	int left_or_right;	
	int anti_left_or_right;
	
	int value;
    int word_index = index >> counter_index_size;
    int offset = index % (1 << counter_index_size);


	T high_value = 0;

	for(int i = 1; i < pyramid_depth; i++)
	{
		
		left_or_right = word_index & 1;
		anti_left_or_right = (left_or_right ^ 1);

		word_index >>= 1;

		value = (counter[i][word_index] >> (offset << lg_used_bits)) & MY_MASK;

		if(((value >> ((1 << lg_used_bits) - 2 + left_or_right)) & 1) == 0)
			return high_value;

		high_value += ((value & HIGHER_MASK) - ((value >> ((1 << lg_used_bits) - 2 + anti_left_or_right)) & 1)) * (1 << ((1 << lg_used_bits) - 2 + 2 * i));
	}
    throw std::overflow_error(
            "Query overflow at the last layer in PCM.");
}

} // namespace OmniSketch::Sketch
