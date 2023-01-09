/**
 * @file SALSACU.h
 * @author dromniscience (you@domain.com)
 * @brief Implementation of SALSA(Tango algorithm)
 *
 * @copyright Copyright (c) 2022
 *
 */
#pragma once
// #define DEBUG

#include <common/hash.h>
#include <common/sketch.h>

namespace OmniSketch::Sketch {
/**
 * @brief SALSA CU Sketch
 *
 * @tparam key_len  length of flowkey
 * @tparam T        type of the counter
 * @tparam hash_t   hashing class
 */
template <int32_t key_len, typename T, typename hash_t = Hash::AwareHash>
class SALSACU : public SketchBase<key_len, T> {
private:
  static int32_t maxCounterLen;
  static uint32_t upperBoundofUnitCounter;

  int32_t depth;
  int32_t width;
  hash_t *hash_fns;
  uint8_t **counter;
  uint8_t **bitMap;

  #ifdef DEBUG
  int32_t oriLen;
  T **oriCounter;
  #endif

  SALSACU(const SALSACU &) = delete;
  SALSACU(SALSACU &&) = delete;

  bool getFlagBit(const int32_t& rowIdx, const int32_t& idx) const;
  void setFlagBit(const int32_t& rowIdx, const int32_t& idx);
  void getBoundary(const int32_t& rowIdx, const int32_t& idx, int32_t& left, int32_t& right) const;
  void changeCounterVal(const int32_t& rowIdx, const int32_t& idx, T val);
  T getCounterVal(const int32_t& rowIdx, const int32_t& left, const int32_t& right) const;
  void setCounterVal(const int32_t& rowIdx, const int32_t& left, const int32_t& right, T val);
  void mergeCounter(const int32_t& rowIdx, const int32_t& left, const int32_t& right, T val);

public:
  /**
   * @brief Construct by specifying depth and width
   *
   */
  SALSACU(int32_t depth_, int32_t width_);
  /**
   * @brief Release the pointer
   *
   */
  ~SALSACU();
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
};

} // namespace OmniSketch::Sketch

//-----------------------------------------------------------------------------
//
///                        Implementation of templated methods
//
//-----------------------------------------------------------------------------

namespace OmniSketch::Sketch {

template <int32_t key_len, typename T, typename hash_t>
int32_t SALSACU<key_len, T, hash_t>::maxCounterLen = sizeof(T) / sizeof(uint8_t);

template <int32_t key_len, typename T, typename hash_t>
uint32_t SALSACU<key_len, T, hash_t>::upperBoundofUnitCounter = (1 << 8);

template <int32_t key_len, typename T, typename hash_t>
bool SALSACU<key_len, T, hash_t>::getFlagBit(const int32_t& rowIdx, const int32_t& idx) const{
    return (bitMap[rowIdx][idx / 8] >> (idx % 8)) & 1;
}

template <int32_t key_len, typename T, typename hash_t>
void SALSACU<key_len, T, hash_t>::setFlagBit(const int32_t& rowIdx, const int32_t& idx){
    bitMap[rowIdx][idx / 8] |= (1 << (idx % 8));
}

template <int32_t key_len, typename T, typename hash_t>
void SALSACU<key_len, T, hash_t>::getBoundary(const int32_t& rowIdx, const int32_t& idx, int32_t& left, int32_t& right) const {
    left = idx;
    right = idx;
    while(left >= 1 && getFlagBit(rowIdx, left - 1)){
        left--;
    }
    while(right < width - 1 && getFlagBit(rowIdx, right)){
        right++;
    }
    assert(right - left + 1 <= maxCounterLen);
}

template <int32_t key_len, typename T, typename hash_t>
void SALSACU<key_len, T, hash_t>::changeCounterVal(const int32_t& rowIdx, const int32_t& idx, T val){
    int32_t left, right;
    getBoundary(rowIdx, idx, left, right);
    int32_t len = right - left + 1;
    T counterVal = val;
    if(counterVal < (1 << (len << 3))){
        // set counter val
        setCounterVal(rowIdx, left, right, counterVal);
    } else{
        // need to merge counters
        mergeCounter(rowIdx, left, right, counterVal);
    }
}

template <int32_t key_len, typename T, typename hash_t>
void SALSACU<key_len, T, hash_t>::setCounterVal(const int32_t& rowIdx, const int32_t& left, const int32_t& right, T val){
    // if(left <= 450795 && right >= 450795){
    //     printf("falty link caught! update value is %d, row is %d, left is %d, right is %d, left bit is %d, right bit is %d\n", val, rowIdx, left, right, (int)getFlagBit(rowIdx, left), (int)getFlagBit(rowIdx, right));
    // }
    int32_t runner = right;
    while(runner >= left){
        counter[rowIdx][runner] = val % upperBoundofUnitCounter;
        val >>= 8;
        runner--;
    }
    assert(val == 0);
}

template <int32_t key_len, typename T, typename hash_t>
T SALSACU<key_len, T, hash_t>::getCounterVal(const int32_t& rowIdx, const int32_t& left, const int32_t& right) const{
    uint64_t ret = counter[rowIdx][left];
    int32_t runner = left + 1;
    while(runner <= right){
        ret = (ret << 8) + counter[rowIdx][runner];
        runner++;
    }
    return (T)ret;
}

template <int32_t key_len, typename T, typename hash_t>
void SALSACU<key_len, T, hash_t>::mergeCounter(const int32_t& rowIdx, const int32_t& left, const int32_t& right, T val){
    int32_t leftPower = (left) & (- (left));
    int32_t rightPower = (right + 1) & (-(right + 1));
    if(left == 0 || leftPower > rightPower){
        // merge with the counter on the right
        int32_t neighborLeft, neighborRight;
        getBoundary(rowIdx, right + 1, neighborLeft, neighborRight);
        T newVal = std::max(val, getCounterVal(rowIdx, neighborLeft, neighborRight));
        int32_t newLeft = left, newRight = neighborRight;
        setFlagBit(rowIdx, right);
        setCounterVal(rowIdx, newLeft, newRight, newVal);
    } else{
        assert(leftPower != rightPower);
        // merge with the counter on the left
        int32_t neighborLeft, neighborRight;
        getBoundary(rowIdx, left - 1, neighborLeft, neighborRight);
        T newVal = std::max(val, getCounterVal(rowIdx, neighborLeft, neighborRight));
        int32_t newLeft = neighborLeft, newRight = right;
        setFlagBit(rowIdx, left - 1);
        setCounterVal(rowIdx, newLeft, newRight, newVal);
    }
}

template <int32_t key_len, typename T, typename hash_t>
SALSACU<key_len, T, hash_t>::SALSACU(int32_t depth_, int32_t width_)
    : depth(depth_), width(maxCounterLen * Util::NextPrime((int32_t)(width_ / maxCounterLen))){

    hash_fns = new hash_t[depth];
    counter = new uint8_t* [depth];
    counter[0] = new uint8_t[depth * width]();
    for(int32_t i = 1; i < depth; ++i){
        counter[i] = counter[i - 1] + width;
    }
    int32_t counterNo = (width + 7) / 8;
    bitMap = new uint8_t* [depth];
    bitMap[0] = new uint8_t[depth * counterNo]();
    for(int32_t i = 1; i < depth; ++i){
        bitMap[i] = bitMap[i - 1] + counterNo;
    }
    #ifdef DEBUG
    oriLen = width / maxCounterLen;
    assert(oriLen * maxCounterLen == width);
    oriCounter = new T*[depth];
    oriCounter[0] = new T[oriLen * depth]();
    for(int32_t i = 1; i < depth; ++i){
        oriCounter[i] = oriCounter[i - 1] + oriLen;
    }
    #endif
}

template <int32_t key_len, typename T, typename hash_t>
SALSACU<key_len, T, hash_t>::~SALSACU(){
    delete[] hash_fns;
    delete[] counter[0];
    delete[] counter;
    delete[] bitMap[0];
    delete[] bitMap;
    #ifdef DEBUG
    delete[] oriCounter[0];
    delete[] oriCounter;
    #endif
}

template <int32_t key_len, typename T, typename hash_t>
void SALSACU<key_len, T, hash_t>::update(const FlowKey<key_len> &flowkey,
                                         T val){
    std::vector<T> allEst(depth);
    std::vector<int32_t> allIdx(depth);
    T min_val = std::numeric_limits<T>::max();
    #ifdef DEBUG
    T ori_min_val = std::numeric_limits<T>::max();
    std::vector<int32_t> allOriIdx(depth);
    #endif
    for(int32_t i = 0; i < depth; ++i){
        int32_t index = hash_fns[i](flowkey) % width;
        allIdx[i] = index;
        int32_t left, right;
        getBoundary(i, index, left, right);
        T est = getCounterVal(i, left, right);
        allEst[i] = est;
        min_val = std::min(min_val, val + est);
        #ifdef DEBUG
        int32_t oriIdx = index / maxCounterLen;
        allOriIdx[i] = oriIdx;
        ori_min_val = std::min(ori_min_val, oriCounter[i][oriIdx] + val);
        #endif
    }
    for(int32_t i = 0; i < depth; ++i){
        if(min_val > allEst[i]){
            changeCounterVal(i, allIdx[i], min_val);
        }
        #ifdef DEBUG
        int32_t oriIdx = allOriIdx[i];
        oriCounter[i][oriIdx] = std::max(oriCounter[i][oriIdx], ori_min_val);
        #endif
    }
}

template <int32_t key_len, typename T, typename hash_t>
T SALSACU<key_len, T, hash_t>::query(const FlowKey<key_len> &flowkey) const {
    #ifdef DEBUG
    // printf("query called!\n");
    // fflush(stdout);
    #endif
    T min_val = std::numeric_limits<T>::max();
    for(int32_t i = 0; i < depth; ++i){
        int32_t index = hash_fns[i](flowkey) % width;
        int32_t left, right;
        getBoundary(i, index, left, right);
        #ifndef DEBUG
        min_val = std::min(min_val, getCounterVal(i, left, right));
        #else
        T est = getCounterVal(i, left, right);
        int32_t oriIdx = index / maxCounterLen;
        if(est > oriCounter[i][oriIdx]){
            printf("idx: %d\n"
                   "oriIdx: %d\n"
                   "left: %d\n"
                   "right: %d\n"
                   "est: %d\n"
                   "CM: %d\n"
                   "row: %d\n", 
                   index, oriIdx, left, right, est, oriCounter[i][oriIdx], i);
            fflush(stdout);
        }
        assert(est <= oriCounter[i][oriIdx]);
        min_val = std::min(min_val, est);
        #endif
    }
    return min_val;
}

template <int32_t key_len, typename T, typename hash_t>
size_t SALSACU<key_len, T, hash_t>::size() const{
    int32_t counterNo = (width + 7) / 8;
    return sizeof(*this) + 
           sizeof(hash_t) * depth + 
           sizeof(uint8_t) * depth * width + 
           sizeof(uint8_t) * depth * counterNo;
}

template <int32_t key_len, typename T, typename hash_t>
void SALSACU<key_len, T, hash_t>::clear() {
    int32_t counterNo = (width + 7) / 8;
    std::fill(counter[0], counter[0] + depth * width, 0);
    std::fill(bitMap[0], bitMap[0] + depth * counterNo, 0);
}

} // namespace OmniSketch::Sketch
