/**
 * @file CHHeavyKeeper.h
 * @author XierLabber (you@domain.com)
 * @brief CHHeavy Keeper
 *
 * @copyright Copyright (c) 2022
 *
 */
#pragma once

#include <common/hash.h>
#include <common/hierarchy.h>
#include <common/sketch.h>

// #define DEBUG

namespace OmniSketch::Sketch {
/**
 * @brief CHHeavy Keeper
 *
 * @tparam key_len  length of flowkey
 * @tparam T        type of the counter
 * @tparam hash_t   hashing class
 */
template <int32_t key_len, int32_t no_layer, typename T,
          typename hash_t>
class CHHeavyKeeper : public SketchBase<key_len, T> {
private:

  static const int32_t MINBUCKET = -1;
  static const int32_t MAXBUCKET = 0xfffffff;

  static const int32_t TOMB = 0xfffffff;

  int32_t depth_;
  int32_t width_;

  double b_;

  hash_t *sketch_hash_fun_;
  hash_t fingerprint_hash_fun_;

  class hash_table_elem;
  class bucket_list_elem;

  T n_min_;

  std::vector<size_t> no_cnt;
  std::vector<size_t> width_cnt;
  std::vector<size_t> no_hash;

  CounterHierarchy<no_layer, T, hash_t> *ch;
  int16_t** FP;

  class StreamSummary;

  double hash_table_alpha;


  StreamSummary StreamSummary_;
  int32_t num_threshold_;

  CHHeavyKeeper(const CHHeavyKeeper &) = delete;
  CHHeavyKeeper(CHHeavyKeeper &&) = delete;
  CHHeavyKeeper &operator=(CHHeavyKeeper) = delete;

public:
  /**
   * @brief Construct by specifying depth, width, threshold size
   *        and the base used to calculate the probability of reduction
   *
   */
  CHHeavyKeeper(int32_t depth, int32_t width, int32_t num_threshold, double b, double hash_table_alpha, 
    double cnt_no_ratio, const std::vector<size_t> &width_cnt, const std::vector<size_t> &no_hash, 
    const int32_t ch_cm_r, const int32_t ch_cm_w);
  /**
   * @brief Release the pointer
   *
   */
  ~CHHeavyKeeper();
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
  int32_t getCHIdx(int32_t i, int32_t j) const{
    return i * width_ + j;
  }
};

} // namespace OmniSketch::Sketch

//-----------------------------------------------------------------------------
//
///                        Implementation of template methods
//
//-----------------------------------------------------------------------------

namespace OmniSketch::Sketch {

template <int32_t key_len, int32_t no_layer, typename T,
          typename hash_t>
class CHHeavyKeeper<key_len, no_layer, T, hash_t>::hash_table_elem{
public:
    bucket_list_elem* parent;
    hash_table_elem* next;
    hash_table_elem* prev;
    hash_table_elem* next_hash_elem;
    FlowKey<key_len> first;
    T second;
    hash_table_elem(): second(MINBUCKET), parent(NULL), next(NULL), next_hash_elem(NULL){};
};

template <int32_t key_len, int32_t no_layer, typename T,
          typename hash_t>
class CHHeavyKeeper<key_len, no_layer, T, hash_t>::bucket_list_elem{
public:
    T value;
    bucket_list_elem* prev;
    bucket_list_elem* next;
    hash_table_elem* child;
    bucket_list_elem(): value(MINBUCKET), prev(NULL), next(NULL), child(NULL){};
};

template <int32_t key_len, int32_t no_layer, typename T,
          typename hash_t>
class CHHeavyKeeper<key_len, no_layer, T, hash_t>::StreamSummary{
public:
    hash_table_elem** hash_table;
    bucket_list_elem* bucket_list_head;
    bucket_list_elem* bucket_list_tail;
    int32_t hash_table_length;
    hash_t hash_func;
    int32_t size_;
    void init(int32_t num_threshold, double hash_table_alpha)
    {
        hash_table_length = Util::NextPrime((int32_t)(num_threshold * hash_table_alpha));
        hash_table = new hash_table_elem* [hash_table_length];
        for(int i = 0; i < hash_table_length; i++)
        {
          hash_table[i] = NULL;
        }
        bucket_list_head = new bucket_list_elem;
        bucket_list_tail = new bucket_list_elem;
        bucket_list_head->next = bucket_list_tail;
        bucket_list_tail->prev = bucket_list_head;
        bucket_list_tail->value = MAXBUCKET;
        bucket_list_head->value = MINBUCKET;
        size_ = 0;
    }
    void clear()
    {
        for(int i = 0; i < hash_table_length; i++)
        {
          hash_table_elem* ptr = hash_table[i];
          while(ptr != NULL)
          {
            hash_table_elem* tmp = ptr;
            ptr = ptr->next_hash_elem;
            delete tmp;
          }
        }
        delete[] hash_table;
        bucket_list_elem* ptr = bucket_list_head;
        while(ptr != bucket_list_tail)
        {
          ptr = ptr->next;
          delete ptr->prev;
        }
        delete bucket_list_tail;
    }
    hash_table_elem* find(FlowKey<key_len> key)
    {
        int32_t index = hash_func(key) % hash_table_length;
        hash_table_elem* ptr = hash_table[index];
        while(ptr != NULL)
        {
          if(ptr->first == key)
          {
            return ptr;
          }
          ptr = ptr->next_hash_elem;
        }
        return NULL;
    }
    void create_new_bucket(bucket_list_elem* prev, hash_table_elem* h)
    {
      #ifdef DEBUG
      assert(prev != NULL);
      assert(h != NULL);
      #endif

        bucket_list_elem* tmp = new bucket_list_elem;
        tmp->next = prev->next;
        tmp->prev = prev;
        prev->next->prev = tmp;
        prev->next = tmp;
        tmp->child = h;
        tmp->value = h->second;
        h->next = h;
        h->prev = h;
        h->parent = tmp;
        #ifdef DEBUG
        assert(h->parent != NULL);
        #endif
        if(tmp->value == 1){
          // printf"VAL 1 CREATE! %p, h: %p\n",tmp, h);
        }
    }
    void insert_bucket_list_with_hash_table(hash_table_elem* h)
    {
      #ifdef DEBUG
      assert(h != NULL);
      #endif

        bucket_list_elem* ptr = bucket_list_head;
        while(ptr->value < h->second)
        {
            ptr = ptr->next;
        }
        #ifdef DEBUG
        assert(ptr != NULL);
        #endif
        if(ptr->value == h->second)
        {
            h->parent = ptr;
            h->next = ptr->child->next;
            h->prev = ptr->child;
            h->next->prev = h;
            ptr->child->next = h;
        }
        else
        {
            create_new_bucket(ptr->prev, h);
        }
      #ifdef DEBUG
      assert(h->parent != NULL);
      assert(h->parent->child != NULL);
      assert(h->parent->child->parent != NULL);
      #endif
    }
    void emplace(FlowKey<key_len> key, T val)
    {
        int32_t index = hash_func(key) % hash_table_length;
        hash_table_elem* ptr = hash_table[index];
        if(ptr == NULL)
        {
          hash_table[index] = new hash_table_elem;
          ptr = hash_table[index];
        }
        else
        {
          while(ptr->next_hash_elem != NULL)
          {
            ptr = ptr->next_hash_elem;
          }
          ptr->next_hash_elem = new hash_table_elem;
          ptr = ptr->next_hash_elem;
        }
        ptr->first = key;
        ptr->second = val;
        insert_bucket_list_with_hash_table(ptr);
        size_++;
        #ifdef DEBUG
        assert(ptr->parent != NULL);
        assert(ptr->parent->child->parent != NULL);
        #endif
    }
    void delete_bucket(bucket_list_elem* b)
    {
      if(b->value == 1){
        // printf"delete bucket: %p\n",b);
      }
      #ifdef DEBUG
      assert(b != NULL);
      assert(b->prev != NULL);
      assert(b->next != NULL);
      #endif
        // printf("B IS %p\n",b);
        // printf("B PREV IS %p\n",b->prev);
        // printf("B NEXT IS %p\n",b->next);
        b->next->prev = b->prev;
        b->prev->next = b->next;
        delete b;
    }
    void push_forward(hash_table_elem* h){
      if(h->parent->value == 1){
        // printf"push called, pushed: %p, parent: %p\n",h,h->parent);
      }

      bucket_list_elem* need_to_delete = NULL;
      // printf("PUSH: REACH HERE! h is %p\n",h);
      // printf("PUSH: REACH HERE! h parent is %p\n",h->parent);
      #ifdef DEBUG
      assert(h != NULL);
      assert(h->parent != NULL);
      assert(h->prev != NULL);
      assert(h->next != NULL);
      assert(h->parent->child != NULL);
      #endif

      if(h->parent->child == h){
        if(h->next != h){
          hash_table_elem* ptr = h->prev;
          ptr->next = h->next;
          h->next->prev = ptr;
          h->parent->child = ptr;
          #ifdef DEBUG
          assert(ptr->parent != NULL);
          assert(ptr->parent == h->parent);
          #endif
          if(h->parent->value == 1)
          {
            // printf"%p NEW CHILD: %p\n", ptr->parent, ptr);
          }
        }
        else{
          need_to_delete = h->parent;
          if(need_to_delete->value == 1){
            // printf"VAL 1 BUT DELETE! %p, cur h: %p\n",need_to_delete,h);
          }
        }
      }
      else{
        hash_table_elem* ptr = h->prev;
        h->next->prev = ptr;
        ptr->next = h->next;
      }

      #ifdef DEBUG
      if(need_to_delete != NULL){
        assert(h->parent->child->parent != NULL);
      }
      #endif

      // printf("PUSH: REACH HERE2!\n");
      // printf("H->PARENT->NEXT IS %p\n",h->parent->next);
      // printf("H->PARENT->NEXT->CHILD IS %p\n",h->parent->next->child);
      if(h->parent->next->value == h->second){
        bucket_list_elem* ptr = h->parent->next;
        #ifdef DEBUG
        assert(ptr != NULL);
        assert(ptr->child != NULL);
        assert(ptr->child->next != NULL);
        #endif
        h->parent = ptr;
        h->next = ptr->child->next;
        h->prev = ptr->child;
        h->next->prev = h;
        ptr->child->next = h;
      }
      else{
        // printf("HERE!?\n");
        create_new_bucket(h->parent, h);
      }
      if(need_to_delete != NULL){
        delete_bucket(need_to_delete);
      }

      #ifdef DEBUG
      assert(h->parent != NULL);
      #endif

      // printf("PUSH: REACH HERE3!\n");
    }
    void erase(hash_table_elem* h)
    {
      if(h->parent->value == 1){
        // printf"erase called! %p is erased, its parent is %p\n", h, h->parent);
      }
        if(h->parent == NULL)
        {
          // printf"H: %p, H->parent: %p,H->VAL: %d\n",h,h->parent,h->second);
        }
        #ifdef DEBUG
        assert(h->parent != NULL);
        #endif
        if(h == NULL)
        {
            return;
        }
        if(h->next == h)
        {
        // printf("REACH HERE! 0\n");
            delete_bucket(h->parent);
        }
        else
        {
        // printf("REACH HERE! 1\n");
            hash_table_elem* ptr = h->prev;
            h->next->prev = ptr;
            ptr->next = h->next;
            h->parent->child = h->next;
        }
        // printf("REACH HERE! 2\n");
        int32_t idx = hash_func(h->first) % hash_table_length;
        hash_table_elem* ptr = hash_table[idx];
        // printf("REACH HERE! 3\n");
        if(ptr == h)
        {
        // printf("REACH HERE! 4\n");
          hash_table[idx] = h->next_hash_elem;
        }
        else
        {
        // printf("REACH HERE! 5\n");
          while(ptr->next_hash_elem != h)
          {
            ptr = ptr->next_hash_elem;
          }
          ptr->next_hash_elem = h->next_hash_elem;
        }
        // printf("REACH HERE! 6\n");
        delete h;
        size_--;
        // printf("REACH HERE! 7\n");
    }
    hash_table_elem* get_least_elem()
    {
        bucket_list_elem* ptr = bucket_list_head->next;
        if(ptr->value == 1){
          // printf"get called! parent: %p\n",ptr);
        }
        if(ptr->value == MAXBUCKET)
        {
            return NULL;
        }
        #ifdef DEBUG
        assert(ptr->child != NULL);
        if(ptr->child->parent == NULL){
          // printf"err ptr: %p\nerr ptr val: %d, err child val: %d\nerror ptr->child: %p\n", ptr, ptr->value, ptr->child->second, ptr->child);
        }
        assert(ptr->child->parent != NULL);
        #endif
        return ptr->child;
    }
    int32_t size() const
    {
        return size_;
    }
};

template <int32_t key_len, int32_t no_layer, typename T,
          typename hash_t>
CHHeavyKeeper<key_len, no_layer, T, hash_t>::CHHeavyKeeper(
    int32_t depth, int32_t width, int32_t num_threshold, double b, double hash_table_alpha, 
    double cnt_no_ratio, const std::vector<size_t> &width_cnt, const std::vector<size_t> &no_hash, 
    const int32_t ch_cm_r, const int32_t ch_cm_w)
    : depth_(depth), width_(Util::NextPrime(width)), hash_table_alpha(hash_table_alpha), 
      num_threshold_(num_threshold), b_(b), n_min_(0),
      width_cnt(width_cnt), no_hash(no_hash) {

  sketch_hash_fun_ = new hash_t[depth_];

  // check ratio
  if (cnt_no_ratio <= 0.0 || cnt_no_ratio >= 1.0) {
    throw std::out_of_range("Out of Range: Ratio of #counters of adjacent "
                            "layers in CH should be in (0, 1), but got " +
                            std::to_string(cnt_no_ratio) + " instead.");
  }
  // prepare no_cnt
  no_cnt.push_back(static_cast<size_t>(this->depth_) * this->width_);
  for (int32_t i = 1; i < no_layer; ++i) {
    size_t last_layer = no_cnt.back();
    no_cnt.push_back(Util::NextPrime(std::ceil(last_layer * cnt_no_ratio)));
  }

  // CH
  ch = new CounterHierarchy<no_layer, T, hash_t>(no_cnt, this->width_cnt,
                                                 this->no_hash, false, true,
                                                 ch_cm_r, ch_cm_w);


  FP = new int16_t *[this->depth_];
  FP[0] = new int16_t[this->depth_ * this->width_];
  for(int i = 1; i < depth_; i++){
    FP[i] = FP[i - 1] + this->width_;
  }

  StreamSummary_.init(num_threshold, hash_table_alpha);

}

template <int32_t key_len, int32_t no_layer, typename T,
          typename hash_t>
CHHeavyKeeper<key_len, no_layer, T, hash_t>::~CHHeavyKeeper() {
  delete[] sketch_hash_fun_;

  delete[] FP[0];
  delete[] FP;

  delete[] ch;
  StreamSummary_.clear();
}

template <int32_t key_len, int32_t no_layer, typename T,
          typename hash_t>
size_t CHHeavyKeeper<key_len, no_layer, T, hash_t>::size() 
  const{
  ch->print_rate("");
  return depth_ * width_ * (sizeof(uint16_t))
         + ch->size()
         + depth_ * sizeof(hash_t)
         + sizeof(CHHeavyKeeper<key_len, no_layer, hash_t, T>)
         + Util::NextPrime((int)(num_threshold_ * hash_table_alpha)) * sizeof(hash_table_elem*)
         + num_threshold_ * (4 * sizeof(void *) + sizeof(T) + sizeof(FlowKey<key_len>));
}

template <int32_t key_len, int32_t no_layer, typename T,
          typename hash_t>
void CHHeavyKeeper<key_len, no_layer, T, hash_t>::update(
    const FlowKey<key_len> &flowkey, T val) {
  auto iter = StreamSummary_.find(flowkey);
  bool flag = (iter != NULL);

  int32_t EstimatedCount = MINBUCKET;
  int16_t FlowFP = fingerprint_hash_fun_(flowkey);

  bool done = false;

  for (int i = 0; i < depth_; i++) {
      int32_t index = sketch_hash_fun_[i](flowkey) % width_;
      int32_t chIdx = getCHIdx(i, index);
      int32_t counterC = ch->getEstCnt(chIdx);
      if (counterC > 0 && (flag || counterC < n_min_) && FP[i][index] == FlowFP) {
        counterC += val;
        EstimatedCount = (counterC < EstimatedCount) ? EstimatedCount : counterC;
        ch->updateCnt(chIdx, val);
        done = true;
      }
  }

  if (!done) {
    for (int i = 0; i < depth_; i++) {
      int32_t index = sketch_hash_fun_[i](flowkey) % width_;
      int32_t chIdx = getCHIdx(i, index);
      if (ch->getEstCnt(chIdx) == 0) {
        ch->updateCnt(chIdx, val);
        FP[i][index] = FlowFP;
        EstimatedCount = 1;
        done = true;
        break;
      }
    }
  }

  if (!done) {
    int32_t minC = ch->getEstCnt(getCHIdx(0, sketch_hash_fun_[0](flowkey) % width_));
    int32_t minCounterID = 0;
    for (int i = 1; i < depth_; i++) {
      int32_t index = sketch_hash_fun_[i](flowkey) % width_;
      int32_t chIdx = getCHIdx(i, index);
      int32_t counterC = ch->getEstCnt(chIdx);
      if (counterC < minC) {
        minC = counterC;
        minCounterID = i;
      }
    }
    int32_t minIndex = sketch_hash_fun_[minCounterID](flowkey) % width_;
    if (((double)rand() / RAND_MAX) < (pow(b_, -minC))) {
      int32_t chIdx = getCHIdx(minCounterID, minIndex);
      bool tmp = (ch->getEstCnt(chIdx) <= val);
      ch->updateCnt(chIdx, -val);
      if (tmp) {
        ch->resetCnt(chIdx, val);
        FP[minCounterID][minIndex] = FlowFP;
        EstimatedCount = 1;
      }
    }
  }

  if (EstimatedCount > 0) {
    if (flag) {
      if(EstimatedCount > iter->second){
        iter->second = EstimatedCount;
        StreamSummary_.push_forward(iter);
      }
      if (n_min_ == iter->second) {
        auto min_iter = StreamSummary_.get_least_elem();
        n_min_ = (min_iter == NULL)? 0 : min_iter->second;
      }
    } else if (StreamSummary_.size() < num_threshold_) {
      StreamSummary_.emplace(flowkey, EstimatedCount);
      n_min_ = std::min(EstimatedCount, n_min_);
    } else if (EstimatedCount > n_min_) {
      auto min_iter = StreamSummary_.get_least_elem();
      StreamSummary_.emplace(flowkey, EstimatedCount);
      StreamSummary_.erase(min_iter);
      min_iter = StreamSummary_.get_least_elem();
      n_min_ = (min_iter == NULL)? 0 : min_iter->second;
    }
  }
}

template <int32_t key_len, int32_t no_layer, typename T,
          typename hash_t>
T CHHeavyKeeper<key_len, no_layer, T, hash_t>::findKth(
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

template <int32_t key_len, int32_t no_layer, typename T,
          typename hash_t>
Data::Estimation<key_len, T>
CHHeavyKeeper<key_len, no_layer, T, hash_t>::getTopK(int32_t k) const{
  T ValueArray[num_threshold_];
  int32_t ArrayLength = 0;
  for(int i=0; i < StreamSummary_.hash_table_length; i++)
  {
    T val = StreamSummary_.hash_table[i].second;
    if(val != MINBUCKET && val != TOMB)
    {
        ValueArray[ArrayLength++] = val;
    }
  }
  int32_t KthValue = findKth(ValueArray, k, ArrayLength);
  Data::Estimation<key_len, T> TopK;
  for(int i=0; i < StreamSummary_.hash_table_length; i++)
  {
    T val = StreamSummary_.hash_table[i].second;
    if(val != MINBUCKET && val != TOMB && val >= KthValue)
    {
        TopK[StreamSummary_.hash_table[i].first] = StreamSummary_.hash_table[i].second;
    }
  }
  return TopK;
}

template <int32_t key_len, int32_t no_layer, typename T,
          typename hash_t>
Data::Estimation<key_len, T>
CHHeavyKeeper<key_len, no_layer, T, hash_t>::getHeavyHitter(
    double val_threshold) const {
  Data::Estimation<key_len, T> heavy_Hitter;
  int found = 0;
  for(int i=0; i < StreamSummary_.hash_table_length; i++)
  {
    hash_table_elem* ptr = StreamSummary_.hash_table[i];
    while(ptr != NULL)
    {
      T val = ptr->second;
      if(val >= val_threshold)
      {
          heavy_Hitter[ptr->first] = ptr->second;
          found++;
      }
      ptr = ptr->next_hash_elem;
    }
  }
  return heavy_Hitter;
}

} // namespace OmniSketch