/**
 * @file CounterTree.h
 * @author XierLabber (you@domain.com)
 * @brief Heavy Keeper
 *
 * @copyright Copyright (c) 2022
 *
 */
#pragma once

#include <common/hash.h>
#include <common/sketch.h>

#include<stdlib.h>
#include<time.h>
#include<string.h>
#include<algorithm>

namespace OmniSketch::Sketch {
/**
 * @brief Counter Tree
 *
 * @tparam key_len  length of flowkey
 * @tparam T        type of the counter
 * @tparam hash_t   hashing class
 */
template <int32_t key_len, typename T, typename hash_t = Hash::AwareHash>
class CounterTree : public SketchBase<key_len, T> {
private:
    typedef uint16_t counter_t;
    int32_t b; // used bit num
    int32_t h; // height of the tree
    int32_t d; // degree of tree nodes
    int32_t m; // number of virtual leaves
    int32_t r; // number of buckets a flow has
    T n;       // number of packets
    int32_t bound_of_counters;
    int32_t counter_num;
    int32_t true_m; // number of true leaves

    counter_t* counter;
    T* real_val;
    uint8_t* flag;
#ifdef USE_RANDKEY
    int8_t** rand_key;
    hash_t hash_func;
#else
    hash_t* hash_func;
#endif

    void update_counter(int32_t idx, T val);
    int32_t get_parent(int32_t idx) const;
    int32_t get_left_son(int32_t idx) const;
    double get_estimated_val(int32_t idx) const;
    T count_sub_tree(int32_t root_idx, T pow) const;

public:
  /**
   * @brief Construct by specifying b, h, d, r and m
   *
   */
  CounterTree(int32_t b_, int32_t h_, int32_t d_, int32_t r_, int32_t m_);
  /**
   * @brief Release the pointer
   *
   */
  ~CounterTree();
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
};

}// namespace OmniSketch::Sketch

//-----------------------------------------------------------------------------
//
///                        Implementation of template methods
//
//-----------------------------------------------------------------------------

namespace OmniSketch::Sketch {

template <int32_t key_len, typename T, typename hash_t>
CounterTree<key_len, T, hash_t>::CounterTree(
    int32_t b_, int32_t h_, int32_t d_, int32_t r_, int32_t m_):
    b(b_), h(h_), d(d_), r(r_), m(Util::NextPrime(m_)), n(0){
    int32_t dh = 1;
    for(int i = 1; i < h; i++)
    {
        dh *= d;
    }
    int32_t running_m = (m + dh - 1) / dh;
    running_m *= dh;
    true_m = running_m;
    counter_num = running_m;
    for(int i = 1; i < h; i++)
    {
        running_m = (running_m + d - 1) / d;
        counter_num += running_m;
    }
    bound_of_counters = (1 << b);
    counter = new counter_t[counter_num];
    flag = new uint8_t[counter_num];
    real_val = new T[true_m];
    memset(counter, 0, counter_num * sizeof(counter_t));
    memset(flag, 0, counter_num * sizeof(uint8_t));
    memset(real_val, 0, true_m * sizeof(T));
#ifdef USE_RANDKEY
    rand_key = new int8_t *[r];
    rand_key[0] = new int8_t[r * key_len];
    for(int i = 1; i < r; i++)
    {
        rand_key[i] = rand_key[i - 1] + key_len;
    }

    std::srand((unsigned)time(NULL));

    for(int i = 0; i < r; i++)
    {
        for(int j = 0; j < key_len; j++)
        {
            rand_key[i][j] = rand() % (1 << (8 * sizeof(int8_t)));
        }
    }
#else
    hash_func = new hash_t[r];
#endif
}

template <int32_t key_len, typename T, typename hash_t>
CounterTree<key_len, T, hash_t>::~CounterTree(){
    delete[] counter;
    delete[] flag;
#ifdef USE_RANDKEY
    delete[] rand_key[0];
    delete[] rand_key;
#else
    delete[] hash_func;
#endif
    delete[] real_val;
}

template <int32_t key_len, typename T, typename hash_t>
int32_t CounterTree<key_len, T, hash_t>::get_parent(
    int32_t idx) const{
    int32_t ans = (idx / d) + true_m;
    return (ans >= counter_num)? -1 : ans;
}

template <int32_t key_len, typename T, typename hash_t>
int32_t CounterTree<key_len, T, hash_t>::get_left_son(
    int32_t idx) const{
    return (idx - true_m) * d;
}

template <int32_t key_len, typename T, typename hash_t>
double CounterTree<key_len, T, hash_t>::get_estimated_val(
    int32_t idx) const {
    T pow = 1;
    int32_t k = 1;
    int32_t root_idx = idx;
    while(flag[root_idx] == 1)
    //while(true)
    {
        int32_t fa = get_parent(root_idx);
        if(fa == -1)
        {
            break;
        }
        k = k * d;
        pow = pow * bound_of_counters;
        root_idx = fa;
    }
    return r * count_sub_tree(root_idx, pow) - (double)(n * k * r) / m;
    //return std::max(0.0, r * count_sub_tree(root_idx, pow) - (double)(n * k * r) / m);
}

template <int32_t key_len, typename T, typename hash_t>
void CounterTree<key_len, T, hash_t>::update_counter(int32_t idx, T val){
    T result = counter[idx] + val;
    if(result >= bound_of_counters)
    {
        counter[idx] = result % bound_of_counters;
        int32_t parent = get_parent(idx);
        if(parent != -1)
        {
            update_counter(parent, (result - counter[idx]) / bound_of_counters);
        }
        else
        {
            printf("WARNING! The last layer of counters overflows\n");
            exit(-1);
        }
        flag[idx] = 1;
    }
    else
    {
        counter[idx] = result;
    }
}

template <int32_t key_len, typename T, typename hash_t>
void CounterTree<key_len, T, hash_t>::update(
    const FlowKey<key_len> &flowkey, T val){
    n += val;
#ifdef USE_RANDKEY
    FlowKey<key_len> tmp_key = flowkey;
    tmp_key ^= (FlowKey<key_len>) (rand_key[rand() % r]);
    int32_t idx = hash_func(tmp_key) % m;
#else
    int32_t idx = hash_func[rand() % r](flowkey) % m;
#endif
    update_counter(idx, val);
    real_val[idx] += val;
}

template <int32_t key_len, typename T, typename hash_t>
T CounterTree<key_len, T, hash_t>::count_sub_tree(
    int32_t root_idx, T pow) const{
    if(root_idx < true_m)
    {
        if(pow != 1)
        {
            printf("WARNING! Wrong pow has been passed to count_sub_tree!\n");
            exit(0);
        }
        return counter[root_idx];
    }
    T ans = counter[root_idx];
    /*
    if(pow == -1)
    {
        pow = 1;
        int32_t running_idx = root_idx;
        while(running_idx >= m)
        {
            running_idx = get_left_son(running_idx);
            pow *= bound_of_counters;
        }
    }
    */
    ans *= pow;
    T reduced_pow = pow / bound_of_counters;
    int32_t ls = get_left_son(root_idx);
    for(int i = 0; i < d; i++)
    {
        ans += count_sub_tree(ls + i, reduced_pow);
    }
    return ans;
}

template <int32_t key_len, typename T, typename hash_t>
T CounterTree<key_len, T, hash_t>::query(const FlowKey<key_len> &flowkey) const{
    double ans = 0;
    for(int i = 0; i < r; i++)
    {
    #ifdef USE_RANDKEY
        FlowKey<key_len> tmp_key = flowkey;
        tmp_key ^= (FlowKey<key_len>)(rand_key[i]);
        int32_t idx = hash_func(tmp_key) % m;
    #else
        int32_t idx = hash_func[i](flowkey) % m;
    #endif
        ans += get_estimated_val(idx);
    }
    return ans / r;
}

template <int32_t key_len, typename T, typename hash_t>
size_t CounterTree<key_len, T, hash_t>::size() const{
#ifdef MY_DEBUG
    int32_t overflow_num = 0;
    int32_t correct_num = 0;
    double total_est_val = 0;
    T total_val = 0;
    double are = 0;
    for(int i = 0; i < m; i++)
    {
        if(flag[i] == 1)
        {
            overflow_num++;
        }
        double est_tmp = get_estimated_val(i);
        if((T)est_tmp == real_val[i])
        {
            correct_num++;
        }
        total_est_val += est_tmp;
        total_val += real_val[i];
        double r = (double)est_tmp / real_val[i];
        are += (r > 1)? (r - 1) : (1 - r);
    }
    T total_rec_val = 0;
    for(int i = 0; i < counter_num; i++)
    {
        int idx = i;
        int pow = 1;
        while(idx >= true_m)
        {
            idx = get_left_son(idx);
            pow = pow * bound_of_counters;
        }
        total_rec_val += counter[i] * pow;
    }
    int32_t my_k = 1;
    for(int32_t i = 1; i < h; i++)
    {
        my_k *= d;
    }
    are /= m;
    printf("RECORD %d PACKETS, THERE ARE %d PACKETS FOR REAL\n", total_rec_val, total_val);
    printf("OVERFLOW RATE: %lf\n",(double)overflow_num / m);
    printf("CORRECT RATE: %lf\n",(double)correct_num / m);
    printf("AVERAGE RATIO: %lf\n",(double)total_est_val / total_val);
    printf("TOTAL EST VAL: %lf\n",total_est_val);
    printf("n * k = %d\n", n * my_k);
    printf("EST ARE: %lf\n",are);

    printf("=============================================\n");
    int32_t l = 0;
    int32_t r = true_m - 1;
    int32_t k = 1;
    for(int j = l; j <= r; j++)
    {
        printf("EST[%d]: %lf\n", j, get_estimated_val(j));
    }
    for(int j = l; j <= r; j++)
    {
        for(int t = 0; t < k; t++)
        {
            printf(" ");
        }
        printf("%d",real_val[j]);
        for(int t = 0; t < k; t++)
        {
            printf(" ");
        }
    }
    printf("\n");
    for(int i = 0; i < h; i++)
    {
        for(int j = l; j <= r; j++)
        {
            for(int t = 0; t < k; t++)
            {
                printf(" ");
            }
            printf("%d",counter[j]);
            for(int t = 0; t < k; t++)
            {
                printf(" ");
            }
        }
        k *= 2;
        l = get_parent(l);
        r = get_parent(r);
        printf("\n");
    }
#endif

#ifdef USE_RANDKEY
    return (b + 1) * counter_num / 8
           + r * key_len * sizeof(int8_t)
           + sizeof(*this);
#else
    return (b + 1) * counter_num / 8
           + r * sizeof(hash_t)
           + sizeof(*this);
#endif
}

} // namespace OmniSketch