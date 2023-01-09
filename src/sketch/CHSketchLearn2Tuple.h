/**
 * @file CHSketchLearn2Tuple.h
 * @author XierLabber<yangshibo@stu.pku.edu.cn>
 * @brief Sketch Learn
 *
 * @copyright Copyright (c) 2022
 *
 */
#pragma once

// #define MY_DEBUG

#include <common/hash.h>
#include <common/hierarchy.h>
#include <common/sketch.h>
#include <vector>
#include <algorithm>
#include <cmath>

namespace OmniSketch::Sketch {
/**
 * @brief Sketch Learn
 *
 * @tparam key_len  length of flowkey
 * @tparam hash_t   hashing class
 * @tparam T        type of the counter
 */
template <int32_t key_len, int32_t no_layer, typename T,
          typename hash_t>
class CHSketchLearn2Tuple : public SketchBase<key_len, T> {

private:

  /**
   * @brief sketch learn
   *
   * @param POSSIBLE_THRESHOLD  
   *         The threshold of hat_p, which is 0.99 as provided in the paper
   * 
   * @param STAR_THRESHOLD        
   *         If there are more than this many * in a regular expression, we 
   *         assume there is no large flow in the corresponding stack
   * 
   * @param MY_ERROR_THRESHOLD_SKETCH
   *         If the valuation is so many times higher than that of the value 
   *         in the smallest Sketch, that flow is considered likely to be a 
   *         false positive flow 
   * 
   * @param MY_ERROR_THRESHOLD_V0
   *         If the valuation is so many times higher than that of the total 
   *         value caught in V0 Sketch, that flow is considered likely to be 
   *         a false positive flow 
   * 
   * @param STEP
   *         The step size of the rate, which should decreases as theta 
   *         increases.
   *         Rate is used to determine whether to terminate the learning
   *         process.
   * 
   * @param START_THETA
   *         The initial value of theta
   * 
   * @param l
   *         The bit-level length of a flowkey
   *         We need l + 1 sketches for every bit of given flows and the 
   *         overall information of the given flows
   * 
   * @param r
   *         The depth of sketches
   * 
   * @param c 
   *         The width of sketches
   * 
   * @param V
   *         The sketches
   * 
   * @param p
   *         The array of mean values of Vk/V0
   * 
   * @param sigma
   *         The array of standard deviations of Vk/V0
   * 
   * @param updated
   *         To indicate whether some new information has been caught by
   *         the sketches.
   * 
   */
  
  const double POSSIBLE_THRESHOLD = 0.99;           
  const int32_t STAR_THRESHOLD = 11;                  
  const double MY_ERROR_THRESHOLD_SKETCH = 2.0;     
  const double MY_ERROR_THRESHOLD_V0 = 0.95;        
  const double STEP = 0.01;                        
  const double START_THETA = 0.5;                   

  static const int32_t l = 8 * key_len;                      
  int32_t r;                                                 
  int32_t c;                                                 

  hash_t* hash_function;

  std::vector<size_t> no_cnt;
  std::vector<size_t> width_cnt;
  std::vector<size_t> no_hash;
  CounterHierarchy<no_layer, T, hash_t> *ch;
  T ***V;

  double* p;
  double* sigma;
  bool updated;
  
  /**
   * @brief The information learned by CHSketchLearn2Tuple, which provides
   *        the estimated value and the probability array of the 
   *        specific flow
   *
   */
  class ans_t
  {
  public:
    char bit_flow[l + 2];        // the string describes the string, which consists of '0' and '1'
    char flow[(l + 7) / 8 + 1];  // the FlowKey represented by char array
    uint32_t size;                  // the size
    double prob_vector[l + 2];    // the probability array
    ans_t(char* bbit_flow, char* fflow, uint32_t ssize = 0, double* pprob_vector = NULL)
    {
        strcpy((char *)bit_flow, (char *)bbit_flow);
        for (int i = 0; i < key_len; i++)
        {
            flow[i] = fflow[i];
        }
        if (pprob_vector != NULL)
        {
            for (int i = 1; i <= l; i++)
            {
                prob_vector[i] = pprob_vector[i];
            }
        }
        size = ssize;
    }
  };

  /**
   * @brief Two ways to represent a specific flow
   *
   */
  class two_types_of_flow
  {
  public:
    char bit_flow[l + 2];       // the string describes the string, which consists of '0' and '1'
    char flow[(l + 7) / 8 + 1]; // the FlowKey represented by char array
    two_types_of_flow(char* bbit_flow, char* fflow)
    {
        strcpy(bit_flow, bbit_flow);
        for (int32_t i = 0; i < (l + 7) / 8; i++)
        {
            flow[i] = fflow[i];
        }
        flow[(l + 7) / 8] = '\0';
    }
  };

  char* current_string;
  int32_t num_of_star;
  std::vector<two_types_of_flow> possible_flows;
  std::vector<ans_t> large_flows;
  std::vector<ans_t> extracted_large_flows;
  std::vector<ans_t> flows_to_remove;

  int32_t get_bit(char* a, int32_t pos);
  void set_bit(char* a, int32_t pos, int32_t v);
  double normalCFD(double value);
  bool my_cmp(char* s1, char* s2);

  int32_t chIdx(int32_t k, int32_t i, int32_t j) const;

  void chUpdateCnt(int32_t k, int32_t i, int32_t j, T val);

  T chGetCnt(int32_t k, int32_t i, int32_t j) const;

  /**
   * @brief Find the flows which has a very high probability
   *        to be large flows
   *
   */
  void find_possible_flows(int32_t i, int32_t j, 
                           int32_t k, char* my_T);
  /**
   * @brief Assume that a large flow is hash into V[k][i][j],
   *        this function calculate the possibility of that, 
   *        the k-th bit of the flow is 1
   *
   */
  double cal_hat_p(double theta, int32_t i, int32_t j,
                     T*** V, double* p, 
                     double* sigma, int32_t k);
  /**
   * @brief Filter out fake streams, using the calculated prob_vector.
   *        This function should be implemented by the network 
   *        operators, so it is not implemented here.
   *
   */
  void large_flow_filter();
  /**
   * @brief Calculate p and sigma
   *
   */
  void Sketch2N_p_sigma();
  /**
   * @brief Ectract large flows in the stack (i,j)
   *        with the threshold theta
   *
   */
  void ExtractLargeFlows(double theta, int32_t i, int32_t j,
                         T*** V, double* p, double* sigma);
  /**
   * @brief Remove the extracted large flows from the sketches
   *        The large flows should be provided in flows_to_remove
   *
   */
  void RemoveFlows();
  /**
   * @brief Determine wether or not all the large flows is 
   *        extracted, if so, learning should be terminated
   *
   */
  bool Terminate(double theta);
  /**
   * @brief To learn from the sketches and extract all the large
   *        flows.
   *
   */
  void Sketch_Learning();

public:
  /**
   * @brief Construct by specifying depth and width
   *
   */
  CHSketchLearn2Tuple(int32_t depth_, int32_t width_, double cnt_no_ratio,
                const std::vector<size_t> &width_cnt,
                const std::vector<size_t> &no_hash);
  /**
   * @brief Release the pointer
   *
   */
  ~CHSketchLearn2Tuple();
  /**
   * @brief Update a flowkey with certain value
   *
   */
  void update(const FlowKey<key_len> &flowkey, T val) override;
  /**
   * @brief Get Heavy Hitter
   *
   */
  Data::Estimation<key_len, T> getHeavyHitter(double threshold) const override;
  /**
   * @brief Query a flowkey
   *
   */
  T query(const FlowKey<key_len> &flowkey) const override;
  /**
   * @brief Get the size of the sketch
   *
   */
  size_t size() const;

  void clear();
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
int32_t CHSketchLearn2Tuple<key_len, no_layer, T, hash_t>::chIdx(int32_t k, int32_t i, int32_t j) const{
    return k * c * r + i * c + j - 1;
}

template <int32_t key_len, int32_t no_layer, typename T,
          typename hash_t>
void CHSketchLearn2Tuple<key_len, no_layer, T, hash_t>::chUpdateCnt(int32_t k, int32_t i, int32_t j, T val){
    ch->updateCnt(chIdx(k, i, j), val);
}

template <int32_t key_len, int32_t no_layer, typename T,
          typename hash_t>
T CHSketchLearn2Tuple<key_len, no_layer, T, hash_t>::chGetCnt(int32_t k, int32_t i, int32_t j) const{
    return ch->getCnt(chIdx(k, i, j));
}

template <int32_t key_len, int32_t no_layer, typename T, typename hash_t>
int32_t CHSketchLearn2Tuple<key_len, no_layer, T, hash_t>::get_bit(char* a, int32_t pos){
      int32_t byte = pos / 8;
      int32_t bit = pos % 8;
      if ((a[byte] & (1 << (bit))) == 0)
      {
          return 0;
      }
      else
      {
          return 1;
      }
    }

template <int32_t key_len, int32_t no_layer, typename T, typename hash_t>
void CHSketchLearn2Tuple<key_len, no_layer, T, hash_t>::set_bit(char* a, int32_t pos, int32_t v){
    int32_t byte = pos / 8;
    int32_t bit = pos % 8;
    if (v == 1)
    {
        a[byte] = a[byte] | (1 << (bit));
    }
    else
    {
        a[byte] = a[byte] & ~(1 << (bit));
    }

}

template <int32_t key_len, int32_t no_layer, typename T, typename hash_t>
void CHSketchLearn2Tuple<key_len, no_layer, T, hash_t>::Sketch2N_p_sigma(){
    double sum;
    double square_sum;
    double tmp_r;
    int32_t total_times = r * c;
    for (size_t k = 0; k <= key_len * 8; k++)
    {
        sum = 0;
        square_sum = 0;
        for (size_t i = 0; i < r; i++)
        {
            for (size_t j = 1; j <= c; j++)
            {
                if(V[0][i][j] != 0)
                {
                    tmp_r = (double)(V[k][i][j]) / (double)(V[0][i][j]);
                    sum += tmp_r;
                    square_sum += tmp_r * tmp_r;
                }
                else
                {
                    total_times--;
                }
            }
        }
        p[k] = (double)sum / (double)total_times;
        double sigma2 = square_sum / (double)total_times- p[k] * p[k];
        sigma[k] = (sigma2 >= 0)? sqrt(sigma2) : 0;
    }
}

template <int32_t key_len, int32_t no_layer, typename T, typename hash_t>
double CHSketchLearn2Tuple<key_len, no_layer, T, hash_t>::normalCFD(double value){
    return 0.5 * erfc(-value / sqrt(2));
}

template <int32_t key_len, int32_t no_layer, typename T, typename hash_t>
void CHSketchLearn2Tuple<key_len, no_layer, T, hash_t>::find_possible_flows
  (int32_t i, int32_t j, int32_t k, char* candidate_string){
    if (k == l + 1)
    {
        char ans[(l + 7) / 8 + 1];

        for (int32_t kk = 1; kk <= l; kk++)
        {
            set_bit((char*)ans, kk - 1, current_string[kk] == '1' ? 1 : 0);
        }

        ans[(l + 7) / 8] = '\0';
        if ((hash_function[i]((FlowKey<key_len>)((int8_t *)ans)) % c + 1) == j)
        {
            int32_t flag = 0;
            possible_flows.push_back(two_types_of_flow(current_string, ans));
        }
        return;
    }
    else
    {
        if (candidate_string[k] != '*')
        {
            current_string[k] = candidate_string[k];
            find_possible_flows(i, j, k + 1, candidate_string);
        }
        else
        {
            current_string[k] = '0';
            find_possible_flows(i, j, k + 1, candidate_string);
            current_string[k] = '1';
            find_possible_flows(i, j, k + 1, candidate_string);
        }
    }
    return;
  }

template <int32_t key_len, int32_t no_layer, typename T, typename hash_t>
void CHSketchLearn2Tuple<key_len, no_layer, T, hash_t>::ExtractLargeFlows
  (double theta, int32_t i, int32_t j,T*** V, double* p, double* sigma){
    
    extracted_large_flows.clear();
    
    // 第一步，计算每个bit的概率估值
    double hat_p[l + 1];
    for (int32_t k = 1; k <= l; k++)
    {
        hat_p[k] = cal_hat_p(theta, i, j, V, p, sigma, k);
    }
    
    //  第二步，找到所有候选的大流，存在possible_flows里面
    num_of_star = 0;
    char candidate_string[l + 2];
    for (int32_t k = 1; k <= l; k++)
    {
        if (hat_p[k] > POSSIBLE_THRESHOLD)
        {
            candidate_string[k] = '1';
        }
        else if (1 - hat_p[k] > POSSIBLE_THRESHOLD)
        {
            candidate_string[k] = '0';
        }
        else
        {
            candidate_string[k] = '*';
            num_of_star++;
        }
    }
    candidate_string[l + 1] = '\0';
    candidate_string[0] = '#';
    if (num_of_star > STAR_THRESHOLD)
    {
        return;
    }
    current_string[l + 1] = '\0';
    current_string[0] = '#';
    possible_flows.clear();
    find_possible_flows(i, j, 1, candidate_string);
    
    //  第三步，估计大流的频率和可能性向量
    double estimated_frequency[l + 1];
    double estimated_p[l + 1];
    for (auto item = possible_flows.begin();
        item != possible_flows.end(); item++)
    {
        int32_t min_sketch = 0xfffffff;
        for (int32_t k = 1; k <= l; k++)
        {
            if (item->bit_flow[k] == '1')
            {
                min_sketch = (V[k][i][j] < min_sketch)? V[k][i][j] : min_sketch;
                double rate = (double)V[k][i][j] / V[0][i][j];
                estimated_frequency[k] = ((rate - p[k]) / (1 - p[k])) * V[0][i][j];
                estimated_p[k] = hat_p[k];
            }
            else
            {
                min_sketch = (V[0][i][j] - V[k][i][j] < min_sketch)? V[0][i][j] - V[k][i][j] : min_sketch;
                double rate = (double)V[k][i][j] / V[0][i][j];
                estimated_frequency[k] = (1 - rate / p[k]) * V[0][i][j];

                estimated_p[k] = 1 - hat_p[k];
            }
        }
        std::sort(estimated_frequency + 1, estimated_frequency + 1 + l);
        double ans_estimated_frequency = estimated_frequency[l / 2];
        if(ans_estimated_frequency > min_sketch)
        {
            if(ans_estimated_frequency > MY_ERROR_THRESHOLD_SKETCH * min_sketch && 
               ans_estimated_frequency > MY_ERROR_THRESHOLD_V0 * V[0][i][j])
            {
                break;
            }
            ans_estimated_frequency = min_sketch;
        }
        extracted_large_flows.push_back(ans_t(item->bit_flow, item->flow, ans_estimated_frequency, estimated_p));
    }
    
    //  第四步，去sketch里查候选流的数据，删掉过小的
    for (auto item = extracted_large_flows.begin(); item != extracted_large_flows.end(); )
    {
        for (int32_t ii = 0; ii < r; ii++)
        {
            if (ii == i)
            {
                continue;
            }
            int32_t jj = hash_function[ii](FlowKey<key_len>((const int8_t *)(item->flow))) % c + 1;
            for (int32_t k = 1; k <= l; k++)
            {
                if (item->bit_flow[k] == '0' && V[0][ii][jj] - V[k][ii][jj] < item->size)
                {
                    item->size = V[0][ii][jj] - V[k][ii][jj];
                }
                else if (item->bit_flow[k] == '1' && V[k][ii][jj] < item->size)
                {
                    item->size = V[k][ii][jj];
                }
            }
        }
        if (item->size < theta * V[0][i][j])
        {
            item = extracted_large_flows.erase(item);
            if (item == extracted_large_flows.end())break;
        }
        else 
        {
            item++;
        }
    }
    return;
  }

template <int32_t key_len, int32_t no_layer, typename T, typename hash_t>
double CHSketchLearn2Tuple<key_len, no_layer, T, hash_t>::cal_hat_p
  (double theta, int32_t i, int32_t j, T*** V, 
  double* p, double* sigma, int32_t k){
    double rate = (double)V[k][i][j] / V[0][i][j];
    if (rate < theta)
    {
        return 0;
    }
    if (1 - rate < theta)
    {
        return 1;
    }
    double ans = 0;
    double prob_1 = (V[k][i][j] - theta * V[0][i][j]) /
        (V[0][i][j] - theta * V[0][i][j]);
    double prob_0 = (V[k][i][j]) /
        (V[0][i][j] - theta * V[0][i][j]);
    double normal_val1 = normalCFD((prob_1 - p[k]) / sigma[k]);
    double normal_val0 = normalCFD((prob_0 - p[k]) / sigma[k]);
    return normal_val1 * p[k] + (1 - normal_val0) * (1 - p[k]);
  }

template <int32_t key_len, int32_t no_layer, typename T, typename hash_t>
void CHSketchLearn2Tuple<key_len, no_layer, T, hash_t>::RemoveFlows(){
    std::vector<ans_t> FF = flows_to_remove;
    uint32_t tmp_hash[r + 1];
    for (int32_t it = 0; it < FF.size(); it++)
    {
        char* ans = FF[it].flow;

        for (size_t i = 0; i < r; i++)
        {
            tmp_hash[i] = hash_function[i]((FlowKey<key_len>)((int8_t*)ans)) % c + 1;
            V[0][i][tmp_hash[i]] -= FF[it].size;
        }
        for (size_t k = 1; k <= key_len * 8; k++)
        {
            if (0 != get_bit((char*)ans, k - 1))
            {
                for (size_t i = 0; i < r; i++)
                {
                    V[k][i][tmp_hash[i]] -= FF[it].size;
                }
            }
        }
    }
}

template <int32_t key_len, int32_t no_layer, typename T, typename hash_t>
bool CHSketchLearn2Tuple<key_len, no_layer, T, hash_t>::Terminate(double theta){
  
    double RATE1 = 0.6826 + STEP * log2(theta);
    double RATE2 = 0.9544 + STEP * log2(theta);
    double RATE3 = 0.9973 + STEP * log2(theta);

    for (size_t k = 1; k <= key_len * 8; k++)
    {
        T sumV = 0;
        T sumV0 = 0;
        size_t sigma_num1 = 0, sigma_num2 = 0, sigma_num3 = 0;
        for (int32_t i = 0; i < r; i++)
        {
            for (int32_t j = 1; j <= c; j++)
            {
                sumV += V[k][i][j];
                sumV0 += V[0][i][j];
                double rate = (double)V[k][i][j] / V[0][i][j];
                if(sigma[k] != 0)
                {
                    if (rate <= p[k] + 3.0 * sigma[k] && rate >= p[k] - 3.0 * sigma[k])
                        sigma_num3++;
                    if (rate <= p[k] + 2.0 * sigma[k] && rate >= p[k] - 2.0 * sigma[k])
                        sigma_num2++;
                    if (rate <= p[k] + 1.0 * sigma[k] && rate >= p[k] - 1.0 * sigma[k])
                        sigma_num1++;
                }
                else
                {
                    sigma_num1++;
                    sigma_num2++;
                    sigma_num3++;
                }
            }
        }
        /*
        if(sigma_num1 == 0 && sigma_num2 == 0 && sigma_num3 == 0)
        {
            printf("I WONDER WHY PROGRAME RUNNING REACH HERE, k IS %ld, sumV IS %d, sumV0 IS %d, p[k] is %lf, sigma[k] is %lf\n", k, sumV, sumV0, p[k], sigma[k]);
        }
        */
        double rate1 = (double)sigma_num1 / (double)(r * c);
        double rate2 = (double)sigma_num2 / (double)(r * c);
        double rate3 = (double)sigma_num3 / (double)(r * c);
        
        if (rate1 < RATE1)
            return false;
        if (rate2 < RATE2)
            return false;
        if (rate3 < RATE3)
            return false;
    }
    return true;
}

template <int32_t key_len, int32_t no_layer, typename T, typename hash_t>
bool CHSketchLearn2Tuple<key_len, no_layer, T, hash_t>::my_cmp(char* s1, char* s2){
    for (size_t i = 0; i < key_len; i++)
    {
        if (s1[i] != s2[i])
        {
            return false;
        }
    }
    return true;
}

template <int32_t key_len, int32_t no_layer, typename T, typename hash_t>
CHSketchLearn2Tuple<key_len, no_layer, T, hash_t>::CHSketchLearn2Tuple(int32_t depth_, int32_t width_, 
             double cnt_no_ratio,
             const std::vector<size_t> &width_cnt_,
             const std::vector<size_t> &no_hash_)
    : r(depth_), c(Util::NextPrime(width_)), width_cnt(width_cnt_), no_hash(no_hash_){
    hash_function = new hash_t[r];

    // check ratio
    if (cnt_no_ratio <= 0.0 || cnt_no_ratio >= 1.0) {
        throw std::out_of_range("Out of Range: Ratio of #counters of adjacent "
                                "layers in CH should be in (0, 1), but got " +
                                std::to_string(cnt_no_ratio) + " instead.");
    }
    // prepare no_cnt
    no_cnt.push_back(static_cast<size_t>(this->l + 1) * this->r * this->c);
    for (int32_t i = 1; i < no_layer; ++i) {
        size_t last_layer = no_cnt.back();
        no_cnt.push_back(Util::NextPrime(std::ceil(last_layer * cnt_no_ratio)));
    }

    ch = new CounterHierarchy<no_layer, T, hash_t>(no_cnt, this->width_cnt,
                                                 this->no_hash);

    V = new T **[l + 1];
    for(int32_t i = 0; i < l + 1; i++)
    {
      V[i] = new T *[r];
      V[i][0] = new T[r * (c + 1)]();
      for(int32_t j = 1; j < r; j++)
      {
        V[i][j] = V[i][j-1] + (c + 1);
      }
    }
    p = new double[l + 1]();
    sigma = new double[l + 1]();
    current_string = new char[l + 2]();
    num_of_star = 0;
    updated = true;
}

template <int32_t key_len, int32_t no_layer, typename T, typename hash_t>
CHSketchLearn2Tuple<key_len, no_layer, T, hash_t>::~CHSketchLearn2Tuple(){
   delete[] hash_function;
   for(int32_t i = 0; i < l + 1; i++)
   {
     delete[] V[i][0];
     delete[] V[i];
   }
   delete[] V;
   delete[] ch;
   delete[] p;
   delete[] sigma;
   delete[] current_string;

   possible_flows.clear();
   large_flows.clear();
   extracted_large_flows.clear();
   flows_to_remove.clear();
}

template <int32_t key_len, int32_t no_layer, typename T, typename hash_t>
size_t CHSketchLearn2Tuple<key_len, no_layer, T, hash_t>::size() const{
   ch->print_rate("");
   return sizeof(*this)
          + r * sizeof(hash_t)
          + ch->size();
}

template <int32_t key_len, int32_t no_layer, typename T, typename hash_t>
void CHSketchLearn2Tuple<key_len, no_layer, T, hash_t>::update(const FlowKey<key_len> &flowkey, T val){

    #ifdef MY_DEBUG
    static int debug = 0;
    debug++;
    #endif

    uint32_t tmp_hash[r + 1];
    for (size_t i = 0; i < r; i++)
    {
        tmp_hash[i] = hash_function[i](flowkey) % c + 1;
        chUpdateCnt(0, i, tmp_hash[i], val);
        #ifdef MY_DEBUG
        V[0][i][tmp_hash[i]] += val;
        if(tmp_hash[0] == 1){
            static int tot = 0;
            tot++;
            if(tot % 10 == 0)
                printf("#: %d, %d, tot: %d, V: %d\n", debug, tmp_hash[0], tot, V[0][0][1]);
        }
        #endif
    }
    for (size_t k = 1; k <= key_len * 8; k++)
    {
        // 0 则对 V[k] 无影响
        if (0 != flowkey.getBit(k - 1))
        {
            for (size_t i = 0; i < r; i++)
            {
                chUpdateCnt(k, i, tmp_hash[i], val);
                #ifdef MY_DEBUG
                V[k][i][tmp_hash[i]] += val;
                #endif
            }
        }
    }
}

template <int32_t key_len, int32_t no_layer, typename T, typename hash_t>
void CHSketchLearn2Tuple<key_len, no_layer, T, hash_t>::Sketch_Learning(){
    printf("LEARNING START!\n");
    for(int i = 1; i < 10; i ++){
        printf("[%d] %d, V: %d\n", i, chGetCnt(0, 0, i), V[0][0][i]);
    }
    for(int k = 0; k <= l; k++){
        for(int i = 0; i < r; i++){
            for(int j = 1; j <= c; j++){
                #ifndef MY_DEBUG
                V[k][i][j] += chGetCnt(k, i, j);
                #else
                if(V[k][i][j] != chGetCnt(k, i, j)){
                    printf("WRONG ANSWER: ORI: %d, GOT: %d\n", V[k][i][j], chGetCnt(k, i, j));
                }
                V[k][i][j] = chGetCnt(k, i, j);
                #endif
            }
        }
    }

    double theta = START_THETA;
    int32_t nnnn = 0;
    large_flows.clear();
    while (1)
    {
        int32_t my_flow_num = 0;
        std::vector<ans_t> FF;
        for (int32_t i = 0; i < r; i++)
        {
            for (int32_t j = 1; j <= c; j++)
            {
                if (0 == V[0][i][j])
                {
                    continue;
                }
                ExtractLargeFlows(theta, i, j, V, p, sigma);
                std::vector<ans_t> temp_F = extracted_large_flows;
                if (!temp_F.empty())
                {
                    my_flow_num++;
                    for (auto it = temp_F.begin(); it < temp_F.end(); it++)
                    {
                        bool temp_Fin = false;
                        if (!FF.empty())
                        {
                            for (auto iter = FF.begin(); iter < FF.end(); iter++)
                            {
                                if (strcmp(iter->bit_flow, it->bit_flow) == 0)
                                {
                                    temp_Fin = true;
                                    break;
                                }
                            }
                        }
                        if (!temp_Fin)
                            FF.push_back(*it);
                    }
                }
            }
        }

        //本次循环找出大流时，剔除大流，重新计算期望、方差
        if (!FF.empty())
        {
            for (auto it = FF.begin(); it < FF.end(); it++)
            {
                bool FF_in = false;
                auto temp_pos = large_flows.end();
                if (!large_flows.empty())
                {
                    for (auto iter = large_flows.begin(); iter < large_flows.end(); iter++)
                    {
                        if (strcmp(iter->bit_flow, it->bit_flow) == 0)
                        {
                            FF_in = true;
                            temp_pos = iter;
                            break;
                        }
                    }
                }
                if (!FF_in)
                    large_flows.push_back(*it);
                else if (FF_in)
                {
                    temp_pos->size += it->size;
                }
            }
            flows_to_remove.clear();
            flows_to_remove = FF;
            RemoveFlows();
            Sketch2N_p_sigma();
        }
        // printf("%d loop is completed______________, theta = %lf\n\n", nnnn, theta);
        nnnn++;

        if (Terminate(theta))
            break;
        //没有找出大流，theta减半
        if (FF.empty())
            theta /= 2;
    }
    large_flow_filter();
    printf("LEARNING END!\n");
}

template <int32_t key_len, int32_t no_layer, typename T, typename hash_t>
void CHSketchLearn2Tuple<key_len, no_layer, T, hash_t>::large_flow_filter(){
  return;
}

template <int32_t key_len, int32_t no_layer, typename T, typename hash_t>
Data::Estimation<key_len, T> CHSketchLearn2Tuple<key_len, no_layer, T, hash_t>::getHeavyHitter(double threshold) const {
    if(updated || large_flows.size() == 0)
    {
      const_cast<CHSketchLearn2Tuple<key_len, no_layer, T, hash_t>*>(this)->Sketch_Learning();
      const_cast<CHSketchLearn2Tuple<key_len, no_layer, T, hash_t>*>(this)->updated = false;
    }
    Data::Estimation<key_len, T> heavy_hitters;
    for(auto it : large_flows)
    {
      if(it.size >= threshold)
      {
        heavy_hitters[FlowKey<key_len>((const int8_t *)it.flow)] = it.size;
      }
    }
    return heavy_hitters;
}

template <int32_t key_len, int32_t no_layer, typename T, typename hash_t>
T CHSketchLearn2Tuple<key_len, no_layer, T, hash_t>::query(const FlowKey<key_len> &flowkey) const{

    static bool cnt_distrib = true;
    if (cnt_distrib) {
        std::vector<double> distrib(32);
        for (int32_t i = 0; i < no_cnt[0]; ++i) {
        T val = ch->getOriginalCnt(i);
        for (int32_t k = 0; k < 32; ++k) {
            if (std::abs(val) >= (1 << k))
            distrib[k] += 1.0;
            else
            break;
        }
        }
        for (int32_t k = 0; k < 32; ++k) {
        std::cout << k << ": " << distrib[k] / no_cnt[0] << " \n";
        if (distrib[k] == 0.0) {
            std::cout << std::endl;
            break;
        }
        }
        cnt_distrib = false;
    }

    if(updated || large_flows.size() == 0)
    {
      const_cast<CHSketchLearn2Tuple<key_len, no_layer, T, hash_t>*>(this)->Sketch_Learning();
      const_cast<CHSketchLearn2Tuple<key_len, no_layer, T, hash_t>*>(this)->updated = false;
    }
    for(auto it : large_flows)
    {
      if(FlowKey<key_len>((const int8_t *)it.flow) == flowkey)
      {
        return it.size;
      }
    }
    T result = 0xfffffff;
    for (int32_t ii = 0; ii < r; ii++)
    {
        int32_t jj = hash_function[ii](flowkey) % c + 1;
        for (int32_t k = 1; k <= l; k++)
        {
            if ( flowkey.getBit(k - 1) == 0 && V[0][ii][jj] - V[k][ii][jj] < result)
            {
                result = V[0][ii][jj] - V[k][ii][jj];
            }
            else if (flowkey.getBit(k - 1) == 1 && V[k][ii][jj] < result)
            {
                result = V[k][ii][jj];
            }
        }
    }
    return (T)result;
}

template <int32_t key_len, int32_t no_layer, typename T, typename hash_t>
void CHSketchLearn2Tuple<key_len, no_layer, T, hash_t>::clear(){

  ch->clear();
  for(int32_t i = 0; i < l + 1; i++)
  {
    for(int32_t j = 0; j < r; j++)
    {
      for(int32_t k = 0; k < c + 1; k++)
      {
        V[i][j][k] = 0;
      }
    }
  }
  for(int32_t i = 0; i < l + 1; i++)
  {
    p[i] = sigma[i] = 0;
  }
  updated = false;
  possible_flows.clear();
  large_flows.clear();
  extracted_large_flows.clear();
  flows_to_remove.clear();
  num_of_star = 0;
}

}
