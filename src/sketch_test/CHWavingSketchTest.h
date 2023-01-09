/**
 * @file CHWavingSketchTest.h
 * @author XierLabber (you@domain.com)
 * @brief Test CH-optimized Waving Sketch
 *
 * @copyright Copyright (c) 2022
 *
 */
#include <common/test.h>
#include <sketch/CHWavingSketch.h>

#define CHWS_PARA_PATH "WS.para"
#define CHWS_TEST_PATH "WS.test"
#define CHWS_DATA_PATH "WS.data"
#define CHWS_CH_PATH "WS.ch"

namespace OmniSketch::Test {

/**
 * @brief Testing class for Waving Sketch
 *
 */
template <int32_t key_len, int32_t no_layer, typename T,
          typename hash_t = Hash::AwareHash>
class CHWavingSketchTest : public TestBase<key_len, T> {
  using TestBase<key_len, T>::config_file;

public:
  /**
   * @brief Constructor
   * @details Names from left to right are
   * - show name
   * - config file
   * - path to the node that contains metrics of interest (concatenated with
   * '.')
   */
  CHWavingSketchTest(const std::string_view config_file)
      : TestBase<key_len, T>("Waving Sketch with CH", config_file, CHWS_TEST_PATH) {
  }

  /**
   * @brief Test CH-optimized Waving Sketch
   * @details An overriden method
   */
  void runTest() override;
};

} // namespace OmniSketch::Test

//-----------------------------------------------------------------------------
//
///                        Implementation of templated methods
//
//-----------------------------------------------------------------------------

namespace OmniSketch::Test {

template <int32_t key_len, int32_t no_layer, typename T, typename hash_t>
void CHWavingSketchTest<key_len, no_layer, T, hash_t>::runTest() {
  /**
   * @brief shorthand for convenience
   *
   */
  using StreamData = Data::StreamData<key_len>;

  /// Part I.
  ///   Parse the config file
  double num_heavy_hitter;
  int32_t bucket_num, heavy_part_length;
  double counter_cnt_no_ratio;
  std::vector<size_t> counter_width_cnt, counter_no_hash;
  double heavy_cnt_no_ratio;
  std::vector<size_t> heavy_width_cnt, heavy_no_hash;
  size_t counter_cm_r, counter_cm_w, heavy_cm_r, heavy_cm_w;

  std::string data_file; // data config
  toml::array arr;       // shortly we will convert it to format

  /// Step ii. Open the config file
  Util::ConfigParser parser(config_file);
  if (!parser.succeed()) {
    return;
  }
  /// Step iii. Set the working node of the parser.
  parser.setWorkingNode(
      CHWS_PARA_PATH); // do not forget to to enclose it with braces
  /// Step iv. Parse num_bits and num_hash
  if (!parser.parseConfig(bucket_num, "bucket_num"))
    return;
  if (!parser.parseConfig(heavy_part_length, "heavy_part_length"))
    return;
  /// Step v. Move to the data node
  parser.setWorkingNode(CHWS_DATA_PATH);
  /// Step vi. Parse data and format
  if (!parser.parseConfig(num_heavy_hitter, "threshold_heavy_hitter"))
    return;
  if (!parser.parseConfig(data_file, "data"))
    return;
  if (!parser.parseConfig(arr, "format"))
    return;
  /// Step vi. Move to the CH node
  parser.setWorkingNode(CHWS_CH_PATH);
  if (!parser.parseConfig(counter_cnt_no_ratio, "counter_cnt_no_ratio"))
    return;
  if (!parser.parseConfig(counter_width_cnt, "counter_width_cnt"))
    return;
  if (!parser.parseConfig(counter_no_hash, "counter_no_hash"))
    return;
  if (!parser.parseConfig(counter_cm_r, "counter_cm_r"))
    return;
  if (!parser.parseConfig(counter_cm_w, "counter_cm_w"))
    return;
  if (!parser.parseConfig(heavy_cnt_no_ratio, "heavy_cnt_no_ratio"))
    return;
  if (!parser.parseConfig(heavy_width_cnt, "heavy_width_cnt"))
    return;
  if (!parser.parseConfig(heavy_no_hash, "heavy_no_hash"))
    return;
  if (!parser.parseConfig(heavy_cm_r, "heavy_cm_r"))
    return;
  if (!parser.parseConfig(heavy_cm_w, "heavy_cm_w"))
    return;

  Data::DataFormat format(arr); // conver from toml::array to Data::DataFormat
  /// [Optional] User-defined rules
  ///
  /// Step vii. Parse Cnt Method.
  parser.setWorkingNode(CHWS_DATA_PATH);
  std::string method;
  Data::HXMethod hx_method = Data::TopK;
  if (!parser.parseConfig(method, "hx_method"))
    return;
  if (!method.compare("Percentile")) {
    hx_method = Data::Percentile;
  }
  Data::CntMethod cnt_method = Data::InLength;
  if (!parser.parseConfig(method, "cnt_method"))
    return;
  if (!method.compare("InPacket")) {
    cnt_method = Data::InPacket;
  }

  /// Part II.
  ///   Prepare sketch and data
  ///
  /// Step i. Initialize a sketch
  std::unique_ptr<Sketch::SketchBase<key_len, T>> ptr(
      new Sketch::CHWavingSketch<key_len, no_layer, T, hash_t>(
          bucket_num, heavy_part_length, counter_cnt_no_ratio, 
          counter_width_cnt, counter_no_hash, heavy_cnt_no_ratio, 
          heavy_width_cnt, heavy_no_hash, counter_cm_r, counter_cm_w,
          heavy_cm_r, heavy_cm_w));
  /// remember that the left ptr must point to the base class in order to call
  /// the methods in it

  this->testSize(ptr);
  this->show();

  /// Step ii. Get ground truth
  ///
  ///       1. read data
  StreamData data(data_file, format); // specify both data file and data format
  if (!data.succeed())
    return;
  Data::GndTruth<key_len, T> gnd_truth, gnd_truth_heavy_hitters;
  gnd_truth.getGroundTruth(data.begin(), data.end(), cnt_method);
  gnd_truth_heavy_hitters.getHeavyHitter(gnd_truth, num_heavy_hitter,
                                         hx_method);
  ///       2. [optional] show data info
  fmt::print("DataSet: {:d} records with {:d} keys ({})\n", data.size(),
             gnd_truth.size(), data_file);
  /// Step iii. Insert the samples and then look up all the flows
  ///
  ///        1. update records into the sketch
  this->testUpdate(ptr, data.begin(), data.end(),
                   cnt_method); // metrics of interest are in config file
  ///        2. query for all the flowkeys
  if (hx_method == Data::TopK) {
    this->testHeavyHitter(
        ptr, gnd_truth_heavy_hitters.min(),
        gnd_truth_heavy_hitters); // metrics of interest are in config file
  } else {
    this->testHeavyHitter(
        ptr, std::floor(gnd_truth.totalValue() * num_heavy_hitter + 1),
        gnd_truth_heavy_hitters); // gnd_truth_heavy_hitter: >, yet WavingSketch: >=
  }
  ///        3. size
  this->testSize(ptr);
  ///        3. show metrics
  this->show();

  printf("\n  WS BUCKET NUM: %d\n  WS HEAVY PART LENGTH: %d\n", bucket_num, heavy_part_length);
  printf("  COUNTER_WIDTH_CNT: [");
  for(int i = 0; i < counter_width_cnt.size();i++)
  {
    printf("%ld", counter_width_cnt[i]);
    if(i != counter_width_cnt.size() - 1)
    {
      printf(", ");
    }
  }
  printf("]\n");
  printf("  COUNTER CM: %ld * %ld\n", counter_cm_r, counter_cm_w);
  printf("  COUNTER RATIO: %lf\n\n", counter_cnt_no_ratio);

  printf("  HEAVY_WIDTH_CNT: [");
  for(int i = 0; i < heavy_width_cnt.size();i++)
  {
    printf("%ld", heavy_width_cnt[i]);
    if(i != heavy_width_cnt.size() - 1)
    {
      printf(", ");
    }
  }
  printf("]\n");
  printf("  HEAVY CM: %ld * %ld\n", heavy_cm_r, heavy_cm_w);
  printf("  HEAVY RATIO: %lf\n\n", heavy_cnt_no_ratio);
  printf("============================================\n");

  return;
}

} // namespace OmniSketch::Test

#undef CHWS_PARA_PATH
#undef CHWS_TEST_PATH
#undef CHWS_DATA_PATH
#undef CHWS_CH_PATH

// Driver instance:
//      AUTHOR: dromniscience
//      CONFIG: sketch_config.toml  # with respect to the `src/` directory
//    TEMPLATE: <13, 2, int32_t, Hash::AwareHash>