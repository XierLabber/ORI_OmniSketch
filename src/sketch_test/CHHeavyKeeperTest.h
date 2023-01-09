/**
 * @file CHHeavyKeeperTest.h
 * @author XierLabber (you@domain.com)
 * @brief Test CH-optimized Heavy Keeper
 *
 * @copyright Copyright (c) 2022
 *
 */
#include <common/test.h>
#include <sketch/CHHeavyKeeper.h>

#define CHHK_PARA_PATH "HK.para"
#define CHHK_TEST_PATH "HK.test"
#define CHHK_DATA_PATH "HK.data"
#define CHHK_CH_PATH "HK.ch"

namespace OmniSketch::Test {

/**
 * @brief Testing class for Heavy Keeper
 *
 */
template <int32_t key_len, int32_t no_layer, typename T,
          typename hash_t = Hash::AwareHash>
class CHHeavyKeeperTest : public TestBase<key_len, T> {
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
  CHHeavyKeeperTest(const std::string_view config_file)
      : TestBase<key_len, T>("Heavy Keeper with CH", config_file, CHHK_TEST_PATH) {
  }

  /**
   * @brief Test CH-optimized Heavy Keeper
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
void CHHeavyKeeperTest<key_len, no_layer, T, hash_t>::runTest() {
  /**
   * @brief shorthand for convenience
   *
   */
  using StreamData = Data::StreamData<key_len>;

  /// Part I.
  ///   Parse the config file
  double num_heavy_hitter;
  int32_t depth, width, num_threshold;
  double b, hash_table_alpha;
  double cnt_no_ratio;
  std::vector<size_t> width_cnt, no_hash;
  size_t cm_r, cm_w;

  std::string data_file; // data config
  toml::array arr;       // shortly we will convert it to format

  /// Step ii. Open the config file
  Util::ConfigParser parser(config_file);
  if (!parser.succeed()) {
    return;
  }
  /// Step iii. Set the working node of the parser.
  parser.setWorkingNode(
      CHHK_PARA_PATH); // do not forget to to enclose it with braces
  /// Step iv. Parse num_bits and num_hash
  if (!parser.parseConfig(depth, "depth"))
    return;
  if (!parser.parseConfig(width, "width"))
    return;
  if (!parser.parseConfig(num_threshold, "num_threshold"))
    return;
  if (!parser.parseConfig(b, "b"))
    return;
  if (!parser.parseConfig(hash_table_alpha, "hash_table_alpha"))
    return;
  /// Step v. Move to the data node
  parser.setWorkingNode(CHHK_DATA_PATH);
  /// Step vi. Parse data and format
  if (!parser.parseConfig(num_heavy_hitter, "threshold_heavy_hitter"))
    return;
  if (!parser.parseConfig(data_file, "data"))
    return;
  if (!parser.parseConfig(arr, "format"))
    return;
  /// Step vi. Move to the CH node
  parser.setWorkingNode(CHHK_CH_PATH);
  if (!parser.parseConfig(cnt_no_ratio, "cnt_no_ratio"))
    return;
  if (!parser.parseConfig(width_cnt, "width_cnt"))
    return;
  if (!parser.parseConfig(no_hash, "no_hash"))
    return;
  if (!parser.parseConfig(cm_r, "cm_r"))
    return;
  if (!parser.parseConfig(cm_w, "cm_w"))
    return;

  Data::DataFormat format(arr); // conver from toml::array to Data::DataFormat
  /// [Optional] User-defined rules
  ///
  /// Step vii. Parse Cnt Method.
  parser.setWorkingNode(CHHK_DATA_PATH);
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
      new Sketch::CHHeavyKeeper<key_len, no_layer, T, hash_t>(
          depth, width, num_threshold, b, hash_table_alpha, 
          cnt_no_ratio, width_cnt, no_hash, cm_r, cm_w));
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
        gnd_truth_heavy_hitters); // gnd_truth_heavy_hitter: >, yet HeavyKeeper: >=
  }
  ///        3. size
  this->testSize(ptr);
  ///        3. show metrics
  this->show();

  printf("\n  HK DEPTH: %d\n  HK WIDTH: %d\n", depth, width);
  printf("  WIDTH_CNT: [");
  for(int i = 0; i < width_cnt.size();i++)
  {
    printf("%ld", width_cnt[i]);
    if(i != width_cnt.size() - 1)
    {
      printf(", ");
    }
  }
  printf("]\n");
  printf("  CM: %ld * %ld\n", cm_r, cm_w);
  printf("  RATIO: %lf\n\n", cnt_no_ratio);
  printf("============================================\n");

  return;
}

} // namespace OmniSketch::Test

#undef CHHK_PARA_PATH
#undef CHHK_TEST_PATH
#undef CHHK_DATA_PATH
#undef CHHK_CH_PATH

// Driver instance:
//      AUTHOR: dromniscience
//      CONFIG: sketch_config.toml  # with respect to the `src/` directory
//    TEMPLATE: <13, 2, int32_t, Hash::AwareHash>