/**
 * @file CHMVSketchTest.h
 * @author XierLabber (you@domain.com)
 * @brief Test CH-optimized MV Sketch Sketch
 *
 * @copyright Copyright (c) 2022
 *
 */
#include <common/test.h>
#include <sketch/CHMVSketch.h>

#define CHMV_PARA_PATH "MV.para"
#define CHMV_TEST_PATH "MV.test"
#define CHMV_DATA_PATH "MV.data"
#define CHMV_CH_PATH "MV.ch"

namespace OmniSketch::Test {

/**
 * @brief Testing class for MV Sketch Sketch
 *
 */
template <int32_t key_len, int32_t no_layer, typename T,
          typename hash_t = Hash::AwareHash>
class CHMVSketchTest : public TestBase<key_len, T> {
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
  CHMVSketchTest(const std::string_view config_file)
      : TestBase<key_len, T>("MV Sketch with CH", config_file, CHMV_TEST_PATH) {
  }

  /**
   * @brief Test CH-optimized MV Sketch Sketch
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
void CHMVSketchTest<key_len, no_layer, T, hash_t>::runTest() {
  /**
   * @brief shorthand for convenience
   *
   */
  using StreamData = Data::StreamData<key_len>;

  /// Part I.
  ///   Parse the config file
  int32_t depth, width;
  double V_cnt_no_ratio;
  std::vector<size_t> V_width_cnt, V_no_hash;
  double C_cnt_no_ratio;
  std::vector<size_t> C_width_cnt, C_no_hash;
  int32_t guess_negative_weight, guess_positive_weight;

  double num_heavy_hitter;
  std::string data_file; // data config
  toml::array arr;       // shortly we will convert it to format

  /// Step ii. Open the config file
  Util::ConfigParser parser(config_file);
  if (!parser.succeed()) {
    return;
  }
  /// Step iii. Set the working node of the parser.
  parser.setWorkingNode(
      CHMV_PARA_PATH); // do not forget to to enclose it with braces
  /// Step iv. Parse num_bits and num_hash
  if (!parser.parseConfig(depth, "depth"))
    return;
  if (!parser.parseConfig(width, "width"))
    return;
  /// Step v. Move to the data node
  parser.setWorkingNode(CHMV_DATA_PATH);
  /// Step vi. Parse data and format
  if (!parser.parseConfig(num_heavy_hitter, "threshold_heavy_hitter"))
    return;
  if (!parser.parseConfig(data_file, "data"))
    return;
  if (!parser.parseConfig(arr, "format"))
    return;
  /// Step vi. Move to the CH node
  parser.setWorkingNode(CHMV_CH_PATH);
  if (!parser.parseConfig(V_cnt_no_ratio, "V_cnt_no_ratio"))
    return;
  if (!parser.parseConfig(V_width_cnt, "V_width_cnt"))
    return;
  if (!parser.parseConfig(V_no_hash, "V_no_hash"))
    return;
  if (!parser.parseConfig(C_cnt_no_ratio, "C_cnt_no_ratio"))
    return;
  if (!parser.parseConfig(C_width_cnt, "C_width_cnt"))
    return;
  if (!parser.parseConfig(C_no_hash, "C_no_hash"))
    return;
  if (!parser.parseConfig(guess_negative_weight,  "guess_negative_weight"))
    return;
  if (!parser.parseConfig(guess_positive_weight, "guess_positive_weight"))
    return;

  Data::DataFormat format(arr); // conver from toml::array to Data::DataFormat
  /// [Optional] User-defined rules
  ///
  /// Step vii. Parse Cnt Method.
  parser.setWorkingNode(CHMV_DATA_PATH);
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
      new Sketch::CHMVSketch<key_len, no_layer, T, hash_t>(
          depth, width, V_cnt_no_ratio, V_width_cnt, V_no_hash,
          C_cnt_no_ratio, C_width_cnt, C_no_hash, guess_negative_weight,
          guess_positive_weight));
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
  this->testQuery(ptr, gnd_truth); // metrics of interest are in config file
  ///        3. query for heavy hitters
  if (hx_method == Data::TopK) {
    this->testHeavyHitter(
        ptr, gnd_truth_heavy_hitters.min(),
        gnd_truth_heavy_hitters); // metrics of interest are in config file
  } else {
    this->testHeavyHitter(
        ptr, std::floor(gnd_truth.totalValue() * num_heavy_hitter + 1),
        gnd_truth_heavy_hitters); // gnd_truth_heavy_hitter: >, yet HashPipe: >=
  }
  ///        3. size
  this->testSize(ptr);
  ///        3. show metrics
  this->show();

  printf("\n  CHMV DEPTH: %d\n  CHMV WIDTH: %d\n", depth, width);
  printf("  V_WIDTH_CNT: [");
  for(int i = 0; i < V_width_cnt.size();i++)
  {
    printf("%ld", V_width_cnt[i]);
    if(i != V_width_cnt.size() - 1)
    {
      printf(", ");
    }
  }
  printf("]\n");
  printf("  V_RATIO: %lf\n\n", V_cnt_no_ratio);
  printf("  C_WIDTH_CNT: [");
  for(int i = 0; i < C_width_cnt.size();i++)
  {
    printf("%ld", C_width_cnt[i]);
    if(i != C_width_cnt.size() - 1)
    {
      printf(", ");
    }
  }
  printf("]\n");
  printf("  C_RATIO: %lf\n\n", C_cnt_no_ratio);
  printf("============================================\n");

  return;
}

} // namespace OmniSketch::Test

#undef CHMV_PARA_PATH
#undef CHMV_TEST_PATH
#undef CHMV_DATA_PATH
#undef CHMV_CH_PATH

// Driver instance:
//      AUTHOR: XierLabber
//      CONFIG: sketch_config.toml  # with respect to the `src/` directory
//    TEMPLATE: <13, 2, int32_t, Hash::AwareHash>