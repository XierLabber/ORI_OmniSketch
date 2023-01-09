/**
 * @file CHUnivMonTest.h
 * @author dromniscience (you@domain.com)
 * @brief Test CH-optimized UnivMon Sketch
 *
 * @copyright Copyright (c) 2022
 *
 */
#include <common/test.h>
#include <sketch/CHUnivMon.h>

#define CHUM_PARA_PATH "UM.para"
#define CHUM_TEST_PATH "UM.test"
#define CHUM_DATA_PATH "UM.data"
#define CHUM_CH_PATH "UM.ch"

namespace OmniSketch::Test {

/**
 * @brief Testing class for UnivMon Sketch
 *
 */
template <int32_t key_len, int32_t no_layer, typename T,
          typename hash_t = Hash::AwareHash>
class CHUnivMonTest : public TestBase<key_len, T> {
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
  CHUnivMonTest(const std::string_view config_file)
      : TestBase<key_len, T>("UnivMon with CH", config_file, CHUM_TEST_PATH) {
  }

  /**
   * @brief Test CH-optimized UnivMon Sketch
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
void CHUnivMonTest<key_len, no_layer, T, hash_t>::runTest() {
  /**
   * @brief shorthand for convenience
   *
   */
  using StreamData = Data::StreamData<key_len>;

  /// Part I.
  ///   Parse the config file
  int32_t depth, width;
  double cnt_no_ratio;
  std::vector<size_t> width_cnt, no_hash;

  std::string data_file; // data config
  toml::array arr;       // shortly we will convert it to format

  /// Step ii. Open the config file
  Util::ConfigParser parser(config_file);
  if (!parser.succeed()) {
    return;
  }
  /// Step iii. Set the working node of the parser.
  parser.setWorkingNode(
      CHUM_PARA_PATH); // do not forget to to enclose it with braces
  /// Step iv. Parse num_bits and num_hash
  if (!parser.parseConfig(depth, "depth"))
    return;
  if (!parser.parseConfig(width, "width"))
    return;
  /// Step v. Move to the data node
  parser.setWorkingNode(CHUM_DATA_PATH);
  /// Step vi. Parse data and format
  if (!parser.parseConfig(data_file, "data"))
    return;
  if (!parser.parseConfig(arr, "format"))
    return;
  /// Step vi. Move to the CH node
  parser.setWorkingNode(CHUM_CH_PATH);
  if (!parser.parseConfig(cnt_no_ratio, "cnt_no_ratio"))
    return;
  if (!parser.parseConfig(width_cnt, "width_cnt"))
    return;
  if (!parser.parseConfig(no_hash, "no_hash"))
    return;

  Data::DataFormat format(arr); // conver from toml::array to Data::DataFormat
  /// [Optional] User-defined rules
  ///
  /// Step vii. Parse Cnt Method.
  std::string method;
  Data::CntMethod cnt_method = Data::InLength;
  parser.setWorkingNode(CHUM_DATA_PATH);
  if (!parser.parseConfig(method, "cnt_method"))
    return;
  if (!method.compare("InPacket")) {
    cnt_method = Data::InPacket;
  }

  /// Step ii. Get ground truth
  ///
  ///       1. read data
  StreamData data(data_file, format); // specify both data file and data format
  if (!data.succeed())
    return;
  Data::GndTruth<key_len, T> gnd_truth;
  gnd_truth.getGroundTruth(data.begin(), data.end(), cnt_method);
  ///       2. [optional] show data info
  fmt::print("DataSet: {:d} records with {:d} keys ({})\n", data.size(),
             gnd_truth.size(), data_file);
  /// Step iii. Insert the samples and then look up all the flows
  ///
  ///        1. update records into the sketch

  int32_t logn = static_cast<int32_t>(std::log2(gnd_truth.size()));

  /// Part II.
  ///   Prepare sketch and data
  ///
  /// Step i. Initialize a sketch
  OmniSketch::Hash::AwareHash(1);
  std::unique_ptr<Sketch::SketchBase<key_len, T>> ptr(
      new Sketch::CHUnivMon<key_len, no_layer, T, hash_t>(
          depth, width, logn, cnt_no_ratio, width_cnt, no_hash));
  /// remember that the left ptr must point to the base class in order to call
  /// the methods in it

  this->testSize(ptr);
  this->show();

  this->testUpdate(ptr, data.begin(), data.end(),
                   cnt_method); // metrics of interest are in config file
  ///        2. query for all the flowkeys
  this->testQuery(ptr, gnd_truth); // metrics of interest are in config file
  ///        3. size
  this->testSize(ptr);
  ///        3. show metrics
  this->show();

  printf("\n  UM DEPTH: %d\n  UM WIDTH: %d\n", depth, width);
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
  printf("  RATIO: %lf\n\n", cnt_no_ratio);
  printf("============================================\n");

  return;
}

} // namespace OmniSketch::Test

#undef CHUM_PARA_PATH
#undef CHUM_TEST_PATH
#undef CHUM_DATA_PATH
#undef CHUM_CH_PATH

// Driver instance:
//      AUTHOR: dromniscience
//      CONFIG: sketch_config.toml  # with respect to the `src/` directory
//    TEMPLATE: <13, 2, int32_t, Hash::AwareHash>