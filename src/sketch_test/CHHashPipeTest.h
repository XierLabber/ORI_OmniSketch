/**
 * @file CHHashPipeTest.h
 * @author XierLabber (lsmfttb@gmail.com)
 * @brief Test CHHashPipe
 *
 * @copyright Copyright (c) 2022
 *
 */
#pragma once

#include <common/test.h>
#include <sketch/CHHashPipe.h>

#define CHHP_PARA_PATH "HP.para"
#define CHHP_TEST_PATH "HP.test"
#define CHHP_DATA_PATH "HP.data"
#define CHHP_CH_PATH "HP.ch"

namespace OmniSketch::Test {

/**
 * @brief Testing class for CHHashPipe
 *
 */
template <int32_t key_len, int32_t no_layer, typename T,
          typename hash_t = Hash::AwareHash>
class CHHashPipeTest : public TestBase<key_len, T> {
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
  CHHashPipeTest(const std::string_view config_file)
      : TestBase<key_len, T>("HashPipe with CH", config_file, CHHP_TEST_PATH) {}

  /**
   * @brief Test CHHashPipe
   * @details An overriden method
   */
  void runTest() override;
};

} // namespace OmniSketch::Test

//-----------------------------------------------------------------------------
//
///                        Implementation of template methods
//
//-----------------------------------------------------------------------------

namespace OmniSketch::Test {

template <int32_t key_len, int32_t no_layer, typename T, typename hash_t>
void CHHashPipeTest<key_len, no_layer, T, hash_t>::runTest() {
  /**
   * @brief shorthand for convenience
   *
   */
  using StreamData = Data::StreamData<key_len>;

  /// Part I.
  ///   Parse the config file
  ///
  /// Step i.  First we list the variables to parse, namely:
  ///
  int32_t depth, ch_depth, width, chcm_r, chcm_c; // sketch config
  double cnt_no_ratio;
  double num_heavy_hitter;
  std::vector<size_t> width_cnt, no_hash;

  std::string data_file;       // data config
  toml::array arr;             // shortly we will convert it to format
  
  /// Step ii. Open the config file
  Util::ConfigParser parser(config_file);
  if (!parser.succeed()) {
    return;
  }
  /// Step iii. Set the working node of the parser. In this case, [CHHashPipe].
  parser.setWorkingNode(
      CHHP_PARA_PATH); // do not forget to to enclose it with braces
  /// Step iv. Parse num_hash, num_group and threshold
  if (!parser.parseConfig(depth, "depth"))
    return;
  if (!parser.parseConfig(width, "width"))
    return;
  /// Step v. To know about the data, we move to the [CHHashPipe.data] node.
  parser.setWorkingNode(CHHP_DATA_PATH);
  /// Step vi. Parse data and format
  if (!parser.parseConfig(num_heavy_hitter, "threshold_heavy_hitter"))
    return;
  if (!parser.parseConfig(data_file, "data"))
    return;
  if (!parser.parseConfig(arr, "format"))
    return;

  Data::DataFormat format(arr); // conver from toml::array to Data::DataFormat
  /// [Optional] User-defined rules
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

  /// Step vi. Move to the CH node
  parser.setWorkingNode(CHHP_CH_PATH);
  if (!parser.parseConfig(ch_depth, "ch_depth"))
    return;
  if (!parser.parseConfig(cnt_no_ratio, "cnt_no_ratio"))
    return;
  if (!parser.parseConfig(width_cnt, "width_cnt"))
    return;
  if (!parser.parseConfig(no_hash, "no_hash"))
    return;
  if (!parser.parseConfig(chcm_r, "chcm_r"))
    return;
  if (!parser.parseConfig(chcm_c, "chcm_c"))
    return;

  /// Part II.
  ///   Prepare sketch and data
  ///
  /// Step i. Initialize a sketch
  std::unique_ptr<Sketch::SketchBase<key_len, T>> ptr(
      new Sketch::CHHashPipe<key_len, no_layer, T, hash_t>(depth, width, cnt_no_ratio, width_cnt, no_hash, chcm_r, chcm_c, ch_depth));
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
        gnd_truth_heavy_hitters); // gnd_truth_heavy_hitter: >, yet HashPipe: >=
  }
  ///        4. size
  this->testSize(ptr);
  ///        5. show metrics
  this->show();

  printf("\n  DEPTH: %d\n  CH DEPTH: %d  WIDTH: %d\n", depth, (ch_depth == -1)? depth : ch_depth, width);
  printf("  WIHPH_CNT: [");
  for(int i = 0; i < width_cnt.size();i++)
  {
    printf("%ld", width_cnt[i]);
    if(i != width_cnt.size() - 1)
    {
      printf(", ");
    }
  }
  printf("]\n");
  printf("  CHCM: %d * %d\n", chcm_r, chcm_c);
  printf("  RATIO: %lf\n\n", cnt_no_ratio);
  printf("============================================\n");

  return;
}

} // namespace OmniSketch::Test

#undef CHHP_PARA_PATH
#undef CHHP_TEST_PATH
#undef CHHP_DATA_PATH
#undef CHHP_CH_PATH

// Driver instance:
//      AUTHOR: XierLabber
//      CONFIG: sketch_config.toml  # with respect to the `src/` directory
//    TEMPLATE: <13, 2, int32_t, Hash::AwareHash>