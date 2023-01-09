/**
 * @file CHQueryingCountingBloomFilterTest.h
 * @author XierLabber (you@domain.com)
 * @brief Testing CHQCBF
 *
 * @copyright Copyright (c) 2022
 *
 */
#pragma once

#include <common/test.h>
#include <sketch/CHQueryingCountingBloomFilter.h>

#define CHQCBF_PARA_PATH "QCBF.para"
#define CHQCBF_TEST_PATH "QCBF.test"
#define CHQCBF_DATA_PATH "QCBF.data"
#define CHQCBF_CH_PATH "QCBF.ch"

namespace OmniSketch::Test {
/**
 * @brief Testing class for Bloom Filter
 *
 */
template <int32_t key_len, int32_t no_layer,  typename T, typename hash_t = Hash::AwareHash>
class CHQueryingCountingBloomFilterTest : public TestBase<key_len, T> {
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
  CHQueryingCountingBloomFilterTest(const std::string_view config_file)
      : TestBase<key_len, T>("CHQCBF", config_file, CHQCBF_TEST_PATH) {}

  /**
   * @brief Test CHCounting Bloom Filter
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

template <int32_t key_len, int32_t no_layer,  typename T, typename hash_t>
void CHQueryingCountingBloomFilterTest<key_len, no_layer, T, hash_t>::runTest() {
  /**
   * @brief shorthand for convenience
   *
   */
  using StreamData = Data::StreamData<key_len>;

  int32_t ncnt, nhash; // sketch config
  std::string data_file;     // data config
  toml::array arr;           // shortly we will convert it to format
  double cnt_no_ratio;
  std::vector<size_t> width_cnt, no_hash;
  size_t cm_r, cm_w;

  Util::ConfigParser parser(config_file);
  if (!parser.succeed()) {
    return;
  }
  parser.setWorkingNode(CHQCBF_PARA_PATH);
  if (!parser.parseConfig(ncnt, "num_cnt"))
    return;
  if (!parser.parseConfig(nhash, "num_hash"))
    return;

  parser.setWorkingNode(CHQCBF_DATA_PATH);
  if (!parser.parseConfig(data_file, "data"))
    return;
  if (!parser.parseConfig(arr, "format"))
    return;
  Data::DataFormat format(arr); // conver from toml::array to Data::DataFormat
  /// [Optional] User-defined rules
  ///
  /// Step vii. Parse Cnt Method.
  std::string method;
  Data::CntMethod cnt_method = Data::InLength;
  if (!parser.parseConfig(method, "cnt_method"))
    return;
  if (!method.compare("InPacket")) {
    cnt_method = Data::InPacket;
  }

  /// Step vi. Move to the CH node
  parser.setWorkingNode(CHQCBF_CH_PATH);
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


  std::unique_ptr<Sketch::SketchBase<key_len, T>> ptr(
      new Sketch::CHQueryingCountingBloomFilter<key_len, no_layer, T, hash_t>(ncnt, nhash, cnt_no_ratio,
                                                                   width_cnt, no_hash, cm_r,
                                                                   cm_w));

  StreamData data(data_file, format);
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
  this->testUpdate(ptr, data.begin(), data.end(),
                   cnt_method); // metrics of interest are in config file
  ///        2. query for all the flowkeys
  this->testQuery(ptr, gnd_truth); // metrics of interest are in config file
  ///        3. size
  this->testSize(ptr);
  ///        3. show metrics
  this->show();

  return;
}

} // namespace OmniSketch::Test

#undef CHQCBF_PARA_PATH
#undef CHQCBF_TEST_PATH
#undef CHQCBF_DATA_PATH

// Driver instance:
//      AUTHOR: dromniscience
//      CONFIG: sketch_config.toml  # with respect to the `src/` directory
//    TEMPLATE: <13, 2, int32_t, Hash::AwareHash>
