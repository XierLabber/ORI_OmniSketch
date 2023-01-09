/**
 * @file CHFlowRadarTest.h
 * @author XierLabber
 * @brief Test CH Flow Radar
 *
 * @copyright Copyright (c) 2022
 *
 */
#pragma once

#include <common/test.h>
#include <common/hash.h>
#include <common/hierarchy.h>
#include <common/sketch.h>
#include <sketch/CHFlowRadar.h>

#define CHFR_PARA_PATH "FlowRadar.para"
#define CHFR_TEST_PATH "FlowRadar.test"
#define CHFR_DATA_PATH "FlowRadar.data"
#define CHFR_CH_PATH "FlowRadar.ch"

namespace OmniSketch::Test {

/**
 * @brief Testing class for CH Flow Radar
 *
 */
template <int32_t key_len, int32_t no_layer, typename T,
          typename hash_t>
class CHFlowRadarTest : public TestBase<key_len, T> {
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
  CHFlowRadarTest(const std::string_view config_file)
      : TestBase<key_len, T>("Flow Radar with CH", config_file, CHFR_TEST_PATH) {}

  /**
   * @brief Test CH Flow Radar
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

template <int32_t key_len, int32_t no_layer, typename T,
          typename hash_t>
void CHFlowRadarTest<key_len, no_layer, T, hash_t>::runTest(){
  // for convenience only
  using StreamData = Data::StreamData<key_len>;

  // parse config
  int32_t flow_filter_bit, flow_filter_hash, count_table_num,
      count_table_hash;  // sketch config
  double flow_cnt_no_ratio, packet_cnt_no_ratio;
  std::vector<size_t> flow_width_cnt, flow_no_hash, packet_width_cnt, packet_no_hash;
  std::string data_file; // data config
  toml::array arr;       // shortly we will convert it to format

  Util::ConfigParser parser(config_file);
  if (!parser.succeed()) {
    return;
  }

  parser.setWorkingNode(CHFR_PARA_PATH);
  if (!parser.parseConfig(flow_filter_bit, "flow_filter_bit"))
    return;
  if (!parser.parseConfig(flow_filter_hash, "flow_filter_hash"))
    return;
  if (!parser.parseConfig(count_table_num, "count_table_num"))
    return;
  if (!parser.parseConfig(count_table_hash, "count_table_hash"))
    return;

  // prepare data
  parser.setWorkingNode(CHFR_DATA_PATH);
  if (!parser.parseConfig(data_file, "data"))
    return;
  if (!parser.parseConfig(arr, "format"))
    return;

  parser.setWorkingNode(CHFR_CH_PATH);
  if (!parser.parseConfig(flow_cnt_no_ratio, "flow_cnt_no_ratio"))
    return;
  if (!parser.parseConfig(flow_width_cnt, "flow_width_cnt"))
    return;
  if (!parser.parseConfig(flow_no_hash, "flow_no_hash"))
    return;
  if (!parser.parseConfig(packet_cnt_no_ratio, "packet_cnt_no_ratio"))
    return;
  if (!parser.parseConfig(packet_width_cnt, "packet_width_cnt"))
    return;
  if (!parser.parseConfig(packet_no_hash, "packet_no_hash"))
    return;

  Data::DataFormat format(arr);
  StreamData data(data_file, format);
  if (!data.succeed())
    return;
  Data::GndTruth<key_len, T> gnd_truth;
  gnd_truth.getGroundTruth(data.begin(), data.end(), Data::InPacket);
  fmt::print("DataSet: {:d} records with {:d} keys ({})\n", data.size(),
             gnd_truth.size(), data_file);

  Data::GndTruth<key_len, T> heavy_part;
  heavy_part.getHeavyHitter(gnd_truth, gnd_truth.size() * 0.3, Data::TopK);
  printf("SIZE: %ld\n", heavy_part.size());

  std::unique_ptr<Sketch::SketchBase<key_len, T>> ptr(
      new Sketch::CHFlowRadar<key_len, no_layer, T, hash_t>(
          flow_filter_bit, flow_filter_hash, count_table_num,
          count_table_hash, flow_cnt_no_ratio, flow_width_cnt,
          flow_no_hash, packet_cnt_no_ratio, packet_width_cnt,
          packet_no_hash));

  this->testSize(ptr);
  this->show();

  this->testUpdate(ptr, data.begin(), data.end(), Data::InPacket);
  this->testDecode(ptr, heavy_part);
  this->testSize(ptr);
  // show
  this->show();

  printf("\n  FLOW FILTER BIT SIZE: %d\n  FLOW FILTER HASH NUM: %d\n  COUNT TABLE NUM: %d\n  COUNT HASH NUM: %d\n", flow_filter_bit, flow_filter_hash, count_table_num, count_table_hash);
  printf("  FLOW WIDTH_CNT: [");
  for(int i = 0; i < flow_width_cnt.size();i++)
  {
    printf("%ld", flow_width_cnt[i]);
    if(i != flow_width_cnt.size() - 1)
    {
      printf(", ");
    }
  }
  printf("]\n");
  printf("  FLOW CH RATIO: %lf\n", flow_cnt_no_ratio);
  printf("  PACKET WIDTH_CNT: [");
  for(int i = 0; i < packet_width_cnt.size();i++)
  {
    printf("%ld", packet_width_cnt[i]);
    if(i != packet_width_cnt.size() - 1)
    {
      printf(", ");
    }
  }
  printf("]\n");
  printf("  PACKET CH RATIO: %lf\n\n", packet_cnt_no_ratio);    
  printf("============================================\n");


  return;
}

} // namespace OmniSketch::Test

#undef CHFR_PARA_PATH
#undef CHFR_TEST_PATH
#undef CHFR_DATA_PATH
#undef CHFR_CH_PATH

// Driver instance:
//      AUTHOR: dromniscience
//      CONFIG: sketch_config.toml  # with respect to the `src/` directory
//    TEMPLATE: <13, 2, int32_t, Hash::AwareHash>