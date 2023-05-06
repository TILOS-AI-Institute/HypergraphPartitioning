#pragma once
#include <algorithm>
#include <cassert>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <istream>
#include <iterator>
#include <memory>
#include <sstream>
#include <vector>

namespace optimal_partitioner {

// operation for two vectors +, -, *,  ==, <
static std::vector<float> operator+(const std::vector<float> a,
                                    const std::vector<float> b) {
  assert(a.size() == b.size());
  std::vector<float> result;
  result.reserve(a.size());
  std::transform(a.begin(), a.end(), b.begin(), std::back_inserter(result),
                 std::plus<float>());
  return result;
}

static std::vector<float> operator-(const std::vector<float> a,
                                    const std::vector<float> b) {
  assert(a.size() == b.size());
  std::vector<float> result;
  result.reserve(a.size());
  std::transform(a.begin(), a.end(), b.begin(), std::back_inserter(result),
                 std::minus<float>());
  return result;
}

class Hypergraph {
public:
  Hypergraph() = default;
  Hypergraph(const int num_parts, const int ub_factor)
      : num_parts_(num_parts), ub_factor_(ub_factor) {}
  Hypergraph(const Hypergraph &) = default;
  Hypergraph(Hypergraph &&) = default;
  Hypergraph &operator=(const Hypergraph &) = default;
  Hypergraph &operator=(Hypergraph &&) = default;
  void ReadHypergraph(int num_parts_, float ub_factor,
                      std::string graph_file_name,
                      std::string fixed_file_name = "");
  void SolveCPLEX();
  void SolveOR();
  void WritePartitionToFile(std::string file_name);
  void Evaluator();

private:
  std::vector<float> GetTotalVertexWeights() const;
  std::vector<std::vector<float>> GetVertexBalance();
  std::vector<float> MultiplyFactor(const std::vector<float> a, float factor);
  int num_parts_;
  int vertex_dimensions_ = 1;
  int hyperedge_dimensions_ = 1;
  float ub_factor_;
  int num_vertices_;
  int num_hyperedges_;
  bool fixed_vertex_flag_;
  std::vector<int> partition_;
  std::vector<int> fixed_attr_;
  std::vector<std::vector<float>> hyperedge_weights_;
  std::vector<std::vector<float>> vertex_weights_;
  std::vector<int> eptr_;
  std::vector<int> eind_;
  std::vector<int> vptr_;
  std::vector<int> vind_;
};

using hypergraph_ptr = std::shared_ptr<Hypergraph>;
} // namespace optimal_partitioner