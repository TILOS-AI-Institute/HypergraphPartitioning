#include "IlpReadHypergraph.h"
#include <string>

int main(int argc, char **argv) {
  using namespace optimal_partitioner;
  std::string hypergraph_file = argv[1];
  int num_parts;
  int ub_factor;
  std::string fixed_file = "";
  if (argc > 4) {
    fixed_file = std::stoi(argv[2]);
    num_parts = std::stoi(argv[3]);
    ub_factor = std::stoi(argv[4]);
  } else {
    num_parts = std::stoi(argv[2]);
    ub_factor = std::stoi(argv[3]);
  }
  hypergraph_ptr hgraph = std::make_shared<Hypergraph>(num_parts, ub_factor);
  hgraph->ReadHypergraph(num_parts, ub_factor, hypergraph_file);
  hgraph->Solve();
  hgraph->Evaluator();
  hgraph->WritePartitionToFile(hypergraph_file + ".part." +
                               std::to_string(num_parts));
  return 0;
}