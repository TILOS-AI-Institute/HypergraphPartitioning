#include "IlpReadHypergraph.h"
#include "ilcplex/cplex.h"
#include "ilcplex/ilocplex.h"
#include <ortools/base/commandlineflags.h>
//#include <ortools/base/init_google.h>
#include <ortools/base/logging.h>
#include <ortools/linear_solver/linear_solver.h>
#include <ortools/linear_solver/linear_solver.pb.h>

#include "ortools/base/logging.h"
#include "ortools/sat/cp_model.h"
#include "ortools/sat/cp_model.pb.h"
#include "ortools/sat/cp_model_solver.h"

namespace optimal_partitioner {
using operations_research::MPConstraint;
using operations_research::MPObjective;
using operations_research::MPSolver;
using operations_research::MPVariable;
using operations_research::sat::BoolVar;
using operations_research::sat::CpModelBuilder;
using operations_research::sat::CpSolverResponse;
using operations_research::sat::CpSolverStatus;
using operations_research::sat::LinearExpr;
using operations_research::sat::SolutionBooleanValue;
using operations_research::sat::Solve;


void Hypergraph::ReadHypergraph(int num_parts, float ub_factor,
                                std::string hypergraph_file_name,
                                std::string fixed_file_name)

{
  int hyperedge_dimensions = 1;
  int vertex_dimensions = 1;
  ub_factor_ = ub_factor;
  num_parts_ = num_parts;
  std::vector<std::vector<int>> hyperedges;
  std::ifstream hypergraph_file_input(hypergraph_file_name);
  std::string cur_line;
  std::getline(hypergraph_file_input, cur_line);
  std::istringstream cur_line_buf(cur_line);
  std::vector<int> stats{std::istream_iterator<int>(cur_line_buf),
                         std::istream_iterator<int>()};
  num_hyperedges_ = stats[0];
  num_vertices_ = stats[1];
  bool hyperedge_weight_flag = false;
  bool vertex_weight_flag = false;
  if (stats.size() == 3) {
    if ((stats[2] % 10) == 1)
      hyperedge_weight_flag = true;
    if (stats[2] >= 10)
      vertex_weight_flag = true;
  }
  for (int i = 0; i < num_hyperedges_; i++) {
    std::getline(hypergraph_file_input, cur_line);
    if (hyperedge_weight_flag == true) {
      std::istringstream cur_line_buf(cur_line);
      std::vector<float> hvec{std::istream_iterator<float>(cur_line_buf),
                              std::istream_iterator<float>()};
      std::vector<float>::iterator breakpoint{hvec.begin() +
                                              hyperedge_dimensions};
      std::vector<float> hwts(hvec.begin(), breakpoint);
      std::vector<int> hyperedge(breakpoint, hvec.end());
      for (auto &value : hyperedge)
        value--;
      hyperedge_weights_.push_back(hwts);
      hyperedges.push_back(hyperedge);
    } else {
      std::istringstream cur_line_buf(cur_line);
      std::vector<int> hyperedge{std::istream_iterator<int>(cur_line_buf),
                                 std::istream_iterator<int>()};
      for (auto &value : hyperedge)
        value--;
      std::vector<float> hwts(hyperedge_dimensions, 1.0);
      hyperedge_weights_.push_back(hwts);
      hyperedges.push_back(hyperedge);
    }
  }

  for (int i = 0; i < num_vertices_; i++) {
    if (vertex_weight_flag == true) {
      std::getline(hypergraph_file_input, cur_line);
      std::istringstream cur_line_buf(cur_line);
      std::vector<float> vwts{std::istream_iterator<float>(cur_line_buf),
                              std::istream_iterator<float>()};
      vertex_weights_.push_back(vwts);
    } else {
      std::vector<float> vwts(vertex_dimensions, 1.0);
      vertex_weights_.push_back(vwts);
    }
  }

  if (fixed_file_name.size() > 0) {
    int part_id = -1;
    fixed_vertex_flag_ = true;
    std::ifstream fixed_file_input(fixed_file_name);
    for (int i = 0; i < num_vertices_; i++) {
      fixed_file_input >> part_id;
      fixed_attr_.push_back(part_id);
    }
    fixed_file_input.close();
  }
  num_vertices_ = vertex_weights_.size();
  num_hyperedges_ = hyperedge_weights_.size();

  eptr_.push_back(static_cast<int>(eind_.size()));
  for (auto hyperedge : hyperedges) {
    eind_.insert(eind_.end(), hyperedge.begin(), hyperedge.end());
    eptr_.push_back(static_cast<int>(eind_.size()));
  }
  // add vertex
  // create vertices from hyperedges
  std::vector<std::vector<int>> vertices(num_vertices_);
  for (int i = 0; i < num_hyperedges_; i++)
    for (auto v : hyperedges[i])
      vertices[v].push_back(i); // i is the hyperedge id
  vptr_.push_back(static_cast<int>(vind_.size()));
  for (auto &vertex : vertices) {
    vind_.insert(vind_.end(), vertex.begin(), vertex.end());
    vptr_.push_back(static_cast<int>(vind_.size()));
  }
  std::cout << "[info] vertices in hypergraph " << num_vertices_ << std::endl;
  std::cout << "[info] hyperedges in hypergraph " << num_hyperedges_
            << std::endl;
}

std::vector<float> Hypergraph::GetTotalVertexWeights() const {
  std::vector<float> total_weight(1, 0.0);
  for (auto &weight : vertex_weights_) {
    total_weight[0] = total_weight[0] + weight[0];
  }
  return total_weight;
}

std::vector<float> Hypergraph::MultiplyFactor(const std::vector<float> a,
                                              float factor) {
  std::vector<float> result = a;
  for (auto &value : result)
    value = value * factor;
  return result;
}

std::vector<std::vector<float>> Hypergraph::GetVertexBalance() {
  std::vector<float> vertex_balance = GetTotalVertexWeights();
  if (num_parts_ == 2) {
    vertex_balance =
        MultiplyFactor(vertex_balance, (ub_factor_ + 50.0) / 100.0);
  } else {
    vertex_balance = MultiplyFactor(vertex_balance, (1 + ub_factor_ * 0.01)/static_cast<float>(num_parts_));
  }
  std::cout << "Max balance constraint " << vertex_balance.front() << std::endl;
  return std::vector<std::vector<float>>(num_parts_, vertex_balance);
}

void Hypergraph::SolveOR() {
  auto start_time_stamp_global = std::chrono::high_resolution_clock::now();
  partition_.resize(num_vertices_);
  std::fill(partition_.begin(), partition_.end(), -1);
  std::cout << "[STATUS]Starting OR tools " << std::endl;
  auto max_block_balance = GetVertexBalance();
  std::vector<int> edge_mask;  // store the hyperedges being used.
  std::set<int> vertex_mask;   // store the vertices being used
  // define comp structure to compare hyperedge ( function: >)
  struct comp
  {
    // comparator function
    bool operator()(const std::pair<int, float>& l,
                    const std::pair<int, float>& r) const
    {
      if (l.second != r.second)
        return l.second > r.second;
      return l.first < r.first;
    }
  };

  // Here we blow up the hyperedge cost to avoid the unstability due to
  // conversion accuray We reserve the 6 digits after numeric point
  const float blow_factor = 1000000.0;
   // From here, we are going to create CP-SAT model for the reduced hypergraph
  std::vector<std::vector<int>> x(num_parts_, std::vector<int>(num_vertices_, 0));
  std::vector<std::vector<int>> y(num_parts_, std::vector<int>(num_hyperedges_, 0));
  // Build CP Model
  CpModelBuilder cp_model;
  // Variables
  // var_x[i][j] is an array of Boolean variables
  // var_x[i][j] is true if vertex i to partition j
  std::vector<std::vector<BoolVar>> var_x(num_parts_,
                                          std::vector<BoolVar>(num_vertices_));
  std::vector<std::vector<BoolVar>> var_y(num_parts_,
                                          std::vector<BoolVar>(num_hyperedges_));
  
  std::cout << "[STATUS] Setting variables " << std::endl;

  for (auto i = 0; i < num_parts_; i++) {
    // initialize var_x
    for (auto v = 0; v < num_vertices_; v++) {
      var_x[i][v] = cp_model.NewBoolVar();
    }
    // initialize var_y
    for (auto e = 0; e < num_hyperedges_; e++) {
      var_y[i][e] = cp_model.NewBoolVar();
    }
  }

  std::cout << "[STATUS] Setting balance constraints " << std::endl;

  // define constraints
  // balance constraint
  // check each dimension
  for (auto k = 0; k < vertex_dimensions_; k++) {
    // allowed balance for each dimension
    for (auto i = 0; i < num_parts_; i++) {
      LinearExpr balance_expr;
      for (int v = 0; v < num_vertices_; v++) {
        balance_expr
            += vertex_weights_[v][k] * var_x[i][v];
      }
      cp_model.AddLessOrEqual(balance_expr, max_block_balance[i][k]);
    }
  }

  std::cout << "[STATUS] Setting vertex constraints " << std::endl;

  // each vertex can only belong to one part
  for (auto v = 0; v < num_vertices_; v++) {
    std::vector<BoolVar> possible_partitions;
    for (auto i = 0; i < num_parts_; i++) {
      possible_partitions.push_back(var_x[i][v]);
    }
    cp_model.AddExactlyOne(possible_partitions);
  }

  std::cout << "[STATUS] Setting hyperedge constraints " << std::endl;

  // Hyperedge constraint
  for (auto e = 0; e < num_hyperedges_; e++) {
    const int start_idx = eptr_[e];
    const int end_idx = eptr_[e + 1];
    for (int j = start_idx; j < end_idx; j++) {
      const int v = eind_[j];
      for (int i = 0; i < num_parts_; i++) {
        cp_model.AddLessOrEqual(var_y[i][e], var_x[i][v]);
      }
    }
  }

  std::cout << "[STATUS] Setting objective function " << std::endl;

  // Objective (Maximize objective function -> Minimize cutsize)
  LinearExpr obj_expr;
  for (auto e = 0; e < num_hyperedges_; e++) {
    for (int i = 0; i < num_parts_; ++i) {
      // obj_expr += var_y[i][e] * cost_value ;
      obj_expr += var_y[i][e] * static_cast<int>(hyperedge_weights_[i].front() * blow_factor)
                  * num_vertices_ * num_parts_;
    }
    for (int v = 0; v < num_vertices_; v++) {
      for (int i = 0; i < num_parts_; i++) {
        obj_expr -= var_x[i][v] * (i * num_vertices_ * num_vertices_ + v);
      }
    }
  }
  cp_model.Maximize(obj_expr);

  std::cout << "[STATUS] Solving the ILP " << std::endl;

  // solve
  const CpSolverResponse response = Solve(cp_model.Build());

  std::cout << "[STATUS] Finished solving the ILP " << std::endl;

  // Print solution.
  if (response.status() == CpSolverStatus::OPTIMAL) {
    for (auto v = 0; v < num_vertices_; v++) {
      for (auto i = 0; i < num_parts_; i++) {
        if (SolutionBooleanValue(response, var_x[i][v])) {
          partition_[v] = i;
        }
      }
    }
  } else {
    for (auto &value : partition_)
      value = (value == -1) ? 0 : value;
  }
  auto end_time_stamp_global = std::chrono::high_resolution_clock::now();
  double total_global_time
      = std::chrono::duration_cast<std::chrono::nanoseconds>(
            end_time_stamp_global - start_time_stamp_global)
            .count();
  total_global_time *= 1e-9;
  std::cout << "Total time for ILP solve " << total_global_time << " seconds " << std::endl;
}

void Hypergraph::SolveCPLEX() {
  auto start_time_stamp_global = std::chrono::high_resolution_clock::now();
  auto max_block_balance = GetVertexBalance();
  std::vector<std::vector<float>> min_block_balance(num_parts_, GetTotalVertexWeights());
  if (num_parts_ == 2) {
    for (int i = 0; i < num_parts_; ++i) {
      min_block_balance[i] = GetTotalVertexWeights() - max_block_balance[i];
    }
  } else {
    for (int i = 0; i < num_parts_; ++i) {
       min_block_balance[i] = MultiplyFactor(
        min_block_balance[i],
        (1.0 - ub_factor_ * 0.01) / static_cast<float>(num_parts_));
    }
  }
  std::cout << "Min balance constraint " << min_block_balance[0].front() << std::endl;
  std::cout << "Total wts " << GetTotalVertexWeights().front() << std::endl;
  std::cout << "Num parts " << num_parts_ << std::endl;
  partition_.resize(num_vertices_);
  std::fill(partition_.begin(), partition_.end(), -1);
  // CPLEX-Based Implementation
  // define environment
  IloEnv myenv;
  IloModel mymodel(myenv);
  // Define constraints
  // For each vertex, define a variable x
  // For each hyperedge, define a variable y
  IloArray<IloNumVarArray> x(myenv, num_parts_);
  IloArray<IloNumVarArray> y(myenv, num_parts_);
  for (int i = 0; i < num_parts_; ++i) {
    x[i] = IloNumVarArray(myenv, num_vertices_, 0, 1, ILOINT);
    y[i] = IloNumVarArray(myenv, num_hyperedges_, 0, 1, ILOINT);
  }
  for (int i = 0; i < vertex_dimensions_; ++i) {
    for (int j = 0; j < num_parts_; ++j) {
      IloExpr balance_expr(myenv);
      for (int k = 0; k < num_vertices_; ++k) {
        balance_expr += vertex_weights_[k][i] * x[j][k];
      }
      // mymodel.add(IloRange(myenv, 0.0, balance_expr,
      // max_block_balance[j][i]));
      if (num_parts_ == 2) {
        mymodel.add(IloRange(myenv, min_block_balance[j][i], balance_expr,
                             max_block_balance[j][i]));
      } else {
        mymodel.add(
            IloRange(myenv, min_block_balance[j][i], balance_expr, max_block_balance[j][i]));
      }

      balance_expr.end();
    }
  }
  if (fixed_vertex_flag_ == true) {
    for (int i = 0; i < num_vertices_; ++i) {
      if (fixed_attr_[i] > -1) {
        mymodel.add(x[fixed_attr_[i]][i] == 1);
      }
    }
  }
  for (int i = 0; i < num_vertices_; ++i) {
    IloExpr vertex_expr(myenv);
    for (int j = 0; j < num_parts_; ++j) {
      vertex_expr += x[j][i];
    }
    mymodel.add(vertex_expr == 1);
    vertex_expr.end();
  }
  // Hyperedge constraint
  for (int i = 0; i < num_hyperedges_; ++i) {
    const int start_idx = eptr_[i];
    const int end_idx = eptr_[i + 1];
    for (int j = start_idx; j < end_idx; j++) {
      const int vertex_id = eind_[j];
      for (int k = 0; k < num_parts_; k++) {
        mymodel.add(y[k][i] <= x[k][vertex_id]);
      }
    }
  }
  // Maximize cutsize objective
  std::vector<float> e_wt_factors(hyperedge_dimensions_, 1.0);
  IloExpr obj_expr(myenv); // empty expression
  for (int i = 0; i < num_hyperedges_; ++i) {
    /*const float cost_value = std::inner_product(
        hyperedge_weights_[i].begin(),
        hyperedge_weights_[i].end(), e_wt_factors.begin(), 0.0);*/
    const float cost_value = hyperedge_weights_[i].front();
    for (int j = 0; j < num_parts_; ++j) {
      obj_expr += cost_value * y[j][i];
    }
  }
  mymodel.add(IloMaximize(myenv, obj_expr));
  obj_expr.end();
  // Model Solution
  IloCplex mycplex(myenv);
  mycplex.extract(mymodel);
  mycplex.setOut(myenv.getNullStream());
  mycplex.solve();
  IloBool feasible = mycplex.solve();
  if (feasible == IloTrue) {
    for (int i = 0; i < num_vertices_; ++i) {
      for (int j = 0; j < num_parts_; ++j) {
        if (mycplex.getValue(x[j][i]) == 1.00) {
          partition_[i] = j;
        }
      }
    }
    // some solution may invalid due to the limitation of ILP solver
    for (auto &value : partition_)
      value = (value == -1) ? 0 : value;
  } else {
    for (auto &value : partition_)
      value = (value == -1) ? 0 : value;
  }
  mycplex.clear();
  myenv.end();
  auto end_time_stamp_global = std::chrono::high_resolution_clock::now();
  double total_global_time
      = std::chrono::duration_cast<std::chrono::nanoseconds>(
            end_time_stamp_global - start_time_stamp_global)
            .count();
  total_global_time *= 1e-9;
  std::cout << "Total time for ILP solve " << total_global_time << " seconds " << std::endl;
}

void Hypergraph::Evaluator() {
  float cost = 0.0;
  for (int e = 0; e < num_hyperedges_; ++e) {
    for (int idx = eptr_[e] + 1; idx < eptr_[e + 1]; ++idx) {
      if (partition_[eind_[idx]] != partition_[eind_[idx - 1]]) {
        cost += hyperedge_weights_[e].front();
        break; // this net has been cut
      }
    } // finish hyperedge e
  }

  std::cout << "[info] cutcost recorded " << cost << std::endl;

  std::vector<float> blocks(num_parts_, 0.0);

  for (int v = 0; v < num_vertices_; ++v) {
    const int p = partition_[v];
    blocks[p] += vertex_weights_[v].front();
  }
  std::cout << "[info] blocks [";
  for (int i = 0; i < num_parts_; ++i) {
    std::cout << blocks[i] << " ";
  }
  std::cout << "]" << std::endl;
}

void Hypergraph::WritePartitionToFile(std::string solution_file) {
  std::ofstream solution_file_output;
  solution_file_output.open(solution_file);
  for (auto part_id : partition_)
    solution_file_output << part_id << std::endl;
  solution_file_output.close();
}
} // namespace optimal_partitioner