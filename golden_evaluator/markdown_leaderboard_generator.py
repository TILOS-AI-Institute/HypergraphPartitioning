#!/usr/bin/python3
import argparse
import os
import os.path
import pickle

from utils import Evaluator

algorithms = ["hMetis", "KaHyPar", "KaHyPar_new", "hMetis_SpecPart", "KaHyPar_SpecPart", "MtKaHyPar"]

parser = argparse.ArgumentParser()
parser.add_argument("result_dir", type=str)
parser.add_argument("instance_dir", type=str)

args = parser.parse_args()

best_results = { }
ubfactors = []
instances = []
for algorithm in algorithms:
  if os.path.isdir(args.result_dir + "/" + algorithm):
    for ubfactor in os.listdir(args.result_dir + "/" + algorithm):
      epsilon = int(ubfactor.split('_')[1])
      if epsilon not in ubfactors:
        ubfactors.append(epsilon)
      best_results_algo = { }
      cached_result_file = args.result_dir + "/" + algorithm + "/" + ubfactor + "/cached_results.json"

      # Check if there exists a cached result file in folder
      if os.path.exists(cached_result_file):
        with open(cached_result_file, "rb") as f:
          best_results_algo = pickle.load(f)
      # otherwise evaluate the cutsize of each solution in folder
      else:
        for solution in os.listdir(args.result_dir + "/" + algorithm + "/" + ubfactor):
          instance = solution.split('.')[0]
          instance_file = args.instance_dir + "/" + instance + ".hgr"
          solution_file = args.result_dir + "/" + algorithm + "/" + ubfactor + "/" + solution
          num_cut, blocks_balance, num_vertices, num_hyperedges = Evaluator(instance_file, solution_file, 2, epsilon)
          print("Processing algorithm " + algorithm + " and instance " + instance + " -> Cut: " + str(num_cut) + " - UBFactor " + str(epsilon))
          if blocks_balance[0] <= 0.5 + float(epsilon)/100.0 and blocks_balance[1] <= 0.5 + float(epsilon) / 100.0:
            if not instance in best_results_algo.keys():
              best_results_algo[instance] = {"cut": num_cut, "num_vertices": num_vertices, "num_hyperedges": num_hyperedges}
            elif num_cut < best_results_algo[instance]["cut"]:
              best_results_algo[instance]["cut"] = num_cut
        # Cache Results
        with open(cached_result_file, "wb") as f:
          pickle.dump(best_results_algo, f)

      # Compare results of algorithm with the previously best solution for each instance
      for instance in best_results_algo.keys():
        if not instance in instances:
          instances.append(instance)
        if not instance in best_results.keys():
          best_results[instance] = { "num_vertices": best_results_algo[instance]["num_vertices"],
                                     "num_hyperedges": best_results_algo[instance]["num_hyperedges"] }
        if not str(epsilon) in best_results[instance].keys():
          best_results[instance][str(epsilon)] = { "best_cut": pow(2,31), "best_algo": [] }

        if best_results_algo[instance]["cut"] < best_results[instance][str(epsilon)]["best_cut"]:
          best_results[instance][str(epsilon)]["best_cut"] = best_results_algo[instance]["cut"]
          best_results[instance][str(epsilon)]["best_algo"] = [algorithm]
        elif best_results_algo[instance]["cut"] == best_results[instance][str(epsilon)]["best_cut"]:
          best_results[instance][str(epsilon)]["best_algo"].append(algorithm)

# Print leaderboard table
ubfactors.sort()
instances.sort(key = lambda x: best_results[x]["num_hyperedges"])
columns = 3 + len(ubfactors)
print( "| Instance | Statistics | | Cutsize |" + (" |" * (columns - 4)))
print("|------------|" + (":------------:|" * (columns-1)))
print("|  | # Vertices | # Hyperedges | " + " | ".join(map(lambda x: "ε = "+str(x)+"%", ubfactors)) + " | ")

for instance in instances:
  print("| " + instance + " | " +
        str(best_results[instance]["num_vertices"]) + " | " +
        str(best_results[instance]["num_hyperedges"]) + " | ", end = "")
  for ubfactor in ubfactors:
    if str(ubfactor) in best_results[instance].keys():
      print(str(best_results[instance][str(ubfactor)]["best_cut"]), end = "")
      for algo in best_results[instance][str(ubfactor)]["best_algo"]:
        print("<br><sub>" + algo + "</sub>", end = "")
      print(" | ", end = "")

    else:
      print(" | ", end = "")
  print()
