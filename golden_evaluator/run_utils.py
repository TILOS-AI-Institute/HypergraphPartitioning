import os 

import sys


from utils import Evaluator

benchmark_directory = "/home/fetzfs_projects/SpecPart/testcases"
solutions_directory = "/home/fetzfs_projects/SpecPart/partition_dump"




#Run IBM benchmarks

summary_report = "ibm_summary.csv"
line = "Benchmark, Numparts, Cutsize"
f = open(summary_report, "w")
f.write(line + "\n")
f.close()

ibm_benchmarks = ["ibm01",
                 "ibm02",
                 "ibm03",
                 "ibm04",
                 "ibm05",
                 "ibm06",
                 "ibm07",
                 "ibm08",
                 "ibm09",
                 "ibm10",
                 "ibm11",
                 "ibm12",
                 "ibm13",
                 "ibm14",
                 "ibm15",
                 "ibm16",
                 "ibm17",
                 "ibm18"]

num_parts_list = [2, 3, 4]
ub_factor = 2

for benchmark in ibm_benchmarks:
    print("[INFO] Evaluating ", benchmark)
    for num_parts in num_parts_list:
        hypergraph_file = benchmark_directory + "/ibm/" + benchmark + ".hgr"
        solutions_file = solutions_directory + "/" + benchmark + ".hgr.specpart.ubfactor." + str(ub_factor) + ".part." + str(num_parts)
        print("Hypergraph file: ", hypergraph_file)
        print("Solutions file: ", solutions_file)
        num_cut, blocks_balance, vtxs, hes = Evaluator(hypergraph_file, solutions_file, num_parts, ub_factor)
        line = benchmark + "," + str(num_parts) + "," + str(num_cut)
        f = open(summary_report, "a")
        f.write(line + "\n")
        f.close()


#Run IBM weighted benchmarks

summary_report = "ibm_w_summary.csv"
line = "Benchmark, Numparts, Cutsize"
f = open(summary_report, "w")
f.write(line + "\n")
f.close()

ibm_w_benchmarks = ["ibm01.weight",
                 "ibm02.weight",
                 "ibm03.weight",
                 "ibm04.weight",
                 "ibm05.weight",
                 "ibm06.weight",
                 "ibm07.weight",
                 "ibm08.weight",
                 "ibm09.weight",
                 "ibm10.weight",
                 "ibm11.weight",
                 "ibm12.weight",
                 "ibm13.weight",
                 "ibm14.weight",
                 "ibm15.weight",
                 "ibm16.weight",
                 "ibm17.weight",
                 "ibm18.weight"]

num_parts_list = [2, 3, 4]
ub_factor = 2

for benchmark in ibm_w_benchmarks:
    print("[INFO] Evaluating ", benchmark)
    for num_parts in num_parts_list:
        hypergraph_file = benchmark_directory + "/ibm_w/" + benchmark + ".hgr"
        solutions_file = solutions_directory + "/" + benchmark + ".hgr.specpart.ubfactor." + str(ub_factor) + ".part." + str(num_parts)
        print("Hypergraph file: ", hypergraph_file)
        print("Solutions file: ", solutions_file)
        num_cut, blocks_balance, vtxs, hes = Evaluator(hypergraph_file, solutions_file, num_parts, ub_factor)
        line = benchmark + "," + str(num_parts) + "," + str(num_cut)
        f = open(summary_report, "a")
        f.write(line + "\n")
        f.close()


#Run Titan benchmarks

summary_report = "titan_summary.csv"
line = "Benchmark, Numparts, Cutsize"
f = open(summary_report, "w")
f.write(line + "\n")
f.close()

titan_benchmarks = [
   "sparcT1_core",
   "neuron",
   "stereo_vision",
   "des90",
   "SLAM_spheric",
   "cholesky_mc",
   "segmentation",
   "bitonic_mesh",
   "dart",
   "openCV",
   "stap_qrd",
   "minres",
   "cholesky_bdti",
   "denoise",
   "sparcT2_core",
   "gsm_switch",
   "mes_noc",
   "LU230",
   "LU_Network",
   "sparcT1_chip2",
   "directrf",
   "bitcoin_miner"]

num_parts_list = [2, 3, 4]
ub_factor = 2

for benchmark in titan_benchmarks:
    print("[INFO] Evaluating ", benchmark)
    for num_parts in num_parts_list:
        hypergraph_file = benchmark_directory + "/titan/" + benchmark + ".hgr"
        solutions_file = solutions_directory + "/" + benchmark + ".hgr.specpart.ubfactor." + str(ub_factor) + ".part." + str(num_parts)
        print("Hypergraph file: ", hypergraph_file)
        print("Solutions file: ", solutions_file)
        num_cut, blocks_balance, vtxs, hes = Evaluator(hypergraph_file, solutions_file, num_parts, ub_factor)
        line = benchmark + "," + str(num_parts) + "," + str(num_cut)
        f = open(summary_report, "a")
        f.write(line + "\n")
        f.close()


