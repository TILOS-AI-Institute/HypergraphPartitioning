function IterativeMSTReconstruction(tree::SimpleWeightedGraph, edge_list::Vector{Vector{Int}})
    for edge in edge_list
        src_vertex = edge[1]
        dst_vertex = edge[2]

        tree.weights[src_vertex, dst_vertex] *= 1e07
    end
end

function FindFeasibleEdges(tree::AbstractGraph, tree_cut_token::TreeCutToken, hypergraph::Hypergraph, incidence_struct::Incidence, fixed_vtxs::Pindex, capacities::Vector{Int}, ncuts::Int)
	cutInfo = analyzeCutsOnTree(hypergraph, incidence_struct, hypergraph.vwts, fixed_vtxs, tree, tree_cut_token.hyperedges_flag)
	hyperedge_wts = hypergraph.hwts
	n = hypergraph.n
    vtxCuts = cutInfo.vtxCuts
	edgeCuts = cutInfo.edgeCuts
	pred = cutInfo.pred
	pindex = cutInfo.pindex
	p1 = pindex.p1
	p2 = pindex.p2
	forced_0 = cutInfo.forced_0
	forced_1 = cutInfo.forced_1
	forced_01 = cutInfo.forced_01
	FB0 = cutInfo.FB0
	FB1 = cutInfo.FB1
	edgeCuts0 = cutInfo.edgeCuts0
	edgeCuts1 = cutInfo.edgeCuts1
    tree_cut_token.pred = cutInfo.pred
    tree_cut_token.nforced_0 = sum(hyperedge_wts[forced_0])
    tree_cut_token.nforced_1 = sum(hyperedge_wts[forced_1])
    tree_cut_token.nforced_01 = sum(hyperedge_wts[forced_01])
    total_vwts = sum(hypergraph.vwts)
    feasible_cuts = Vector{Int}()
    feasible_edges = Vector{Vector{Int}}()

    for i in 1:hypergraph.n
		tree_cut_token.exc0[i] = edgeCuts[i] + tree_cut_token.nforced_0 - FB0[i] + edgeCuts1[i] + tree_cut_token.nforced_01
		tree_cut_token.exc1[i] = edgeCuts[i] + tree_cut_token.nforced_1 - FB1[i] + edgeCuts0[i] + tree_cut_token.nforced_01
	end

	tree_cut_token.exc0[1] = tree_cut_token.exc1[1] = hypergraph.e

	for i in 1:length(tree_cut_token.exc0)
		tree_cut_token.cutCost0[i] = tree_cut_token.exc0[i]
		tree_cut_token.cutCost1[i] = tree_cut_token.exc1[i]
		(tree_cut_token.cutCost[i], pol) = findmin([tree_cut_token.cutCost0[i], tree_cut_token.cutCost1[i]])
		tree_cut_token.polarity[i] = pol-1

		if pol == 0
			tree_cut_token.areaUtil0[i] = vtxCuts[i]	
			tree_cut_token.areaUtil1[i] = total_vwts - vtxCuts[i]
            
		else
			tree_cut_token.areaUtil1[i] = vtxCuts[i]
			tree_cut_token.areaUtil0[i] = total_vwts - vtxCuts[i]
		end

        if tree_cut_token.areaUtil0[i] > capacities[1] || tree_cut_token.areaUtil1[i] > capacities[2]
            tree_cut_token.areaCost[i] = 1e09
        else
            push!(feasible_cuts, tree_cut_token.cutCost[i])
            push!(feasible_edges, [i, pred[i]])
        end

		tree_cut_token.ratioCost[i] = tree_cut_token.cutCost[i]/(tree_cut_token.areaUtil0[i] * tree_cut_token.areaUtil1[i])
		tree_cut_token.tCost_v[i] = tree_cut_token.cutCost[i] + tree_cut_token.areaCost[i]
	end

    return feasible_cuts, feasible_edges
end
