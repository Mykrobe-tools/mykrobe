from operator import itemgetter

import anytree

class LineagePredictor:
    def __init__(self, variant_to_lineage):
        """lineage_calls should be a dictionary of variant name -> lineage name.
        eg X42Y -> "lineageM.N"."""
        self.variant_to_lineage = variant_to_lineage
        self._make_tree()
        self.result = {}

    def _make_tree(self):
        self.tree_nodes = {}
        self.tree_root = anytree.Node("root")

        for variant, lineage in self.variant_to_lineage.items():
            pieces = lineage["name"].split(".")
            nodes_to_add = [".".join(pieces[:i+1]) for i in range(len(pieces))]

            for i, node_name in enumerate(nodes_to_add):
                if node_name in self.tree_nodes:
                    continue
                parent = self.tree_root if i == 0 else self.tree_nodes[nodes_to_add[i-1]]
                new_node = anytree.Node(node_name, parent=parent)
                self.tree_nodes[node_name] = new_node


    def _score_each_lineage_node(self, lineage_calls):
        scores = {}

        for lineage_name, var_dict in lineage_calls.items():
            best_score = None
            for var_name, call in var_dict.items():
                if self.variant_to_lineage[var_name]["use_ref_allele"]:
                    wanted_genotype = [0, 0]
                else:
                    wanted_genotype = [1, 1]

                if call["genotype"] == wanted_genotype:
                    multiplier = 1 if len(call["info"]["filter"]) == 0 else 0.5
                elif call["genotype"] == [0, 1]:
                    multiplier = 0.5
                else:
                    multiplier = -1

                score = multiplier * call["info"]["conf"]
                if best_score is None or best_score < score:
                    best_score = score

            scores[lineage_name] = 0 if best_score is None else best_score

        return scores


    def _get_paths_and_scores(self, lineage_calls):
        paths = []
        used_lineages = set()
        node_scores = self._score_each_lineage_node(lineage_calls)

        for leaf_node in self.tree_root.leaves:
            path = leaf_node.path[1:]
            path_scores = [(x.name, node_scores.get(x.name, 0)) for x in path]
            path_leaf = leaf_node
            while len(path_scores) > 0 and path_scores[-1][1] <= 0:
                path_scores.pop()
                path_leaf = path_leaf.parent
            if len(path_scores) > 0 and path_leaf.name not in used_lineages:
                paths.append({
                    "score": sum([x[1] for x in path_scores]),
                    "lineage": path_leaf.name,
                    "scores": {x[0]: x[1] for x in path_scores}})
                used_lineages.add(path_leaf.name)

        paths.sort(key=itemgetter("score"), reverse=True)
        return paths


    # This is just a simple implementation for initial testing. It just finds
    # the single best lineage by looking at the best path in the tree based
    # on conf scores from the genotyping. Expect to replace this with
    # something more sophisitcated, that can identify mixed samples and report
    # a mix of lineages.
    def call_lineage(self, lineage_calls):
        """Calls the lineage from the given lineage calls.  lineage_calls should
        be the dictionary from a Genotype instance lineage_calls_dict. Has
        lineage name -> {dict of variant name -> call}."""
        paths = self._get_paths_and_scores(lineage_calls)

        if len(paths) == 0:
            return None

        paths.sort(key=itemgetter("score"), reverse=True)
        best_paths = [x for x in paths if x["score"] == paths[0]["score"]]
        used_calls = {}
        result = {
            "lineage": [x["lineage"] for x in best_paths],
            "calls": used_calls,
        }

        for path_dict in best_paths:
            used_calls[path_dict["lineage"]] = {}
            for lineage in path_dict["scores"]:
                used_calls[path_dict["lineage"]][lineage] = lineage_calls[lineage]

        return result

