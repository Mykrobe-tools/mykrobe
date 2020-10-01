from operator import itemgetter

import anytree

class LineagePredictor:
    def __init__(self, variant_to_lineage):
        """lineage_calls should be a dictionary of variant name -> lineage name.
        eg X42Y -> "lineageM.N"."""
        self.variant_to_lineage = variant_to_lineage
        self.report_names = {d["name"]: d.get("report_name", d["name"]) for d in variant_to_lineage.values()}
        self._make_tree()
        self.result = {}


    @classmethod
    def _lineage_to_itself_plus_parents(cls, lineage):
        pieces = lineage.split(".")
        return [".".join(pieces[:i+1]) for i in range(len(pieces))]


    def _make_tree(self):
        self.tree_nodes = {}
        self.tree_root = anytree.Node("root")

        for variant, lineage in self.variant_to_lineage.items():
            nodes_to_add = LineagePredictor._lineage_to_itself_plus_parents(lineage["name"])

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
    def call_lineage_using_conf_scores(self, lineage_calls):
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
                used_calls[path_dict["lineage"]][lineage] = lineage_calls.get(lineage, None)

        return result


    def _genotype_each_lineage_node(self, lineage_calls):
        genotypes = {}

        for lineage_name, var_dict in lineage_calls.items():
            best_geno = None
            for var_name, call in var_dict.items():
                if self.variant_to_lineage[var_name]["use_ref_allele"]:
                    wanted_genotype = [0, 0]
                else:
                    wanted_genotype = [1, 1]

                if call["genotype"] == wanted_genotype:
                    geno = 1
                elif call["genotype"] == [0, 1]:
                    geno = 0.5
                else:
                    geno = 0

                if best_geno is None or best_geno < geno:
                    best_geno = geno

            genotypes[lineage_name] = 0 if best_geno is None else best_geno

        return genotypes

    def _get_good_paths_using_genotype_calls(self, lineage_calls, min_frac_called=0.5):
        paths = {}
        node_genos = self._genotype_each_lineage_node(lineage_calls)

        for leaf_node in self.tree_root.leaves:
            path = leaf_node.path[1:]
            path_genos = [(x.name, node_genos.get(x.name, 0)) for x in path]
            path_leaf = leaf_node
            while len(path_genos) > 0 and path_genos[-1][1] == 0:
                path_genos.pop()
                path_leaf = path_leaf.parent

            number_good_nodes = len([x for x in path_genos if x[1] > 0])
            if len(path_genos) == 0 or number_good_nodes / len(path_genos) < min_frac_called or path_leaf.name in paths:
                continue

            # We could end up having say [l1, l1.1, l1.1.2], and also just
            # [l1] because the leaf l1.2 had no support.  Could come across these
            # in either order.
            # Check if new path is a subpath of one we've already found:
            if any(path_leaf.name in x for x in paths):
                continue

            subpaths_to_remove = [x for x in paths if x in path_leaf.name]
            for x in subpaths_to_remove:
                del paths[x]

            paths[path_leaf.name] = {
                "good_nodes": number_good_nodes,
                "tree_depth": len(path_genos),
                "genotypes": {x[0]: x[1] for x in path_genos}}

        return paths


    def replace_dict_keys(self, d):
        return {self.report_names.get(k, k): v for k, v in d.items()}

    # The initial implementation assumed that the user would have consistent
    # lineage names with X.Y being a child of X. Unfortunately, for TB we
    # actually have lineage2.2.8 being a child of lineage2.2.7. This class
    # uses the dots in the names to define the tree. So we now have a
    # "report name", which is what the user will get in the output JSON
    # file. This function goes through all the results and changes all
    # the relevant lineage names to their "report name" where they are different
    # from the internal name that defines the lineage tree.
    def _apply_report_names_to_lineage_call_result(self, result):
        result["lineage"] = [self.report_names.get(x, x) for x in result["lineage"]]
        result["calls"] = self.replace_dict_keys(result["calls"])
        for lineage, d in result["calls"].items():
            result["calls"][lineage] = self.replace_dict_keys(d)
        result["calls_summary"] = self.replace_dict_keys(result["calls_summary"])
        for d in result["calls_summary"].values():
            d["genotypes"] = self.replace_dict_keys(d["genotypes"])


    def call_lineage(self, lineage_calls, min_frac_called=0.5):
        """Calls the lineage from the given lineage calls.  lineage_calls should
        be the dictionary from a Genotype instance lineage_calls_dict. Has
        lineage name -> {dict of variant name -> call}."""
        paths = self._get_good_paths_using_genotype_calls(lineage_calls, min_frac_called=min_frac_called)

        if len(paths) == 0:
            return None

        used_calls = {}
        result = {
            "lineage": sorted(list(paths.keys())),
            "calls_summary": paths,
            "calls": used_calls,
        }

        for lineage, path_dict in paths.items():
            used_calls[lineage] = {}
            for lineage2 in path_dict["genotypes"]:
                used_calls[lineage][lineage2] = lineage_calls.get(lineage2, None)

        self._apply_report_names_to_lineage_call_result(result)
        return result

