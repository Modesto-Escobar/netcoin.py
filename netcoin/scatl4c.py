"""looking4clusters and communities inside a scattergram.

Ports of ``surScatL4C.R`` (``l4cLayout``, ``l4cClusters``, ``l4cRun``, ``addL4C``, ``alignPlane``) and
``surScatCommunities.R`` (``addCommunities``): what :func:`sur_scat` adds, besides its own plane and groups,
when it is given ``layouts`` and ``clusters``.
"""
from __future__ import annotations

import logging
import random
import warnings

import numpy as np
import pandas as pd

from .core import netCoin
from .rotations import ROTATIONS, rotation_name
from .scattools import (_copy, add_axes, add_clusters, current_cluster_means, current_clusters, current_layouts,
                        layouts_align, list_l4c_clusters, sampled_note)

log = logging.getLogger("netcoin")

L4C_REDUCTIONS = ("pca", "tsne", "cmds", "nmf", "umap")
L4C_METHODS = ("kmeans", "pam", "hclust")
L4C_MAX_CASES = 5000


# ===================================================================== what is asked for ====
def l4c_layout(layouts) -> list[str]:
    """The reductions that ``layouts`` asks for, in the order asked (``True`` asks for all, ``"mds"`` is ``"cmds"``)."""
    if layouts is None or layouts is False:
        return []
    if layouts is True:
        return list(L4C_REDUCTIONS)
    if not isinstance(layouts, (str, list, tuple)):
        warnings.warn("surScat has no links, so a network layout does not apply to it and was ignored. layouts takes "
                      "the names of dimensionality reductions: " + ", ".join(L4C_REDUCTIONS))
        return []
    layouts = [layouts] if isinstance(layouts, str) else list(layouts)
    asked = ["cmds" if a == "mds" else a for a in rotation_name(layouts)]
    if "all" in asked:
        layouts = [l for l, a in zip(layouts, asked) if a != "all"] + list(L4C_REDUCTIONS)
        asked = [a for a in asked if a != "all"] + list(L4C_REDUCTIONS)
    own = [a in ROTATIONS for a in asked]
    layouts = [l for l, o in zip(layouts, own) if not o]
    asked = [a for a, o in zip(asked, own) if not o]
    bad = [l for l, a in zip(layouts, asked) if a not in L4C_REDUCTIONS]
    if bad:
        warnings.warn("surScat has no links, so network layouts do not apply to it: "
                      + ", ".join(f'"{b}"' for b in bad) + " ignored. layouts takes the names of dimensionality "
                      "reductions: " + ", ".join(L4C_REDUCTIONS) + ", and the rotations " + ", ".join(ROTATIONS)
                      + ', or "all" for all of them')
    return list(dict.fromkeys(a for a in asked if a in L4C_REDUCTIONS))


def rotation_layout(layouts) -> list[str]:
    """The rotations asked for (``"all"`` stands for the reductions and the rotations)."""
    if not isinstance(layouts, (str, list, tuple)):
        return []
    asked = rotation_name(list(layouts) if not isinstance(layouts, str) else layouts)
    if "all" in asked:
        asked = [a for a in asked if a != "all"] + list(ROTATIONS)
    return list(dict.fromkeys(a for a in asked if a in ROTATIONS))


def l4c_clusters(clusters) -> list[str]:
    """The clustering methods asked for; ``True`` stands for all of them and the latent model that suits the analysis."""
    if clusters is None or clusters is False:
        return []
    if clusters is True:
        return list(L4C_METHODS) + ["latent"]
    known = list(L4C_METHODS) + ["lpa", "lca"]
    if not isinstance(clusters, (str, list, tuple)):
        raise ValueError("clusters must be True, False or the names of clustering methods: " + ", ".join(known))
    given = [clusters] if isinstance(clusters, str) else list(clusters)
    asked = [str(c).lower() for c in given]
    bad = [c for c, a in zip(given, asked) if a not in known]
    if bad:
        raise ValueError("clusters takes the names of clustering methods (" + ", ".join(known) + "), not "
                         + ", ".join(f'"{b}"' for b in bad))
    return list(dict.fromkeys(asked))


def latent_clusters(methods: list[str], kind: str) -> list[str]:
    """Which latent model runs, if any, once the type of analysis is known."""
    from .latent import latent_model

    model = latent_model(kind)
    if "latent" in methods:
        return [model]
    named = [m.upper() for m in methods if m in ("lpa", "lca")]
    if not named:
        return []
    wrong = [m for m in named if m != model]
    if wrong:
        warnings.warn(", ".join(wrong) + " left out: " + ('LCA needs qualitative active variables (type="mca")'
                      if model == "LPA" else 'LPA needs quantitative active variables (type="pca")'))
    return [model] if model in named else []


# ============================================================================== the run ====
def l4c_run(data, reductions: list[str], methods: list[str], nclusters: list[int], seed=None,
            force_execution: bool = False) -> dict:
    """Run looking4clusters on the cases and return the planes and groups asked for, case by case."""
    from .clusters import (_create, run_kmeans, run_mds, run_nmf, run_pam_hclust, run_pca, run_tsne, run_umap)

    frame = pd.DataFrame(data).reset_index(drop=True)
    n = len(frame)
    obj = _create(frame, False)
    big = n > L4C_MAX_CASES and not force_execution
    skipped: list[str] = []

    def attempt(what, function, *args, **kwargs):
        try:
            function(obj, *args, **kwargs)
        except Exception as exc:
            warnings.warn(f"looking4clusters could not compute {what}: {exc}")

    ks = [int(k) for k in dict.fromkeys(nclusters) if 2 <= k < n]
    if methods and ks:
        pam_hclust = [m for m in methods if m in ("pam", "hclust")]
        if big and pam_hclust:
            skipped += pam_hclust
            pam_hclust = []
        rest = sorted(ks)
        while rest:
            sk = rest[0] + 5
            if "kmeans" in methods:
                attempt("kmeans", run_kmeans, selectedk=sk, seed=seed)
            if pam_hclust:
                attempt(" and ".join(pam_hclust), run_pam_hclust, selectedk=sk)
            rest = [k for k in rest if k > sk + 5]

    clusters: dict = {}
    if obj.clusters:
        found = list_l4c_clusters(obj)
        width = len(str(max(ks)))
        for meth in [m for m in methods if m in L4C_METHODS and m not in skipped]:
            for k in ks:
                key = f"{meth}({k})"
                if key in found:
                    clusters[f"{meth}({k:0{width}d})"] = found[key]

    steps = {"pca": lambda o: run_pca(o), "tsne": lambda o: run_tsne(o, seed=seed), "cmds": lambda o: run_mds(o),
             "nmf": lambda o: run_nmf(o, seed=seed), "umap": lambda o: run_umap(o, seed=seed)}
    for r in reductions:
        if big and r == "cmds":
            skipped.append(r)
            continue
        attempt(r, lambda o, r=r: steps[r](o))
    axes = {r: obj.reductions[r].reset_index(drop=True) for r in reductions if r in obj.reductions}
    # the pca of looking4clusters works on the covariances of the variables as they come, not on the
    # standardized ones surScat's own PCA plane comes from, so its plane and axes are told apart
    named = {}
    for r, plane in axes.items():
        pct = obj.pctvar.get(r)
        plane = plane.copy()
        if r == "pca":
            plane.columns = ["cPC1", "cPC2"]
        if pct is not None and len(pct) == 2 and np.all(np.isfinite(pct)):
            plane.columns = [f"{c} ({p:.1f}%)" for c, p in zip(plane.columns, pct)]
        named["cPC" if r == "pca" else r] = plane
    if skipped:
        many = len(skipped) > 1
        warnings.warn(f"with {n} cases, more than {L4C_MAX_CASES}, " + ", ".join(skipped) + " would take too long and "
                      + ("were" if many else "was") + " left out. State force_execution=True to compute "
                      + ("them" if many else "it") + " anyway.")
    return {"axes": named, "clusters": clusters}


def align_plane(align_to, planes: list[str]) -> int:
    """The 1-based index of the plane ``align_to`` points to (1, the main plane, when none)."""
    if align_to is None:
        return 1
    hit: list[int] = []
    if isinstance(align_to, (int, np.integer)) and not isinstance(align_to, bool):
        hit = [int(align_to)] if 1 <= align_to <= len(planes) else []
    else:
        keys = [str(align_to).lower()] + rotation_name(str(align_to))
        low = [p.lower() for p in planes]
        hit = [i + 1 for i, p in enumerate(low) if p in keys]
        if not hit:
            hit = [i + 1 for i, p in enumerate(low) if any(p.startswith(k + " ") for k in keys)]
    if not hit:
        shown = f'"{align_to}"' if isinstance(align_to, str) else align_to
        warnings.warn(f"alignTo={shown} is none of the planes ({', '.join(planes)}), so they are aligned to the "
                      "main plane")
        return 1
    return hit[0]


def add_l4c(net, found: dict, idx=None, kept=None, case_weight=None, sort: bool = True, suffix=None, align_to=None,
            axes_labels=None):
    """Add what :func:`l4c_run` found, case by case, to the scattergram (``addL4C``).

    ``idx`` is the 1-based case-to-pattern map (``None`` when the nodes are cases) and ``kept`` the 1-based
    nodes that ``max_n`` drew.
    """
    if not found["axes"] and not found["clusters"]:
        return net
    net = _copy(net)
    first = next(iter(found["clusters"].values())) if found["clusters"] else next(iter(found["axes"].values()))
    sel = np.arange(len(first))
    mapping = None if idx is None else np.asarray(idx)
    if kept is not None:
        kept = np.asarray(kept)
        if idx is None:
            sel = kept - 1
        else:
            sel = np.where(np.isin(mapping, kept))[0]
            lookup = {int(v): i + 1 for i, v in enumerate(kept)}
            mapping = np.array([lookup[int(v)] for v in mapping[sel]])
    saved = {"case_to_pattern": net.scat.get("case_to_pattern"), "case_weight": net.scat.get("case_weight")}
    if mapping is None:
        net.scat.pop("case_to_pattern", None)
    else:
        net.scat["case_to_pattern"] = mapping
    if mapping is not None and case_weight is not None:
        net.scat["case_weight"] = np.asarray(case_weight)[sel]
    else:
        net.scat.pop("case_weight", None)
    if found["clusters"]:
        frame = pd.DataFrame({k: np.asarray(v)[sel] for k, v in found["clusters"].items()})
        net = add_clusters(net, frame, sort=sort, max_groups=None)
    if found["axes"]:
        net = add_axes(net, {k: v.iloc[sel].reset_index(drop=True) for k, v in found["axes"].items()},
                       axes_labels=axes_labels)
    for key, value in saved.items():
        if value is None:
            net.scat.pop(key, None)
        else:
            net.scat[key] = value
    if net.scat.get("cluster_columns"):
        net.clusters = current_clusters(net)
        net.cluster_means = current_cluster_means(net)
    if suffix is not None and found["axes"]:
        planes = list(current_layouts(net))
        net = layouts_align(net, "procrustes", to=align_plane(align_to, planes), suffix=suffix)
    return net


# ================================================================================ communities ====
COMMUNITY_ALGORITHMS = {
    "ed": ("edgeBetweenness", 2000), "fa": ("fastGreedy", np.inf), "la": ("labelProp", np.inf),
    "le": ("leadingEigen", np.inf), "lo": ("louvain", np.inf), "op": ("optimal", 100),
    "sp": ("spinglass", 2000), "wa": ("walktrap", np.inf),
}
RANDOM_METHODS = ("la", "lo", "sp")


def community_code(method) -> str:
    """The two-letter code of a community method given as code, igraph name or its name here."""
    if not isinstance(method, str):
        raise ValueError("each community method must be a single name")
    key = method.lower()
    if key.startswith("cluster_"):
        key = key[len("cluster_"):]
    full = {v[0].lower(): k for k, v in COMMUNITY_ALGORITHMS.items()}
    if key in full:
        return full[key]
    if key[:2] in COMMUNITY_ALGORITHMS:
        return key[:2]
    raise ValueError(f'"{method}" is not a community detection method. Available: '
                     + ", ".join(f"{k} ({v[0]})" for k, v in COMMUNITY_ALGORITHMS.items()))


def lump_communities(membership, min_size: float, max_communities: int | None = None) -> dict:
    """Gather into one residual group the communities with fewer cases than ``min_size`` (``lumpCommunities``).

    ``min_size`` is a proportion of the cases when below 1 and a count of them otherwise.
    """
    membership = np.asarray(membership)
    n = len(membership)
    threshold = max(2, int(np.ceil(min_size * n))) if min_size < 1 else int(np.ceil(min_size))
    sizes = pd.Series(membership).value_counts().sort_index()
    big = [c for c, s in sizes.items() if s >= threshold]
    if not big:
        raise ValueError(f"every community holds fewer than {threshold} cases, so all of them would be gathered into "
                         "the residual group. Lower min_size, or raise min_link so that the correlation network is "
                         "sparser and its communities larger.")
    ordered = sorted(big, key=lambda c: (-sizes[c], c))
    n_small = len(sizes) - len(big)
    n_ranked = 0
    if max_communities is not None and len(ordered) + (n_small > 0) > max_communities:
        keep = int(max_communities) - 1
        n_ranked = len(ordered) - keep
        ordered = ordered[:keep]
    number = {c: i + 1 for i, c in enumerate(ordered)}
    out = np.array([number.get(c, 0) for c in membership])
    lumped = out == 0
    rest = None
    if lumped.any():
        rest = len(ordered) + 1
        out[lumped] = rest
    return {"clusters": out, "rest": rest, "n_groups": int(out.max()), "n_lumped": int(lumped.sum()),
            "n_small": int(n_small), "n_ranked": int(n_ranked), "threshold": threshold}


def add_communities(net, data: pd.DataFrame, variables=None, community="lo", scale: bool = True,
                    method: str = "pearson", criteria: str = "value", min_link: float = 0.95,
                    max_link: float = np.inf, pairwise: bool = False, min_size: float = 0.01,
                    max_communities: int | None = None, name: str | None = None, sort: bool = True, weight=None,
                    max_groups: int | None = 25, max_nodes: float | None = None, seed: int | None = 2020,
                    weighted: bool = True):
    """Add a cluster column per community-detection method, from the network of correlations between the cases.

    Port of ``addCommunities``. The cases are correlated with one another over ``variables`` (the transposed
    reading of the table); the communities of that network are a grouping of the cases, and those smaller than
    ``min_size`` (a proportion below 1, a count otherwise) are gathered into a residual group.
    ``community`` takes codes (``"lo"``) or names (``"louvain"``): edgeBetweenness, fastGreedy, labelProp,
    leadingEigen, louvain, optimal, spinglass, walktrap. Needs ``python-igraph``.

    The links of the network are weighed by their correlation (what the R code means to do). R's ``toIgraph``
    does not carry that weight over from the network of ``netCorr``, so R finds the communities of the
    unweighted network; ``weighted=False`` does the same here.
    """
    from .igraph_io import to_igraph
    from .layouts import compute_membership, require_igraph
    from .netcorr import net_corr

    require_igraph()
    if not isinstance(net, netCoin):
        raise TypeError("net must be a netCoin object returned by sur_scat")
    if not isinstance(data, pd.DataFrame):
        raise ValueError("data must be the data frame the scattergram was drawn from")
    if min_size <= 0:
        raise ValueError("min_size must be a proportion of the cases (below 1) or a count of them (1 or more)")
    if max_communities is not None:
        if max_communities < 2:
            raise ValueError("max_communities must be a count of two groups or more, the residual one among them, "
                             "or None to leave as many groups as clear min_size")
        max_communities = int(max_communities)
    codes = [community_code(c) for c in ([community] if isinstance(community, str) else list(community))]
    if len(set(codes)) != len(codes):
        raise ValueError("community holds the same method more than once: "
                         + ", ".join(sorted({c for c in codes if codes.count(c) > 1})))
    variables = list(data.columns) if variables is None else list(variables)
    missing = [v for v in variables if v not in data.columns]
    if missing:
        raise ValueError("variables not found in data: " + ", ".join(missing))
    if len(variables) < 3:
        warnings.warn(f"the cases are correlated over {len(variables)} variables, so each correlation rests on "
                      f"{len(variables)} points and says very little. Compare the cases over more variables.")
    columns = {}
    for v in variables:
        x = data[v]
        if isinstance(x.dtype, pd.CategoricalDtype):
            if not x.cat.ordered:
                raise ValueError(f"{v} is a nominal variable, so correlating the cases over it would impose an "
                                 "arbitrary order on its categories. Turn it into an ordered factor or a numeric one.")
            columns[v] = (x.cat.codes + 1).astype(float).where(x.notna())
        elif pd.api.types.is_numeric_dtype(x) and not pd.api.types.is_bool_dtype(x):
            columns[v] = x.astype(float)
        else:
            raise ValueError(f"{v} must be numeric or an ordered factor to compare the cases over it")
    D = pd.DataFrame(columns).dropna().reset_index(drop=True)
    n_nodes = len(net.nodes)
    idx, sampled, origin = net.scat.get("case_to_pattern"), net.scat.get("sampled_nodes"), net.scat.get("sampled_from")
    if sampled is not None:
        if origin != len(D):
            raise ValueError(f"the {n_nodes} nodes are a sample of {origin} patterns, not of the {len(D)} complete "
                             "cases of data, so which cases each node stands for is no longer recorded and the "
                             "communities cannot be lined up with them. Raise max_n in sur_scat so that every "
                             "pattern is drawn, or drop patterns.")
        D = D.iloc[np.asarray(sampled) - 1].reset_index(drop=True)
    elif idx is not None:
        if len(idx) != len(D):
            raise ValueError(f"data holds {len(D)} complete cases but the scattergram was drawn from {len(idx)}. "
                             "Pass the same data frame and the same variables sur_scat was given.")
    elif len(D) != n_nodes:
        raise ValueError(f"data holds {len(D)} complete cases but the scattergram has {n_nodes} nodes. Pass the "
                         "same data frame and the same variables sur_scat was given." + sampled_note(net))
    if len(D) < 3:
        raise ValueError("at least three cases are needed to build a network of them")
    if scale:
        constant = [c for c in D.columns if D[c].std(ddof=1) == 0]
        if constant:
            raise ValueError("these variables take a single value across the cases, so they cannot be "
                             "standardized: " + ", ".join(constant) + ". Drop them, or set scale=False.")
        D = (D - D.mean()) / D.std(ddof=1)
    width = len(str(len(D)))
    transposed = pd.DataFrame(D.to_numpy().T, columns=[f"{i + 1:0{width}d}" for i in range(len(D))])
    corr = net_corr(transposed, method=method, criteria=criteria, min_link=min_link, max_link=max_link,
                    pairwise=pairwise, link_bipolar=False)
    graph = to_igraph(corr)
    if graph.vcount() != len(D):
        raise ValueError(f"the correlation network came back with {graph.vcount()} nodes for {len(D)} cases")
    if not weighted and "weight" in graph.es.attributes():
        del graph.es["weight"]
    weights = graph.es["weight"] if "weight" in graph.es.attributes() else None
    if weights is not None and len(weights) and min(weights) < 0:
        warnings.warn("the network holds links of negative correlation, which the community methods cannot weigh, "
                      "so the communities were found without weighing them. Raise min_link above 0 to keep only "
                      "links of positive correlation.")
        graph.es["weight"] = [1.0] * graph.ecount()
    ids = graph.vs["name"] if "name" in graph.vs.attributes() else None
    order = np.arange(len(D)) if ids is None else np.array([ids.index(f"{i + 1:0{width}d}") for i in range(len(D))])
    reports, ranked = [], False
    out = net
    for code in codes:
        label, limit = COMMUNITY_ALGORITHMS[code]
        limit = limit if max_nodes is None else max_nodes
        if graph.vcount() > limit:
            allowed = [v[0] for v in COMMUNITY_ALGORITHMS.values() if np.isinf(v[1])]
            raise ValueError(f"{label} does not answer on a network of {graph.vcount()} nodes; it bears {limit}. "
                             f"Choose one of {', '.join(allowed)}, cap the nodes with max_n in sur_scat, or raise "
                             "max_nodes to wait for it anyway.")
        if seed is not None:
            random.seed(seed)
        try:
            membership = np.asarray(compute_membership(graph, code))[order]
        except Exception as exc:
            isolated = sum(1 for d in graph.degree() if d == 0)
            tail = (f" It holds {isolated} isolated cases, which some methods refuse. Lower min_link so that fewer "
                    "cases are left without links, or choose another method.") if isolated else ""
            raise ValueError(f"{label} could not find communities in this network: {exc}{tail}") from exc
        lump = lump_communities(membership, min_size, max_communities)
        own = f"{label}({lump['n_groups']})"
        column = own if name is None else (f"{name}.{own}" if len(codes) > 1 else name)
        coded = [f"{c:0{len(str(lump['n_groups']))}d}" for c in lump["clusters"]]
        out = add_clusters(out, pd.DataFrame({column: coded}), sort=sort, weight=weight, max_groups=max_groups,
                           rest=lump["rest"])
        if lump["n_lumped"]:
            why = []
            if lump["n_small"]:
                why.append(f"{lump['n_small']} of fewer than {lump['threshold']} cases")
            if lump["n_ranked"]:
                why.append(f"{lump['n_ranked']} past the {max_communities - 1} largest")
            reports.append(f"{column}: {lump['n_small'] + lump['n_ranked']} communities ({', '.join(why)}), "
                           f"{lump['n_lumped']} cases in all")
            ranked = ranked or lump["n_ranked"] > 0
    if reports:
        message = ("Gathered into the residual group -- " + "; ".join(reports) + ". Set "
                   + ("min_size=1 and max_communities=None " if ranked else "min_size=1 ")
                   + "to keep every community apart.")
        log.info(message)
    return out
