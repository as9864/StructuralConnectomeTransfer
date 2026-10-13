"""E4: do wiring rules fitted on one nervous system explain another?

Generative network model without coordinates (the topological part of the
economic models of Betzel et al. 2016 / Akarca et al. 2021): starting from an
empty graph on n nodes, add m edges one at a time, choosing a non-edge (i, j)
with probability

    P(i~j) ∝ (M_ij + eps)^alpha * ((k_i + 1)(k_j + 1))^beta

where M is the matching index (overlap of the two neighbourhoods, homophily)
and k the current degree (preferential attachment for beta > 0). alpha = beta
= 0 is an Erdos-Renyi graph and serves as the null rule.

Because every target is sampled at the same n and density (as in
run_baselines.py), m is the same for all targets: each rule on a grid of
(alpha, beta) is simulated once and scored against every target. Fit quality
is Betzel's energy, the largest Kolmogorov-Smirnov distance between model and
target distributions of node degree, clustering and betweenness.

Fitting uses half of each target's samples, evaluation the other half, so a
target's own rule is not evaluated on the samples it was fitted to. For a
source S and target G:

    TI(S -> G) = (E_null(G) - E_S(G)) / (E_null(G) - E_self(G))

1 = S's rule explains G as well as G's own rule, 0 = no better than ER.

Usage:
    python -m sctransfer.run_rule_transfer
"""
import argparse
from pathlib import Path

import networkx as nx
import numpy as np
import pandas as pd
from scipy.stats import ks_2samp

from sctransfer import normalized as N
from sctransfer.run_baselines import graphs as baseline_graphs
from sctransfer.run_worm_fly_scope import giant_component, grow

RESULTS_DIR = Path(__file__).resolve().parent.parent / "results"
N_NODES = 150
ALPHAS = np.array([-1.0, -0.5, 0.0, 0.5, 1.0, 2.0, 3.0])
BETAS = np.round(np.arange(-0.5, 2.51, 0.1), 2)  # a coarse first pass put every best beta in [1, 1.5]
EPS = 1e-5


def generate(n, m, alpha, beta, rng):
    A = np.zeros((n, n))
    iu = np.triu_indices(n, 1)
    for _ in range(m):
        k = A.sum(1)
        common = A @ A
        denom = k[:, None] + k[None, :] - common - 2 * A
        M = np.divide(common, denom, out=np.zeros_like(common), where=denom > 0)
        score = (M[iu] + EPS) ** alpha * (np.outer(k + 1, k + 1)[iu]) ** beta
        score[A[iu] > 0] = 0
        score = np.nan_to_num(score, posinf=0.0)
        if score.sum() <= 0:
            score = (A[iu] == 0).astype(float)
        pick = rng.choice(len(score), p=score / score.sum())
        i, j = iu[0][pick], iu[1][pick]
        A[i, j] = A[j, i] = 1
    return A


def node_stats(A):
    G = nx.from_numpy_array(A)
    deg = np.array([d for _, d in G.degree()], float)
    clu = np.array(list(nx.clustering(G).values()))
    bet = np.array(list(nx.betweenness_centrality(G).values()))
    return deg, clu, bet


def energy(model, target):
    """model, target: (deg, clu, bet) arrays pooled over graphs."""
    return max(ks_2samp(model[t], target[t]).statistic for t in range(3))


def pool(stats_list):
    return tuple(np.concatenate([s[t] for s in stats_list]) for t in range(3))


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--density", type=float, default=0.04)
    parser.add_argument("--n-samples", type=int, default=20, help="target samples per graph (half fit, half test)")
    parser.add_argument("--n-reps", type=int, default=4, help="model graphs per rule")
    parser.add_argument("--n-boot", type=int, default=20, help="bootstrap refits over each target's fit samples")
    parser.add_argument("--seed", type=int, default=0)
    args = parser.parse_args()
    RESULTS_DIR.mkdir(exist_ok=True)
    rng = np.random.default_rng(args.seed)
    m = int(round(args.density * N_NODES * (N_NODES - 1) / 2))

    # targets: same sampling as run_baselines.py, binarized at the target density
    targets = {}
    for name, (group, W) in baseline_graphs().items():
        if name == "worm_whole_animal":
            continue  # different scope (whole body vs brain), see experimental_design.md 4d
        adj = ((W + W.T) > 0).astype(np.int8).tocsr()
        giant = np.flatnonzero(giant_component(adj))
        samples = []
        while len(samples) < args.n_samples:
            keep = grow(adj, N_NODES, [rng.choice(giant)], rng)
            if keep is None:
                continue
            A = N.binarize(W[keep][:, keep].toarray(), args.density)
            if A.sum() // 2 >= 0.95 * m:
                samples.append(node_stats(A))
        half = args.n_samples // 2
        targets[name] = {"group": group, "fit": pool(samples[:half]), "fit_list": samples[:half],
                         "test": pool(samples[half:])}
        print(f"[rules] target {name}: {args.n_samples} samples", flush=True)

    # simulate every rule once
    rules = [(a, b) for a in ALPHAS for b in BETAS]
    model_stats = {}
    for r, (a, b) in enumerate(rules):
        model_stats[(a, b)] = pool([node_stats(generate(N_NODES, m, a, b, rng)) for _ in range(args.n_reps)])
        if r % 20 == 0:
            print(f"[rules] simulated {r + 1}/{len(rules)} rules", flush=True)

    # energy of every rule on every target, on fit and test halves
    rows = []
    for name, t in targets.items():
        for rule, ms in model_stats.items():
            rows.append({"target": name, "alpha": rule[0], "beta": rule[1],
                         "energy_fit": energy(ms, t["fit"]), "energy_test": energy(ms, t["test"])})
    E = pd.DataFrame(rows)
    E.to_csv(RESULTS_DIR / f"rule_energy_d{args.density:g}.csv", index=False)

    best = E.loc[E.groupby("target")["energy_fit"].idxmin()].set_index("target")[["alpha", "beta", "energy_fit"]]
    best["group"] = [targets[t]["group"] for t in best.index]
    print("\n=== Best rule per graph (fitted on half the samples) ===")
    print(best.round(3).to_string())

    test = E.set_index(["target", "alpha", "beta"])["energy_test"]
    rule_list = list(model_stats)

    def transfer_table(best_rules):
        out = []
        for src in targets:
            for tgt in targets:
                e_null = test[(tgt, 0.0, 0.0)]
                e_self = test[(tgt, *best_rules[tgt])]
                e_src = test[(tgt, *best_rules[src])]
                gs, gt = targets[src]["group"], targets[tgt]["group"]
                kind = "self" if src == tgt else ("within " + gs if gs == gt else "from " + gs + " to " + gt)
                out.append({"source": src, "target": tgt, "kind": kind, "energy_null": e_null,
                            "energy_self": e_self, "energy_source": e_src,
                            "TI": (e_null - e_src) / (e_null - e_self) if e_null > e_self else np.nan})
        return pd.DataFrame(out)

    T = transfer_table({t: tuple(best.loc[t, ["alpha", "beta"]]) for t in targets})
    T.to_csv(RESULTS_DIR / f"rule_transfer_d{args.density:g}.csv", index=False)

    # bootstrap: refit every target on resampled fit samples, recompute TI on the fixed test halves
    boots = []
    for b in range(args.n_boot):
        best_b = {}
        for name, t in targets.items():
            idx = rng.integers(0, len(t["fit_list"]), len(t["fit_list"]))
            fit_b = pool([t["fit_list"][i] for i in idx])
            best_b[name] = min(rule_list, key=lambda r: energy(model_stats[r], fit_b))
        Tb = transfer_table(best_b)
        Tb["boot"] = b
        Tb["source_beta"] = [best_b[s][1] for s in Tb["source"]]
        boots.append(Tb)
    B = pd.concat(boots)
    B.to_csv(RESULTS_DIR / f"rule_transfer_boot_d{args.density:g}.csv", index=False)

    pd.set_option("display.width", 200)
    print("\n=== Transfer index by kind: point estimate (median over pairs), bootstrap 95% interval ===")
    nonself = T[T.kind != "self"]
    per_boot = B[B.kind != "self"].groupby(["kind", "boot"])["TI"].median().groupby("kind")
    summary = pd.DataFrame({"TI": nonself.groupby("kind")["TI"].median(),
                            "lo": per_boot.quantile(0.025), "hi": per_boot.quantile(0.975),
                            "pairs": nonself.groupby("kind").size()})
    print(summary.round(2).to_string())
    print("\nBootstrap best beta per group (median, 2.5-97.5%):")
    bb = B[B.kind == "self"].assign(group=lambda x: [targets[s]["group"] for s in x["source"]])
    print(bb.groupby("group")["source_beta"].describe(percentiles=[0.025, 0.5, 0.975]).round(2).to_string())
    print("\nEnergies (test half): self vs null, per target")
    print(T[T.kind == "self"][["target", "energy_self", "energy_null"]].round(3).to_string(index=False))
    land = E.pivot_table(index=["alpha", "beta"], columns="target", values="energy_fit")
    print("\n=== Energy-landscape similarity (Spearman correlation across all rules) ===")
    print(land.corr(method="spearman").round(2).to_string())
    print(f"\nSaved results/rule_energy_d{args.density:g}.csv, results/rule_transfer_d{args.density:g}.csv")


if __name__ == "__main__":
    main()
