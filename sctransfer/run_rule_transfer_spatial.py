"""E4 with a spatial term: the economic generative model with neuron positions.

Extends run_rule_transfer.py (topology only) with a wiring-cost term:

    P(i~j) ∝ D_ij^eta * (M_ij + eps)^alpha * ((k_i + 1)(k_j + 1))^beta

D is the Euclidean distance between the two neurons, divided by the mean
pairwise distance of that sample so that eta means the same thing in a worm
(atlas units, ~100 um body) and a fly (nm). eta < 0 penalizes long
connections. Energy adds the KS distance of (normalized) edge lengths to the
degree / clustering / betweenness terms, as in Betzel et al. (2016).

Positions differ per sample, so every rule is simulated on every target
sample (fit half and test half), not once for all targets as in the
topological version. Neurons without a position are dropped before sampling.

Usage:
    python -m sctransfer.run_rule_transfer_spatial
"""
import argparse
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.spatial.distance import pdist, squareform
from scipy.stats import ks_2samp

from sctransfer import coords
from sctransfer import normalized as N
from sctransfer.run_baselines import graphs as baseline_graphs
from sctransfer.run_rule_transfer import EPS, N_NODES, node_stats
from sctransfer.run_worm_fly_scope import giant_component, grow

RESULTS_DIR = Path(__file__).resolve().parent.parent / "results"
ETAS = [-3.0, -2.0, -1.5, -1.0, -0.5, 0.0]
ALPHAS = [0.0, 1.0]  # homophily was never selected on finer grids (alpha = 0 everywhere)
BETAS = [0.5, 0.75, 1.0, 1.25, 1.5, 1.75, 2.0]
FINE_ETAS = [-2.0, -1.5, -1.25, -1.0, -0.75, -0.5, -0.25, 0.0]
NULL_RULE = (0.0, 0.0, 0.0)
D_FLOOR = 0.05  # neurons sharing a position (approximate fly centroids) get a small, finite distance


def generate(D, m, eta, alpha, beta, rng):
    n = D.shape[0]
    iu = np.triu_indices(n, 1)
    dist_term = D[iu] ** eta
    A = np.zeros((n, n))
    for _ in range(m):
        k = A.sum(1)
        common = A @ A
        denom = k[:, None] + k[None, :] - common - 2 * A
        M = np.divide(common, denom, out=np.zeros_like(common), where=denom > 0)
        score = dist_term * (M[iu] + EPS) ** alpha * (np.outer(k + 1, k + 1)[iu]) ** beta
        score[A[iu] > 0] = 0
        score = np.nan_to_num(score, posinf=0.0)
        if score.sum() <= 0:
            score = (A[iu] == 0).astype(float)
        pick = rng.choice(len(score), p=score / score.sum())
        A[iu[0][pick], iu[1][pick]] = A[iu[1][pick], iu[0][pick]] = 1
    return A


def stats(A, D):
    deg, clu, bet = node_stats(A)
    iu = np.triu_indices(A.shape[0], 1)
    return deg, clu, bet, D[iu][A[iu] > 0]


def energy(model, target):
    return max(ks_2samp(model[t], target[t]).statistic for t in range(4))


def pool(lst):
    return tuple(np.concatenate([s[t] for s in lst]) for t in range(4))


def position_lookup(fly_position):
    """fly_position "soma": cell bodies; "synapse": larval synapse centroids
    (exact, from CATMAID connectors) and adult synapse centroids approximated
    from neuropil synapse counts (coords.flywire_synapse_centroids)."""
    pos = {}
    pos.update({("worm", k): v for k, v in coords.worm_positions().items()})
    larva_kind = "soma" if fly_position == "soma" else "synapse"
    pos.update({("larva", k): v for k, v in coords.larva_positions(larva_kind).items()})
    adult = {"soma": lambda: coords.flywire_positions("soma"),
             "synapse": coords.flywire_synapse_centroids,
             "synapse_exact": coords.flywire_exact_synapse_centroids}[fly_position]()
    pos.update({("adult_fly", k): v for k, v in adult.items()})
    return pos


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--density", type=float, default=0.04)
    parser.add_argument("--n-samples", type=int, default=10, help="per graph: half fit, half test")
    parser.add_argument("--n-boot", type=int, default=20)
    parser.add_argument("--n-reps", type=int, default=1,
                        help="model graphs per rule and sample; energy = their mean (removes simulation noise)")
    parser.add_argument("--fly-position", choices=["soma", "synapse", "synapse_exact"], default="synapse",
                        help="synapse = approximate adult centroids; synapse_exact = from the full synapse table")
    parser.add_argument("--fine-eta", action="store_true", help="eta in steps of 0.25 (homophily fixed at 0)")
    parser.add_argument("--reuse", action="store_true",
                        help="keep energies already saved for this tag and simulate only missing rules "
                             "(samples are redrawn identically from the same seed)")
    parser.add_argument("--seed", type=int, default=0)
    args = parser.parse_args()
    RESULTS_DIR.mkdir(exist_ok=True)
    rng = np.random.default_rng(args.seed)
    m = int(round(args.density * N_NODES * (N_NODES - 1) / 2))
    pos = position_lookup(args.fly_position)
    tag = f"{args.fly_position}_d{args.density:g}"
    etas, alphas = (FINE_ETAS, [0.0]) if args.fine_eta else (ETAS, ALPHAS)
    if args.fine_eta:
        tag += "_fine"
    if args.n_reps > 1:
        tag += f"_r{args.n_reps}"
    rules = [(e, a, b) for e in etas for a in alphas for b in BETAS]
    if NULL_RULE not in rules:
        rules.append(NULL_RULE)  # the ER null every transfer index is measured against
    energy_path = RESULTS_DIR / f"rule_energy_spatial_{tag}.csv"
    done = set()
    if args.reuse and energy_path.exists():
        old = pd.read_csv(energy_path)
        done = set(map(tuple, old[["eta", "alpha", "beta"]].drop_duplicates().to_numpy()))

    # samples: (target stats, D) per graph
    targets = {}
    for name, (group, W, ids) in baseline_graphs(with_ids=True).items():
        if ids is None:
            continue
        key = [(group, i.item() if hasattr(i, "item") else i) for i in ids]
        has = np.array([k in pos for k in key])
        W, key = W[has][:, has], [k for k, h in zip(key, has) if h]
        P = np.array([pos[k] for k in key], dtype=float)
        adj = ((W + W.T) > 0).astype(np.int8).tocsr()
        giant = np.flatnonzero(giant_component(adj))
        samples = []
        while len(samples) < args.n_samples:
            keep = grow(adj, N_NODES, [rng.choice(giant)], rng)
            if keep is None:
                continue
            A = N.binarize(W[keep][:, keep].toarray(), args.density)
            if A.sum() // 2 < 0.95 * m:
                continue
            D = squareform(pdist(P[keep]))
            D = np.maximum(D / D[np.triu_indices(N_NODES, 1)].mean(), D_FLOOR)
            np.fill_diagonal(D, 1.0)
            samples.append({"stats": stats(A, D), "D": D})
        targets[name] = {"group": group, "samples": samples}
        print(f"[spatial] {name}: {has.sum():,}/{len(has):,} neurons with positions", flush=True)

    # energy of every rule on every sample (one model graph per rule and sample)
    half = args.n_samples // 2
    rows = []
    for t, (name, tgt) in enumerate(targets.items()):
        for s, smp in enumerate(tgt["samples"]):
            for rule in rules:
                if rule in done:
                    continue
                e_reps = [energy(stats(generate(smp["D"], m, *rule, rng), smp["D"]), smp["stats"])
                          for _ in range(args.n_reps)]
                rows.append({"target": name, "sample": s, "half": "fit" if s < half else "test",
                             "eta": rule[0], "alpha": rule[1], "beta": rule[2],
                             "energy": float(np.mean(e_reps))})
        print(f"[spatial] simulated all rules on {name} ({t + 1}/{len(targets)})", flush=True)
    E = pd.DataFrame(rows)
    if done:
        E = pd.concat([old, E], ignore_index=True)
    E.to_csv(energy_path, index=False)

    rule_cols = ["eta", "alpha", "beta"]
    test = E[E.half == "test"].groupby(["target"] + rule_cols)["energy"].mean()

    def best_rules(fit_rows):
        g = fit_rows.groupby(["target"] + rule_cols)["energy"].mean().reset_index()
        b = g.loc[g.groupby("target")["energy"].idxmin()].set_index("target")
        return {t: tuple(b.loc[t, rule_cols]) for t in b.index}

    def transfer(best):
        out = []
        for src in targets:
            for tgt in targets:
                e_null = test[(tgt, 0.0, 0.0, 0.0)]
                e_self, e_src = test[(tgt, *best[tgt])], test[(tgt, *best[src])]
                gs, gt = targets[src]["group"], targets[tgt]["group"]
                kind = "self" if src == tgt else ("within " + gs if gs == gt else "from " + gs + " to " + gt)
                out.append({"source": src, "target": tgt, "kind": kind, "energy_null": e_null,
                            "energy_self": e_self, "energy_source": e_src,
                            "TI": (e_null - e_src) / (e_null - e_self) if e_null > e_self else np.nan})
        return pd.DataFrame(out)

    fit = E[E.half == "fit"]
    best = best_rules(fit)
    T = transfer(best)
    T.to_csv(RESULTS_DIR / f"rule_transfer_spatial_{tag}.csv", index=False)

    boots = []
    for b in range(args.n_boot):
        parts = []
        for name in targets:
            chosen = rng.integers(0, half, half)
            f = fit[fit.target == name]
            parts += [f[f["sample"] == c] for c in chosen]
        Tb = transfer(best_rules(pd.concat(parts)))
        Tb["boot"] = b
        boots.append(Tb)
    B = pd.concat(boots)
    B.to_csv(RESULTS_DIR / f"rule_transfer_spatial_boot_{tag}.csv", index=False)

    pd.set_option("display.width", 200)
    print("\n=== Best rule per graph (eta: distance, alpha: homophily, beta: preferential attachment) ===")
    print(pd.DataFrame(best, index=rule_cols).T.assign(group=lambda x: [targets[t]["group"] for t in x.index]).to_string())
    print("\n=== Transfer index by kind: point estimate (median over pairs), bootstrap 95% interval ===")
    ns = T[T.kind != "self"]
    pb = B[B.kind != "self"].groupby(["kind", "boot"])["TI"].median().groupby("kind")
    print(pd.DataFrame({"TI": ns.groupby("kind")["TI"].median(), "lo": pb.quantile(0.025),
                        "hi": pb.quantile(0.975), "pairs": ns.groupby("kind").size()}).round(2).to_string())
    print("\nEnergies (test half): self vs null, per target")
    print(T[T.kind == "self"][["target", "energy_self", "energy_null"]].round(3).to_string(index=False))
    # does the spatial term help? best energy with eta fixed at 0 vs free
    no_space = test.reset_index()
    no_space = no_space[no_space.eta == 0.0].groupby("target")["energy"].min()
    with_space = test.groupby(level="target").min()
    print("\nBest test energy without vs with the distance term:")
    print(pd.DataFrame({"eta=0": no_space, "eta free": with_space}).round(3).to_string())
    print(f"\nSaved results/rule_*_spatial_{tag}.csv")


if __name__ == "__main__":
    main()
