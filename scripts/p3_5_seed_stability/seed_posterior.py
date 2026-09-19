"""p3_5 seed stability: seed_posterior_v1.json generator.

重建说明（2026-09-09）：本步此前由未入库的一次性执行生成（results 有产物、仓库无
生成器），全链重跑前补齐此缺口。规则经旧 v22 三种子 run 逐字段反推验证：
  - 每种子取 eval.json 的 pair_direction_heldout 中 direction_r > 0 的边集；
  - 稳定边 = 在 >=2/3 种子为正的边；
  - mean_r = 仅在为正的种子上的 direction_r 均值；
  - jaccard = 三对种子的正边集 Jaccard。
全部输入动态读取（用户 2026-09-09 指令：无旧版硬编码残留）。
"""
import argparse
import itertools
import json
from collections import Counter
from pathlib import Path

A = Path("/public/home/mengxl/dzy/pd_product_assets")


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--ver", default="v23", help="CCWM run version (checkpoints/project/ccwm-<ver>-s<seed>)")
    ap.add_argument("--seeds", default="42,43,44")
    ap.add_argument("--min-count", type=int, default=2, help="正方向种子数门槛（默认>=2；多数派可设 >len(seeds)/2）")
    ap.add_argument("--out", default=str(A / "results/p3_5_seed_stability/seed_posterior_v1.json"))
    args = ap.parse_args()

    seeds = [int(s) for s in args.seeds.split(",")]
    pos_sets, rs = {}, {}
    for s in seeds:
        ev = json.loads((A / f"checkpoints/project/ccwm-{args.ver}-s{s}/eval.json").read_text())
        table = ev["pair_direction_heldout"]
        pos_sets[s] = {(p["blood_state"], p["brain_state"]) for p in table if p["direction_r"] > 0}
        rs[s] = {(p["blood_state"], p["brain_state"]): p["direction_r"] for p in table}

    jaccard = {
        f"s{a}_s{b}": round(len(pos_sets[a] & pos_sets[b]) / len(pos_sets[a] | pos_sets[b]), 3)
        for a, b in itertools.combinations(seeds, 2)
    }
    cnt = Counter()
    for s in seeds:
        for e in pos_sets[s]:
            cnt[e] += 1
    stable_set = {e for e, c in cnt.items() if c >= args.min_count}
    # 输出顺序：按首个种子的 eval 表序（该表已按 -direction_r 排序），其余按出现序
    first_table = json.loads((A / f"checkpoints/project/ccwm-{args.ver}-s{seeds[0]}/eval.json").read_text())["pair_direction_heldout"]
    ordered = [(p["blood_state"], p["brain_state"]) for p in first_table
               if (p["blood_state"], p["brain_state"]) in stable_set]
    ordered += [e for e in stable_set if e not in ordered]
    stable = [
        {
            "blood": b,
            "brain": r,
            "mean_r": round(sum(rs[s][(b, r)] for s in seeds if (b, r) in pos_sets[s])
                            / sum(1 for s in seeds if (b, r) in pos_sets[s]), 4),
        }
        for b, r in ordered
    ]

    out = {"stable_edges": stable, "jaccard": jaccard,
           "rule": f"direction_r>0 in >={args.min_count}/{len(seeds)} seeds; mean over positive seeds",
           "ver": args.ver, "seeds": seeds}
    Path(args.out).parent.mkdir(parents=True, exist_ok=True)
    Path(args.out).write_text(json.dumps(out, indent=2, ensure_ascii=False))
    print(json.dumps({"n_stable": len(stable), "jaccard": jaccard}, ensure_ascii=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
