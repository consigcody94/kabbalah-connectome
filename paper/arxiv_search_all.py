"""Run remaining arxiv searches; save incrementally; UTF-8 safe."""
import json
import os
import sys
import time

# Force UTF-8 stdout on Windows
sys.stdout.reconfigure(encoding="utf-8")

import arxiv

QUERIES = [
    ("1_kabbalah_brain", 'Kabbalah brain'),
    ("1b_sephirot",      'sephirot OR sefirot OR Qliphoth'),
    ("1c_sacred_geo",    'sacred geometry brain'),
    ("2_corpus_callosum",'corpus callosum connectome graph'),
    ("3_bridge_node",    'bridge node centrality brain network'),
    ("4_mythology",      'I Ching graph network structure'),
    ("4b_kinship",       'kinship network graph theory'),
    ("5_hemispheric",    'hemispheric lateralization functional connectivity'),
    ("5b_mcgilchrist",   'left right hemisphere asymmetry attention'),
]

OUT = "data/arxiv_results.json"
results = {}
if os.path.exists(OUT):
    try:
        results = json.load(open(OUT))
    except Exception:
        results = {}

def safe(s):
    return s.encode("ascii", "replace").decode("ascii")[:100]

def run(name, query, n=8, retries=3):
    for attempt in range(retries):
        try:
            client = arxiv.Client(page_size=n, delay_seconds=8, num_retries=3)
            search = arxiv.Search(
                query=query, max_results=n,
                sort_by=arxiv.SortCriterion.Relevance,
            )
            out = []
            for r in client.results(search):
                out.append({
                    "id": r.get_short_id(),
                    "title": r.title,
                    "authors": [a.name for a in r.authors][:3],
                    "year": r.published.year if r.published else None,
                    "category": r.primary_category,
                    "abstract": r.summary[:500],
                })
            return out
        except Exception as e:
            print(f"  attempt {attempt+1}/{retries} FAILED: {e}", flush=True)
            time.sleep(30)
    return []

for name, q in QUERIES:
    if name in results and results[name]:
        print(f"=== {name}: cached ({len(results[name])}) ===", flush=True)
        continue
    print(f"\n=== {name}: '{q}' ===", flush=True)
    out = run(name, q)
    print(f"  -> {len(out)} results", flush=True)
    for r in out[:5]:
        print(f"   [{r['year']}] {safe(r['title'])}", flush=True)
    results[name] = out
    json.dump(results, open(OUT, "w"), indent=2, ensure_ascii=False)
    time.sleep(15)

print(f"\nDone. {sum(len(v) for v in results.values())} total hits across {len(results)} queries.")
