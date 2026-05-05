"""Pull PubMed citations for the brain regions identified in our mapping."""
import json
import sys
import time
import urllib.parse
import urllib.request

sys.stdout.reconfigure(encoding="utf-8")

QUERIES = {
    "thalamus_central_integrator": [
        '"thalamus" AND "central hub" AND "connectome"',
        '"thalamus" AND "rich club" AND "brain network"',
    ],
    "amygdala_bilateral_integration": [
        '"amygdala" AND "anterior commissure" AND "interhemispheric"',
        '"amygdala" AND "bilateral" AND "limbic integration"',
    ],
    "precuneus_dmn_hub": [
        '"precuneus" AND "default mode network" AND "hub"',
        '"precuneus" AND "self-referential"',
    ],
    "thalamic_hub_consciousness": [
        '"thalamus" AND "consciousness" AND "integration"',
        '"thalamic" AND "default mode"',
    ],
    "rich_club_brain": [
        '"rich club" AND "human connectome"',
    ],
    "subcortical_hubs_topology": [
        '"basal ganglia" AND "graph theory" AND "connectome"',
        '"subcortical" AND "hub" AND "structural connectome"',
    ],
}


def fetch(url, timeout=30):
    req = urllib.request.Request(
        url, headers={"User-Agent": "tree-brain-paper/1.0"})
    return urllib.request.urlopen(req, timeout=timeout).read().decode("utf-8")


def pubmed(query, n=5):
    s_url = (f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?"
             f"db=pubmed&term={urllib.parse.quote(query)}&retmax={n}&retmode=json")
    try:
        d = json.loads(fetch(s_url))
    except Exception as e:
        return [], f"search err: {e}"
    ids = (d.get("esearchresult") or {}).get("idlist") or []
    if not ids:
        return [], None
    time.sleep(0.4)
    sm_url = (f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?"
              f"db=pubmed&id={','.join(ids)}&retmode=json")
    try:
        sm = json.loads(fetch(sm_url))
    except Exception as e:
        return [], f"summ err: {e}"
    out = []
    for pid in ids:
        r = (sm.get("result") or {}).get(pid)
        if not r: continue
        authors = [a.get("name") for a in (r.get("authors") or [])][:3]
        year = (r.get("pubdate") or "")[:4]
        doi = ""
        for aid in r.get("articleids") or []:
            if aid.get("idtype") == "doi":
                doi = aid.get("value")
                break
        out.append({
            "pmid": pid,
            "title": r.get("title", ""),
            "year": year,
            "authors": authors,
            "venue": r.get("fulljournalname") or r.get("source"),
            "doi": doi,
        })
    return out, None


results = {}
for topic, queries in QUERIES.items():
    print(f"\n=== {topic} ===")
    bucket = []
    for q in queries:
        print(f"  > {q}")
        out, err = pubmed(q, n=5)
        if err:
            print(f"    ERR: {err}")
            continue
        print(f"    {len(out)} results")
        for r in out[:3]:
            print(f"     [{r['year']}] {r['title'][:90]}")
        bucket.extend(out)
        time.sleep(1.0)
    # Dedupe
    seen = set()
    uniq = []
    for r in bucket:
        if r["pmid"] in seen: continue
        seen.add(r["pmid"])
        uniq.append(r)
    results[topic] = uniq[:5]

with open("data/pubmed_clinical.json", "w") as f:
    json.dump(results, f, indent=2, ensure_ascii=False)
print(f"\nWrote data/pubmed_clinical.json")
print(f"Total unique results: {sum(len(v) for v in results.values())}")
