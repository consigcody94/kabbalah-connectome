"""
Multi-database literature search across OpenAlex, Semantic Scholar, PubMed,
Crossref, Europe PMC.

Six topics x ~3 APIs each. Filters results for relevance.
"""
import json
import os
import sys
import time
import urllib.parse
import urllib.request

sys.stdout.reconfigure(encoding="utf-8")

EMAIL = "rockvillecodymaryland@gmail.com"
HDRS = {"User-Agent": f"academic-search-skill ({EMAIL})"}
OUT = "data/multi_db_results.json"

TOPICS = [
    ("1_kabbalah_brain", [
        '"Kabbalah" AND ("brain" OR "neuroscience" OR "neural")',
        'sephirot AND brain',
        'Hermetic AND cognitive science',
    ]),
    ("2_jung_kabbalah", [
        'Drob AND Kabbalah',
        'Jung AND Kabbalah',
        'Jungian AND sephirot',
    ]),
    ("3_mcgilchrist_followup", [
        'McGilchrist AND hemispheric',
        '"Master and his Emissary"',
        'hemispheric AND attention AND lateralization 2022',
    ]),
    ("4_callosum_graph", [
        'corpus callosum graph theory bottleneck connectome',
        'corpus callosum betweenness centrality brain network',
        'callosotomy connectome graph',
    ]),
    ("5_esoteric_graphs", [
        'I Ching network graph',
        'Tarot graph theory',
        'mandala network analysis',
    ]),
    ("6_sacred_geo_brain", [
        '"Flower of Life" brain',
        'sacred geometry neuroanatomy',
        'Phi golden ratio brain',
    ]),
]


def fetch(url, timeout=30):
    req = urllib.request.Request(url, headers=HDRS)
    try:
        with urllib.request.urlopen(req, timeout=timeout) as r:
            return r.read().decode("utf-8", errors="replace")
    except Exception as e:
        return f"__ERR__:{e}"


def openalex(query, n=10):
    url = (f"https://api.openalex.org/works?"
           f"search={urllib.parse.quote(query)}"
           f"&per_page={n}&mailto={EMAIL}")
    raw = fetch(url)
    if raw.startswith("__ERR__"):
        return [], raw
    try:
        d = json.loads(raw)
        out = []
        for w in d.get("results", []):
            title = w.get("title") or w.get("display_name") or ""
            year = w.get("publication_year")
            doi = w.get("doi")
            authors = [a.get("author", {}).get("display_name", "")
                       for a in w.get("authorships", [])][:3]
            cites = w.get("cited_by_count", 0)
            host = ((w.get("primary_location") or {}).get("source") or {}).get(
                "display_name", "")
            abs_idx = w.get("abstract_inverted_index") or {}
            abstract = ""
            if abs_idx:
                positions = []
                for word, posns in abs_idx.items():
                    for p in posns:
                        positions.append((p, word))
                positions.sort()
                abstract = " ".join(w for _, w in positions)[:600]
            out.append({
                "src": "OpenAlex", "title": title, "year": year,
                "authors": authors, "venue": host, "doi": doi,
                "cites": cites, "abstract": abstract,
            })
        return out, None
    except Exception as e:
        return [], f"parse error: {e}"


def semantic_scholar(query, n=10):
    url = (f"https://api.semanticscholar.org/graph/v1/paper/search?"
           f"query={urllib.parse.quote(query)}"
           f"&limit={n}"
           f"&fields=title,authors,year,abstract,citationCount,"
           f"venue,externalIds,openAccessPdf,url,tldr")
    raw = fetch(url)
    if raw.startswith("__ERR__"):
        return [], raw
    try:
        d = json.loads(raw)
        out = []
        for w in d.get("data") or []:
            doi = (w.get("externalIds") or {}).get("DOI")
            tldr = (w.get("tldr") or {}).get("text") or ""
            out.append({
                "src": "S2", "title": w.get("title"), "year": w.get("year"),
                "authors": [a.get("name") for a in (w.get("authors") or [])][:3],
                "venue": w.get("venue"), "doi": doi,
                "cites": w.get("citationCount", 0),
                "abstract": (w.get("abstract") or "")[:600],
                "tldr": tldr,
            })
        return out, None
    except Exception as e:
        return [], f"parse error: {e}"


def pubmed(query, n=10):
    s_url = (f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?"
             f"db=pubmed&term={urllib.parse.quote(query)}&retmax={n}&retmode=json")
    raw = fetch(s_url)
    if raw.startswith("__ERR__"):
        return [], raw
    try:
        d = json.loads(raw)
        ids = (d.get("esearchresult") or {}).get("idlist") or []
        if not ids:
            return [], None
        time.sleep(0.5)
        sm_url = (f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?"
                  f"db=pubmed&id={','.join(ids)}&retmode=json")
        sm_raw = fetch(sm_url)
        sm = json.loads(sm_raw)
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
                "src": "PubMed", "title": r.get("title"),
                "year": year, "authors": authors,
                "venue": r.get("fulljournalname") or r.get("source"),
                "doi": doi, "cites": None, "abstract": "",
            })
        return out, None
    except Exception as e:
        return [], f"parse error: {e}"


def relevant(r, kws):
    """Loose relevance check: title/abstract contains at least one keyword
    from the relevance list."""
    text = f"{r.get('title','')} {r.get('abstract','')}".lower()
    return any(k.lower() in text for k in kws)


# Relevance keywords per topic — used to filter spam
RELEVANCE = {
    "1_kabbalah_brain":      ["kabbalah", "sephirot", "qabal", "hermetic",
                              "mystical", "esoteric"],
    "2_jung_kabbalah":       ["kabbalah", "sephirot", "jung", "drob"],
    "3_mcgilchrist_followup":["mcgilchrist", "hemispher", "lateral",
                              "divided brain", "right hemisphere",
                              "left hemisphere"],
    "4_callosum_graph":      ["callosum", "interhemispheric", "commissur",
                              "callosotomy"],
    "5_esoteric_graphs":     ["i ching", "yi jing", "hexagram", "tarot",
                              "kabbalah", "mandala", "esoteric", "kinship"],
    "6_sacred_geo_brain":    ["flower of life", "sacred geometry", "golden ratio",
                              "phi", "fibonacci", "platonic", "vesica"],
}


results = {}
if os.path.exists(OUT):
    try:
        results = json.load(open(OUT))
    except Exception:
        results = {}

for topic, queries in TOPICS:
    if topic in results and results[topic]:
        print(f"\n=== {topic}: cached ({len(results[topic])}) ===", flush=True)
        continue
    print(f"\n=== {topic} ===", flush=True)
    rel_kws = RELEVANCE[topic]
    raw_pool = []
    for q in queries:
        print(f"  > query: {q}", flush=True)
        for fn_name, fn in [("OpenAlex", openalex),
                            ("S2",        semantic_scholar),
                            ("PubMed",    pubmed)]:
            r, err = fn(q, n=10)
            if err:
                print(f"    {fn_name}: ERR {err[:80]}", flush=True)
            else:
                print(f"    {fn_name}: {len(r)} raw -> ", end="")
                rel = [x for x in r if relevant(x, rel_kws)]
                print(f"{len(rel)} relevant", flush=True)
                raw_pool.extend(rel)
            time.sleep(1.0)
        time.sleep(2.0)
    # dedupe by DOI or title
    seen = set()
    final = []
    for r in raw_pool:
        key = (r.get("doi") or r.get("title") or "").lower()[:80]
        if key in seen or not key: continue
        seen.add(key)
        final.append(r)
    # sort by citations desc
    final.sort(key=lambda x: -(x.get("cites") or 0))
    print(f"  final relevant unique: {len(final)}", flush=True)
    for r in final[:5]:
        print(f"   [{r.get('year')}] cites={r.get('cites')} {r.get('title','')[:90]}",
              flush=True)
    results[topic] = final
    json.dump(results, open(OUT, "w"), indent=2, ensure_ascii=False)

print(f"\nDONE. Wrote {OUT}.")
