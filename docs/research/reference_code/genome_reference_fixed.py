#!/usr/bin/env python3
"""Resolve a reference genome name, accession or file path to an indexed FASTA.

No API key and no contact detail is sent to NCBI. Only the standard library
plus pysam is required (pandas/openpyxl only for XLSX position files).

Public API
    resolve_genome(name_or_path, cache_dir=None, download=True, ...) -> dict
    parse_positions(source, genome=None, ...) -> list[dict]
    fetch(genome, seqid, start, end) -> str          (1-based inclusive)

Resolution order
    1. existing local file (plain, gzip or bgzip FASTA)
    2. assembly accession      GCF_/GCA_ + 9 digits (+ .version)
    3. nucleotide accession    NC_000962.3, NZ_CP012345.1, MN908947.3, AL123456
    4. organism name -> NCBI Taxonomy -> RefSeq reference genome

Version 2 (2026-11-11). Changes against version 1 (genome_reference_v1.py):
    F1  organism names containing "/" (influenza strains) are no longer taken for file paths
    F2  rows with stray numeric tokens (thousands separators such as 1,234,567 or
        NC_000962.3:1,234) raise an error instead of giving a silently wrong target
    F3  targets wider than max_span (default 50,000 bp) raise an error (numeric name + position
        used to become a megabase range)
    F4  header tables with both pos and end columns honour end
    F5  unreadable, empty or duplicate-id FASTA files raise GenomeResolutionError
    F6  per-query source and notes are kept in queries.json, so a cache hit no longer shows the
        source and notes of whichever query last wrote metadata.json
    F7  fetch() checks bounds (no silent truncation) and accepts aliases
    F8  BED track/browser lines are skipped; .bed rows are read strictly as chrom, start, end, name
    F9  list items separated by spaces are split like file rows
    F10 "ref" is no longer a seqid column name (it collides with the REF allele column)
    F11 downloads take an exclusive lock file (two processes on one cache)
    F12 scaffold counts in the warning are reported by NCBI role (alt-scaffold, fix-patch, ...)
    F13 a missing position file (name ending .txt/.tsv/.csv/.xlsx/.bed) raises FileNotFoundError
"""
from __future__ import annotations

import contextlib
import gzip
import hashlib
import json
import os
import re
import shutil
import sys
import tempfile
import threading
try:
    import fcntl
except ImportError:          # Windows
    fcntl = None
import time
import urllib.error
import urllib.parse
import urllib.request
import zipfile
from pathlib import Path

import pysam

DATASETS = "https://api.ncbi.nlm.nih.gov/datasets/v2"
EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
USER_AGENT = "primerjinn-genome-resolver/1.0"
NAME_TTL_S = 30 * 24 * 3600          # name -> accession cache lifetime
ASSEMBLY_RE = re.compile(r"^GC[AF]_\d{9}(\.\d+)?$", re.I)
# RefSeq (NC_, NZ_, NG_, NT_, NW_, AC_ ...) and INSDC (U00096, MN908947, CP012345, AL123456)
NUCCORE_RE = re.compile(
    r"^(?:[A-Z]{2}_[A-Z]{0,6}\d{5,}|[A-Z]\d{5}|[A-Z]{2}\d{6}|[A-Z]{2}\d{8}"
    r"|[A-Z]{4}\d{8,10}|[A-Z]{6}\d{9,11})(\.\d+)?$", re.I)   # last two: INSDC WGS/TSA contigs


class GenomeResolutionError(Exception):
    """Could not resolve the query. The message says what to try next."""


class AmbiguousGenomeError(GenomeResolutionError):
    """Several genomes match. `.alternatives` lists them."""

    def __init__(self, msg, alternatives):
        super().__init__(msg)
        self.alternatives = alternatives


# ----------------------------------------------------------------------------
# HTTP with per-host rate limiting and retry
# ----------------------------------------------------------------------------
class _Limiter:
    def __init__(self, per_second):
        self.gap = 1.0 / per_second
        self.t = 0.0
        self.lock = threading.Lock()

    def wait(self):
        with self.lock:
            now = time.monotonic()
            delay = self.t + self.gap - now
            if delay > 0:
                time.sleep(delay)
            self.t = time.monotonic()


# Keyless limits: Datasets 5 req/s, E-utilities 3 req/s. Stay under both.
_LIMITERS = {"api.ncbi.nlm.nih.gov": _Limiter(4.0), "eutils.ncbi.nlm.nih.gov": _Limiter(2.5)}


def _open(url, timeout=60, retries=5):
    host = urllib.parse.urlsplit(url).netloc
    lim = _LIMITERS.get(host)
    last = None
    for attempt in range(retries):
        if lim:
            lim.wait()
        req = urllib.request.Request(url, headers={"User-Agent": USER_AGENT})
        try:
            return urllib.request.urlopen(req, timeout=timeout)
        except urllib.error.HTTPError as e:
            last = e
            if (e.code == 429 or e.code >= 500) and attempt < retries - 1:
                retry_after = e.headers.get("Retry-After")
                time.sleep(float(retry_after) if retry_after and retry_after.isdigit()
                           else 1.0 * 2 ** attempt)
                continue
            if e.code == 429 or e.code >= 500:
                break
            raise
        except (urllib.error.URLError, TimeoutError, ConnectionError) as e:
            last = e
            if attempt < retries - 1:
                time.sleep(1.0 * 2 ** attempt)
    raise GenomeResolutionError(f"NCBI request failed after {retries} tries: {url} ({last})")


def _get_json(url):
    with _open(url) as r:
        return json.loads(r.read().decode() or "{}")


def _download(url, dest, timeout=600):
    tmp = f"{dest}.part{os.getpid()}"
    with _open(url, timeout=timeout) as r, open(tmp, "wb") as fh:
        shutil.copyfileobj(r, fh, 1 << 20)
    os.replace(tmp, dest)


@contextlib.contextmanager
def _locked(gdir):
    """Exclusive lock per genome folder so two processes do not download the same genome twice."""
    Path(gdir).mkdir(parents=True, exist_ok=True)
    fh = open(Path(gdir) / ".lock", "w")
    try:
        if fcntl:
            fcntl.flock(fh, fcntl.LOCK_EX)
        yield
    finally:
        if fcntl:
            fcntl.flock(fh, fcntl.LOCK_UN)
        fh.close()


def _q(s):
    return urllib.parse.quote(str(s), safe="")


# ----------------------------------------------------------------------------
# FASTA helpers
# ----------------------------------------------------------------------------
def _is_gzip(path):
    with open(path, "rb") as fh:
        return fh.read(2) == b"\x1f\x8b"


def _write_upper_fasta(lines, dest):
    """Copy FASTA text to dest, upper-casing sequence lines (removes soft-masking)."""
    tmp = f"{dest}.part{os.getpid()}"
    n = 0
    with open(tmp, "wb") as out:
        for line in lines:
            if line.startswith(b">"):
                n += 1
                out.write(line)
            else:
                out.write(line.upper())
    if n == 0:
        os.remove(tmp)
        raise GenomeResolutionError("downloaded data contains no FASTA record")
    os.replace(tmp, dest)
    return n


def _index(fasta, cache_dir):
    """Return the .fai path, creating it next to the FASTA or, if read-only, in the cache."""
    fai = Path(str(fasta) + ".fai")
    if fai.exists() and fai.stat().st_mtime >= Path(fasta).stat().st_mtime:
        return str(fai)
    _check_fasta(fasta)
    writable = os.access(Path(fasta).parent, os.W_OK)
    try:
        if writable:
            pysam.faidx(str(fasta))
            return str(fai)
        alt = Path(cache_dir) / "local_index" / (_file_key(fasta) + ".fai")
        alt.parent.mkdir(parents=True, exist_ok=True)
        if not alt.exists():
            pysam.faidx(str(fasta), "--fai-idx", str(alt))      # index kept outside a read-only folder
        return str(alt)
    except Exception as e:
        raise GenomeResolutionError(
            f"{fasta} cannot be indexed as FASTA (lines of unequal length inside a record, or not "
            f"a FASTA file): {str(e).strip()[:200]}") from e


def _check_fasta(fasta):
    """Reject files that are not FASTA or that repeat a sequence id (faidx would drop the repeat)."""
    seen, first = set(), None
    with open(fasta, "rb") as fh:
        for line in fh:
            if first is None and line.strip():
                first = line
                if not line.startswith(b">"):
                    break
            if line.startswith(b">"):
                sid = line[1:].split(None, 1)[0] if line[1:].strip() else b""
                if not sid:
                    raise GenomeResolutionError(f"{fasta}: FASTA header without a sequence id.")
                if sid in seen:
                    raise GenomeResolutionError(f"{fasta}: duplicate sequence id '{sid.decode()}'.")
                seen.add(sid)
    if first is None or not first.startswith(b">"):
        raise GenomeResolutionError(f"{fasta} is empty or does not start with a FASTA header ('>').")


def _file_key(path):
    st = Path(path).stat()
    h = hashlib.sha1(f"{Path(path).resolve()}|{st.st_size}|{st.st_mtime_ns}".encode()).hexdigest()[:12]
    return f"{Path(path).name.split('.')[0]}-{h}"


def _gc_scan(fasta, fai):
    """Per-sequence length and GC from the file itself. GC = (G+C)/(A+C+G+T), N and IUPAC excluded."""
    out = []
    with pysam.FastaFile(str(fasta), filepath_index=fai) as fa:
        for sid, length in zip(fa.references, fa.lengths):
            gc = at = 0
            for i in range(0, length, 50_000_000):          # chunked: bounded memory on large genomes
                s = fa.fetch(sid, i, min(length, i + 50_000_000)).upper()
                gc += s.count("G") + s.count("C")
                at += s.count("A") + s.count("T")
            out.append({"seqid": sid, "length": length, "acgt": gc + at,
                        "gc_percent": round(100.0 * gc / (gc + at), 2) if gc + at else None,
                        "_gc": gc})
    return out


def _descriptions(fasta):
    d = {}
    opener = gzip.open if _is_gzip(fasta) else open
    with opener(fasta, "rb") as fh:
        for line in fh:
            if line.startswith(b">"):
                parts = line[1:].decode(errors="replace").strip().split(None, 1)
                d[parts[0]] = parts[1] if len(parts) > 1 else ""
    return d


def _finalise(meta, fasta, cache_dir, seq_info=None):
    """Index, scan GC, pick the default sequence, fill meta."""
    fai = _index(fasta, cache_dir)
    scan = _gc_scan(fasta, fai)
    desc = _descriptions(fasta) if seq_info is None else {}
    seq_info = seq_info or {}
    seqs = []
    for s in scan:
        info = seq_info.get(s["seqid"], {})
        d = desc.get(s["seqid"], "")
        mtype = info.get("type") or (
            "Plasmid" if re.search(r"\bplasmid\b", d, re.I) else
            "Organelle" if re.search(r"mitochondri|chloroplast|apicoplast|plastid", d, re.I) else
            "Chromosome")
        seqs.append({"seqid": s["seqid"], "length": s["length"], "gc_percent": s["gc_percent"],
                     "type": mtype, "role": info.get("role", "unknown"),
                     "name": info.get("name", ""), "aliases": _aliases(s["seqid"], info)})
    tot_gc = sum(s["_gc"] for s in scan)
    tot = sum(s["acgt"] for s in scan)
    meta.update(fasta=str(fasta), fai=fai, sequences=seqs,
                n_sequences=len(seqs), total_length=sum(s["length"] for s in seqs),
                gc_percent=round(100.0 * tot_gc / tot, 2) if tot else None)
    _set_default(meta)
    return meta


def _aliases(seqid, info):
    al = {seqid.split(".")[0]}
    for k in ("genbank", "name"):
        v = info.get(k)
        if v:
            al.add(v)
            al.add(v.split(".")[0])
            if k == "name" and info.get("type") == "Chromosome" and re.fullmatch(r"\d+|[XYZWxyzw]", v):
                al.add("chr" + v)
    al.discard(seqid)
    return sorted(al)


def _set_default(meta):
    seqs = meta["sequences"]
    primary = [s for s in seqs if s["type"] == "Chromosome" and s["role"] in ("assembled-molecule", "unknown")]
    pool = primary or seqs
    best = max(pool, key=lambda s: s["length"])
    meta["default_seqid"] = best["seqid"]
    meta["default_unambiguous"] = len(seqs) == 1 or len(primary) == 1
    if len(seqs) > 1:
        kinds = {}
        for s in seqs:
            k = s["type"] if s["role"] in ("assembled-molecule", "unknown") else s["role"]
            kinds[k] = kinds.get(k, 0) + 1
        what = ", ".join(f"{n} {k.lower()}" for k, n in sorted(kinds.items()))
        if meta["default_unambiguous"]:
            meta["warnings"].append(
                f"{len(seqs)} sequences ({what}); bare integer positions are taken on the only "
                f"chromosome {best['seqid']} ({best['length']} bp). Use seqid:pos for the others.")
        else:
            meta["warnings"].append(
                f"{len(seqs)} sequences ({what}); bare integer positions are ambiguous. "
                f"Give positions as seqid:pos (longest sequence is {best['seqid']}).")


# ----------------------------------------------------------------------------
# NCBI lookups
# ----------------------------------------------------------------------------
def _taxon_lookup(name):
    """Exact Taxonomy name/common name/tax id -> node dict, or None."""
    d = _get_json(f"{DATASETS}/taxonomy/taxon/{_q(name)}")
    for node in d.get("taxonomy_nodes", []):
        if "taxonomy" in node:
            return node["taxonomy"]
    return None


def _taxon_suggest(name):
    d = _get_json(f"{DATASETS}/taxonomy/taxon_suggest/{_q(name)}")
    return d.get("sci_name_and_ids", [])


def _tokens(s):
    return set(re.findall(r"[a-z0-9]+", s.lower()))


def _assembly_row(r):
    ai, st, org = r.get("assembly_info", {}), r.get("assembly_stats", {}), r.get("organism", {})
    return {"accession": r["accession"], "paired_accession": r.get("paired_accession"),
            "organism": org.get("organism_name"), "tax_id": org.get("tax_id"),
            "infraspecific": org.get("infraspecific_names", {}),
            "assembly_name": ai.get("assembly_name"), "assembly_level": ai.get("assembly_level"),
            "refseq_category": ai.get("refseq_category"), "assembly_status": ai.get("assembly_status"),
            "total_length": int(st.get("total_sequence_length", 0)),
            "gc_percent_ncbi": st.get("gc_percent"),
            "n_chromosomes": st.get("total_number_of_chromosomes"),
            "n_sequences": st.get("number_of_component_sequences")}


def _genome_reports(path, params):
    url = f"{DATASETS}/genome/{path}/dataset_report?" + urllib.parse.urlencode(params)
    return [_assembly_row(r) for r in _get_json(url).get("reports", [])]


def _fmt_alt(rows):
    return "\n".join(
        f"  {r['accession']}  {r['organism']}  [{r['assembly_name']}; {r['assembly_level']}; "
        f"{r['refseq_category'] or 'not reference'}; {r['total_length']} bp]" for r in rows)


def _name_to_assembly(name, notes):
    tax = _taxon_lookup(name)
    if tax is None:
        sugg = _taxon_suggest(name)
        qt = _tokens(name)
        hits = [s for s in sugg if qt <= _tokens(s["sci_name"]) or qt <= _tokens(s.get("matched_term", ""))]
        if len(hits) == 1:
            tax = _taxon_lookup(hits[0]["tax_id"])
            notes.append(f"'{name}' is not an exact NCBI Taxonomy name; matched "
                         f"'{hits[0]['sci_name']}' (taxid {hits[0]['tax_id']}) through taxon_suggest.")
        else:
            row = _assembly_esearch(name)
            if row:
                notes.append(f"'{name}' matched through an Entrez assembly search.")
                return row
            listing = "\n".join(f"  {s['sci_name']} (taxid {s['tax_id']})" for s in (hits or sugg)[:10])
            if not listing:
                raise GenomeResolutionError(
                    f"'{name}' is not a recognised NCBI Taxonomy name and nothing similar was found. "
                    f"Check the spelling, or pass a taxid, an accession (GCF_..., NC_...) or a FASTA path.")
            raise AmbiguousGenomeError(
                f"'{name}' is not a recognised NCBI Taxonomy name."
                + (f" Closest names:\n{listing}\n" if listing else " ")
                + "Pass one of these names, a taxid, or an accession (GCF_..., NC_...).",
                (hits or sugg)[:10])
    taxid = tax["tax_id"]
    rows = _genome_reports(f"taxon/{taxid}", {"filters.reference_only": "true", "page_size": 20})
    if len(rows) == 1:
        if rows[0]["tax_id"] != taxid:
            notes.append(f"'{tax['organism_name']}' (taxid {taxid}) resolved to the reference genome of "
                         f"'{rows[0]['organism']}' (taxid {rows[0]['tax_id']}).")
        return rows[0]
    if len(rows) > 1:
        raise AmbiguousGenomeError(
            f"'{name}' (taxid {taxid}) has {len(rows)} RefSeq reference genomes. Pass one accession:\n"
            + _fmt_alt(rows), rows)
    # no designated reference genome under this taxon: RefSeq assemblies, complete first
    for extra in ({"filters.assembly_level": "complete_genome"}, {}):
        rows = _genome_reports(f"taxon/{taxid}", {
            "filters.assembly_source": "refseq", "filters.exclude_atypical": "true",
            "filters.assembly_version": "current", "page_size": 20, **extra})
        if len(rows) == 1:
            notes.append(f"No designated reference genome for taxid {taxid}; using its only RefSeq "
                         f"assembly {rows[0]['accession']}.")
            return rows[0]
        if rows:
            raise AmbiguousGenomeError(
                f"'{name}' (taxid {taxid}) has no designated reference genome and several RefSeq "
                f"assemblies (first {len(rows)} shown). Pass one accession:\n" + _fmt_alt(rows), rows)
    raise GenomeResolutionError(
        f"'{name}' (taxid {taxid}) has no RefSeq assembly. Pass a GenBank assembly (GCA_...) "
        f"or nucleotide accession, or a local FASTA.")


def _assembly_esearch(name):
    term = f'{name} AND "reference genome"[filter] AND latest[filter]'
    d = _get_json(f"{EUTILS}/esearch.fcgi?db=assembly&retmode=json&retmax=5&tool=primerjinn&term={_q(term)}")
    ids = d.get("esearchresult", {}).get("idlist", [])
    if len(ids) != 1:
        return None
    s = _get_json(f"{EUTILS}/esummary.fcgi?db=assembly&retmode=json&tool=primerjinn&id={ids[0]}")
    acc = s.get("result", {}).get(ids[0], {}).get("assemblyaccession")
    rows = _genome_reports(f"accession/{acc}", {}) if acc else []
    return rows[0] if rows else None


def _sequence_reports(acc):
    info, token = {}, None
    while True:
        p = {"page_size": 1000}
        if token:
            p["page_token"] = token
        d = _get_json(f"{DATASETS}/genome/accession/{acc}/sequence_reports?" + urllib.parse.urlencode(p))
        for r in d.get("reports", []):
            sid = r.get("refseq_accession") if acc.upper().startswith("GCF") else r.get("genbank_accession")
            sid = sid or r.get("refseq_accession") or r.get("genbank_accession")
            info[sid] = {"name": r.get("chr_name", ""), "type": r.get("assigned_molecule_location_type", ""),
                         "role": r.get("role", ""), "length": r.get("length"),
                         "gc_percent_ncbi": r.get("gc_percent"),
                         "genbank": r.get("genbank_accession") if sid != r.get("genbank_accession") else None}
        token = d.get("next_page_token")
        if not token:
            return info


def _fetch_assembly(row, cache_dir, max_download_mb):
    acc = row["accession"]
    gdir = Path(cache_dir) / acc
    fasta = gdir / f"{acc}.fna"
    if fasta.exists():
        return fasta
    with _locked(gdir):
        return _fetch_assembly_locked(row, cache_dir, max_download_mb, acc, gdir, fasta)


def _fetch_assembly_locked(row, cache_dir, max_download_mb, acc, gdir, fasta):
    if fasta.exists():                      # another process finished while we waited
        return fasta
    summ = _get_json(f"{DATASETS}/genome/accession/{acc}/download_summary?include_annotation_type=GENOME_FASTA")
    size = summ.get("hydrated", {}).get("estimated_file_size_mb", 0)
    if size > max_download_mb:
        raise GenomeResolutionError(
            f"{acc} download is about {size} MB, above max_download_mb={max_download_mb}. "
            f"Raise the limit or pass a local FASTA.")
    gdir.mkdir(parents=True, exist_ok=True)
    zpath = gdir / "ncbi_dataset.zip"
    _download(f"{DATASETS}/genome/accession/{acc}/download?include_annotation_type=GENOME_FASTA", zpath)
    if not zipfile.is_zipfile(zpath):          # unknown accessions return HTTP 200 and a broken zip
        zpath.unlink()
        raise GenomeResolutionError(f"NCBI Datasets returned no valid archive for {acc}.")
    with zipfile.ZipFile(zpath) as z:
        names = [n for n in z.namelist() if n.endswith("_genomic.fna")]
        if not names:
            raise GenomeResolutionError(f"No genomic FASTA inside the Datasets archive for {acc}.")
        md5 = {}
        if "md5sum.txt" in z.namelist():
            for line in z.read("md5sum.txt").decode().splitlines():
                h, _, n = line.partition(" ")
                md5[n.strip()] = h
        h = hashlib.md5()
        with z.open(names[0]) as fh:
            for chunk in iter(lambda: fh.read(1 << 20), b""):
                h.update(chunk)
        if md5.get(names[0]) and md5[names[0]] != h.hexdigest():
            raise GenomeResolutionError(f"md5 mismatch for {names[0]}; delete {gdir} and retry.")
        with z.open(names[0]) as fh:
            _write_upper_fasta(fh, str(fasta))
    zpath.unlink()
    return fasta


def _nuccore_summary(acc):
    d = _get_json(f"{EUTILS}/esummary.fcgi?db=nuccore&retmode=json&tool=primerjinn&id={_q(acc)}")
    res = d.get("result", {})
    if d.get("error") or not res.get("uids"):
        raise GenomeResolutionError(
            f"'{acc}' looks like a nucleotide accession but NCBI nuccore does not know it "
            f"({d.get('error', 'no record')}).")
    return res[res["uids"][0]]


def _fetch_nuccore(accver, slen, cache_dir):
    gdir = Path(cache_dir) / accver
    fasta = gdir / f"{accver}.fna"
    if fasta.exists():
        return fasta
    with _locked(gdir):
        return _fetch_nuccore_locked(accver, slen, gdir, fasta)


def _fetch_nuccore_locked(accver, slen, gdir, fasta):
    if fasta.exists():
        return fasta
    url = f"{EUTILS}/efetch.fcgi?db=nuccore&rettype=fasta&retmode=text&tool=primerjinn&id={_q(accver)}"
    try:
        with _open(url, timeout=600) as r:
            _write_upper_fasta(r, str(fasta))
    except urllib.error.HTTPError as e:          # unknown id -> HTTP 400
        raise GenomeResolutionError(f"efetch failed for {accver}: HTTP {e.code}") from e
    with pysam.FastaFile(str(fasta)) as fa:
        if fa.nreferences != 1 or (slen and fa.lengths[0] != int(slen)):
            got = fa.lengths[0] if fa.nreferences else 0
            fasta.unlink()
            Path(str(fasta) + ".fai").unlink(missing_ok=True)
            raise GenomeResolutionError(
                f"Truncated download for {accver}: got {got} bp, expected {slen} bp. Retry.")
    return fasta


# ----------------------------------------------------------------------------
# Public: resolve_genome
# ----------------------------------------------------------------------------
def default_cache_dir():
    if os.environ.get("PRIMERJINN_CACHE"):
        return Path(os.environ["PRIMERJINN_CACHE"])
    base = os.environ.get("XDG_CACHE_HOME") or (Path.home() / ".cache")
    return Path(base) / "primerjinn" / "genomes"


def resolve_genome(name_or_path, cache_dir=None, download=True, max_download_mb=1500, refresh=False):
    """Resolve a path, accession or organism name to a FASTA plus metadata.

    Returns a dict: query, source (local|assembly|nuccore|taxonomy), accession, organism,
    tax_id, assembly_name, fasta, fai, sequences[{seqid,length,gc_percent,type,role,name,aliases}],
    n_sequences, total_length, gc_percent, default_seqid, default_unambiguous, notes, warnings.
    With download=False nothing is downloaded: fasta is None and lengths/GC are NCBI's values.
    Raises GenomeResolutionError / AmbiguousGenomeError with an actionable message.
    """
    query = str(name_or_path).strip()
    if not query:
        raise GenomeResolutionError("Empty genome name.")
    cache_dir = Path(cache_dir) if cache_dir else default_cache_dir()
    cache_dir.mkdir(parents=True, exist_ok=True)
    meta = {"query": query, "notes": [], "warnings": []}

    # 1. local file
    p = Path(query).expanduser()
    if p.is_file():
        meta.update(source="local", accession=None, organism=None, tax_id=None, assembly_name=None)
        fasta = p.resolve()
        if _is_gzip(fasta):
            # plain gzip cannot be indexed and BGZF is slower than plain: keep a plain copy in the cache
            out = cache_dir / "local" / (_file_key(fasta) + ".fna")
            if not out.exists():
                out.parent.mkdir(parents=True, exist_ok=True)
                with gzip.open(fasta, "rb") as fh:
                    _write_upper_fasta(fh, str(out))
            meta["notes"].append(f"Decompressed {fasta} to {out}.")
            fasta = out
        return _finalise(meta, fasta, cache_dir)
    looks_like_path = (re.search(r"\.(fa|fna|fasta|fas|fsa)(\.(gz|bgz))?$", query, re.I)
                       or query.startswith(("/", "./", "../", "~", ".\\"))
                       or (os.sep in query and " " not in query and "(" not in query))
    if looks_like_path:                      # F1: 'Influenza A virus (A/Puerto Rico/8/1934(H1N1))' is a name
        raise GenomeResolutionError(f"File not found: {query}")

    # cached metadata for this exact query
    qfile = cache_dir / "queries.json"
    qcache = json.loads(qfile.read_text()) if qfile.exists() else {}
    key = re.sub(r"\s+", " ", query.lower())
    hit = qcache.get(key)
    if hit and not refresh and (hit.get("pinned") or time.time() - hit["time"] < NAME_TTL_S):
        mfile = cache_dir / hit["accession"] / "metadata.json"
        if mfile.exists():
            cached = json.loads(mfile.read_text())
            if cached.get("fasta") and Path(cached["fasta"]).exists():
                cached["query"] = query
                cached["source"] = hit.get("source", cached.get("source"))
                cached["notes"] = hit.get("notes", cached.get("notes", [])) + ["Loaded from cache."]
                return cached

    # 2. assembly accession
    seq_info = None
    if ASSEMBLY_RE.match(query):
        rows = _genome_reports(f"accession/{query.upper()}", {})
        if not rows:
            raise GenomeResolutionError(
                f"Assembly accession {query} not found at NCBI (unknown, suppressed, or a "
                f"superseded version; try it without the .version suffix).")
        row = rows[0]
        if row["accession"].upper() != query.upper():
            meta["notes"].append(f"{query} resolved to current version {row['accession']}.")
        if row["assembly_status"] != "current":
            meta["warnings"].append(f"{row['accession']} has status '{row['assembly_status']}'.")
        meta["source"] = "assembly"
    # 3. nucleotide accession
    elif NUCCORE_RE.match(query):
        s = _nuccore_summary(query)
        accver = s["accessionversion"]
        if s.get("status") in ("suppressed", "replaced") or s.get("replacedby"):
            meta["warnings"].append(f"{accver} status: {s.get('status')}, replaced by {s.get('replacedby')}.")
        if accver.upper() != query.upper():
            meta["notes"].append(f"{query} resolved to current version {accver}.")
        meta.update(source="nuccore", accession=accver, organism=s.get("organism"),
                    tax_id=s.get("taxid"), assembly_name=None, title=s.get("title"),
                    moltype=s.get("moltype"), topology=s.get("topology"))
        if not download:
            meta.update(fasta=None, fai=None, n_sequences=1, total_length=int(s["slen"]), gc_percent=None,
                        sequences=[{"seqid": accver, "length": int(s["slen"]), "gc_percent": None,
                                    "type": "Plasmid" if s.get("genome") == "plasmid" else "Chromosome",
                                    "role": "unknown", "name": "", "aliases": [accver.split(".")[0]]}])
            _set_default(meta)
            return meta
        fasta = _fetch_nuccore(accver, s.get("slen"), cache_dir)
        _finalise(meta, fasta, cache_dir)
        return _store(meta, cache_dir, qfile, qcache, key, pinned="." in query)
    # 4. organism name / taxid
    else:
        row = _name_to_assembly(query, meta["notes"])
        meta["source"] = "taxonomy"

    acc = row["accession"]
    meta.update(accession=acc, organism=row["organism"], tax_id=row["tax_id"],
                assembly_name=row["assembly_name"], assembly_level=row["assembly_level"],
                refseq_category=row["refseq_category"], infraspecific=row["infraspecific"],
                gc_percent_ncbi=row["gc_percent_ncbi"])
    seq_info = _sequence_reports(acc)
    if not download:
        seqs = [{"seqid": sid, "length": i["length"], "gc_percent": i["gc_percent_ncbi"], "type": i["type"],
                 "role": i["role"], "name": i["name"], "aliases": _aliases(sid, i)}
                for sid, i in seq_info.items()]
        meta.update(fasta=None, fai=None, sequences=seqs, n_sequences=len(seqs),
                    total_length=row["total_length"], gc_percent=row["gc_percent_ncbi"])
        _set_default(meta)
        return meta
    fasta = _fetch_assembly(row, cache_dir, max_download_mb)
    _finalise(meta, fasta, cache_dir, seq_info)
    return _store(meta, cache_dir, qfile, qcache, key,
                  pinned=meta["source"] == "assembly" and "." in query)


def _store(meta, cache_dir, qfile, qcache, key, pinned):
    mfile = Path(cache_dir) / meta["accession"] / "metadata.json"
    mfile.write_text(json.dumps(meta, indent=1))
    qcache = {**qcache, **(json.loads(qfile.read_text()) if qfile.exists() else {})}   # keep entries other processes added
    qcache[key] = {"accession": meta["accession"], "time": time.time(), "pinned": pinned,
                   "source": meta["source"], "notes": list(meta["notes"])}
    tmp = f"{qfile}.part{os.getpid()}"
    Path(tmp).write_text(json.dumps(qcache, indent=1))
    os.replace(tmp, qfile)
    return meta


def fetch(genome, seqid, start, end):
    """Sequence for 1-based inclusive [start, end], upper case."""
    s = _seq_lookup(genome).get(str(seqid).lower())
    if s is None:
        raise ValueError(f"unknown sequence '{seqid}'")
    if not (1 <= start <= end <= s["length"]):
        raise ValueError(f"{seqid}:{start}-{end} is outside 1..{s['length']}")
    with pysam.FastaFile(genome["fasta"], filepath_index=genome["fai"]) as fa:
        return fa.fetch(s["seqid"], start - 1, end).upper()


# ----------------------------------------------------------------------------
# Public: parse_positions
# ----------------------------------------------------------------------------
_LOCUS_RE = re.compile(r"^(?:(?P<seq>[^\s:]+):)?(?P<a>\d+)(?:\s*(?:-|\.\.)\s*(?P<b>\d+))?$")
_COLS = {"name": {"name", "id", "target", "label", "gene", "locus", "mutation", "variant"},
         "seqid": {"seqid", "chrom", "chr", "chromosome", "contig", "sequence", "accession", "reference"},
         "pos": {"pos", "position", "site", "coordinate", "coord"},
         "start": {"start", "begin", "from", "chromstart"},
         "end": {"end", "stop", "to", "chromend"}}


def _seq_lookup(genome):
    look = {}
    if genome:
        for s in genome["sequences"]:
            for k in [s["seqid"], *s["aliases"]]:
                look.setdefault(k.lower(), s)
    return look


def _rows_from_file(path):
    path = Path(path)
    ext = path.suffix.lower()
    if ext in (".xlsx", ".xls"):
        import pandas as pd
        df = pd.read_excel(path, header=None, dtype=str)
        return [[("" if v is None or v != v else str(v).strip()) for v in row] for row in df.values.tolist()]
    rows = []
    opener = gzip.open if _is_gzip(path) else open
    with opener(path, "rt", encoding="utf-8-sig") as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith(("#", "track ", "browser ")):
                continue
            delim = "\t" if "\t" in line else "," if "," in line else ";" if ";" in line else None
            rows.append([f.strip().strip('"') for f in (line.split(delim) if delim else line.split())])
    return rows


def parse_positions(source, genome=None, zero_based_half_open=None, allow_ambiguous_default=False,
                    max_span=50_000):
    """Parse targets into [{name, seqid, start, end}] with 1-based inclusive coordinates.

    source: a file (txt/tsv/csv/xlsx/bed, optionally gz) or an iterable of ints/strings.
    Accepted per row: 6735 | NC_000962.3:6735 | 761007-761277 | NC_000962.3:761007-761277,
    optionally with a name, or columns named name / seqid|chrom / pos / start / end.
    .bed files are 0-based half-open and converted. Each row is validated against `genome`.
    """
    is_file = isinstance(source, (str, os.PathLike)) and Path(source).is_file()
    if is_file:
        rows = _rows_from_file(source)
        if zero_based_half_open is None:
            zero_based_half_open = str(source).lower().endswith((".bed", ".bed.gz"))
    else:
        if isinstance(source, (str, os.PathLike)) and re.search(
                r"\.(txt|tsv|csv|xlsx|xls|bed)(\.gz)?$", str(source), re.I):
            raise FileNotFoundError(f"Position file not found: {source}")   # F13
        if isinstance(source, (str, int)):
            source = [source]
        def _split(x):
            x = str(x).strip()
            return re.split(r"[\t,]", x) if re.search(r"[\t,]", x) else re.split(r"\s+(?![-.])(?<![-.]\s)", x)
        rows = [_split(x) if not isinstance(x, (list, tuple)) else [str(v) for v in x] for x in source]
    look = _seq_lookup(genome)
    cols = None
    if rows:
        low = [f.lower() for f in rows[0]]
        known = set().union(*_COLS.values())
        if any(f in known for f in low) and not any(_LOCUS_RE.match(f) for f in rows[0]):
            cols = {}
            for i, f in enumerate(low):
                for k, names in _COLS.items():
                    if f in names:
                        cols.setdefault(k, i)
            rows = rows[1:]
            if "pos" not in cols and "start" not in cols:
                raise ValueError(f"Header {low} has no pos/start column.")
    out, seen = [], set()
    for n, row in enumerate(rows, 1):
        row = [f.strip() for f in row]
        if not any(row):
            continue
        name = seq = None
        try:
            if cols:
                g = lambda k: row[cols[k]] if k in cols and cols[k] < len(row) and row[cols[k]] else None
                name, seq = g("name"), g("seqid")
                a = g("pos") or g("start")
                b = g("end")
                m = _LOCUS_RE.match(a or "")
                if not m:
                    raise ValueError(f"cannot read a position from '{a}'")
                seq = seq or m["seq"]
                if b and not b.isdigit():
                    raise ValueError(f"cannot read an end position from '{b}'")
                if b and m["b"] and int(b) != int(m["b"]):
                    raise ValueError(f"range '{a}' and end column '{b}' disagree")
                start, end = int(m["a"]), int(b or m["b"] or m["a"])
            elif zero_based_half_open:           # F8: BED is positional: chrom, start, end[, name ...]
                if len(row) < 3 or not (row[1].isdigit() and row[2].isdigit()):
                    raise ValueError("a BED row needs chrom, start, end")
                seq, start, end = row[0], int(row[1]), int(row[2])
                name = row[3] if len(row) > 3 and row[3] else None
            else:
                ints = [f for f in row if f.isdigit()]
                loci = [f for f in row if not f.isdigit() and _LOCUS_RE.match(f)]
                rest = [f for f in row if f and not f.isdigit() and f not in loci]
                for f in list(rest):
                    if f.lower() in look:
                        seq = f
                        rest.remove(f)
                        break
                if len(loci) > 1 or (loci and ints) or len(ints) > 2:
                    raise ValueError("more numbers than one position or one start and end "
                                     "(thousands separators such as 761,155 are not supported)")
                if loci:
                    m = _LOCUS_RE.match(loci[0])
                    seq = m["seq"] or seq
                    start, end = int(m["a"]), int(m["b"] or m["a"])
                elif ints:
                    start, end = int(ints[0]), int(ints[1] if len(ints) > 1 else ints[0])
                else:
                    raise ValueError("no position found")
                name = rest[0] if rest else None
            if zero_based_half_open:
                start += 1
            if start < 1 or end < start:
                raise ValueError(f"start {start} and end {end} must satisfy 1 <= start <= end")
            if max_span and end - start + 1 > max_span:
                raise ValueError(f"target {start}-{end} spans {end - start + 1} bp, above max_span={max_span}; "
                                 f"if the first number is a name, put a letter in it")
            if genome:
                if seq is None:
                    if not genome["default_unambiguous"] and not allow_ambiguous_default:
                        raise ValueError(
                            f"the genome has {genome['n_sequences']} sequences, so the position needs a "
                            f"sequence id, for example {genome['default_seqid']}:{start}")
                    s = look[genome["default_seqid"].lower()]
                else:
                    s = look.get(seq.lower())
                    if s is None:
                        ids = ", ".join(x["seqid"] for x in genome["sequences"][:8])
                        raise ValueError(f"unknown sequence '{seq}'. Known: {ids}")
                seq = s["seqid"]
                if end > s["length"]:
                    raise ValueError(f"{end} is beyond the end of {seq} ({s['length']} bp)")
        except ValueError as e:
            raise ValueError(f"Position row {n} {row}: {e}") from None
        k = (seq, start, end)
        if k in seen:
            continue
        seen.add(k)
        out.append({"name": name or (f"{seq}_{start}" if start == end else f"{seq}_{start}_{end}")
                    if seq else name or (str(start) if start == end else f"{start}_{end}"),
                    "seqid": seq, "start": start, "end": end})
    if not out:
        raise ValueError("No positions found in the input.")
    return out


if __name__ == "__main__":
    import argparse
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("genome")
    ap.add_argument("--cache_dir")
    ap.add_argument("--no-download", action="store_true")
    ap.add_argument("--positions")
    a = ap.parse_args()
    try:
        g = resolve_genome(a.genome, a.cache_dir, download=not a.no_download)
    except GenomeResolutionError as e:
        sys.exit(f"ERROR: {e}")
    show = dict(g)
    if len(show["sequences"]) > 6:
        show["sequences"] = show["sequences"][:6] + [f"... {len(g['sequences']) - 6} more"]
    print(json.dumps(show, indent=1))
    if a.positions:
        for t in parse_positions(a.positions, g):
            print(t)
