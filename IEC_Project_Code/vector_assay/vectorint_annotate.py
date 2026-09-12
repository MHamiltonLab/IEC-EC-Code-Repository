#!/usr/bin/env python3
"""Parse integration SAMs, deduplicate molecules, cluster junctions, and annotate host-gene features."""

import gzip
import re
from collections import defaultdict

import numpy as np
from _assay import GENE_REFERENCE, sam_path, sample_layout

REFGENE = str(GENE_REFERENCE)
JUNCTION_WINDOW = 10

# Retained source-curated cancer-gene set; no resource release was supplied.
CANCER_GENES = {
    "ABL1",
    "ABL2",
    "ABI1",
    "ACSL3",
    "AFF1",
    "AFF3",
    "AFF4",
    "AKAP9",
    "AKT1",
    "AKT2",
    "AKT3",
    "ALK",
    "AMER1",
    "APC",
    "AR",
    "ARAF",
    "ARHGAP26",
    "ARHGEF12",
    "ARID1A",
    "ARID1B",
    "ARID2",
    "ARNT",
    "ASPSCR1",
    "ASXL1",
    "ASXL2",
    "ATF1",
    "ATIC",
    "ATM",
    "ATP1A1",
    "ATP2B3",
    "ATR",
    "ATRX",
    "AXIN1",
    "AXIN2",
    "BAP1",
    "BARD1",
    "BCL2",
    "BCL2L1",
    "BCL2L2",
    "BCL3",
    "BCL6",
    "BCL7A",
    "BCL9",
    "BCL9L",
    "BCL10",
    "BCL11A",
    "BCL11B",
    "BCOR",
    "BCORL1",
    "BCR",
    "BIRC3",
    "BIRC6",
    "BLM",
    "BMPR1A",
    "BRAF",
    "BRCA1",
    "BRCA2",
    "BRD3",
    "BRD4",
    "BRIP1",
    "BTG1",
    "BTK",
    "BUB1B",
    "CACNA1D",
    "CALR",
    "CAMTA1",
    "CANT1",
    "CARD11",
    "CARS",
    "CASP8",
    "CBFA2T3",
    "CBFB",
    "CBL",
    "CBLB",
    "CBLC",
    "CCDC6",
    "CCNB1IP1",
    "CCND1",
    "CCND2",
    "CCND3",
    "CCNE1",
    "CD274",
    "CD28",
    "CD74",
    "CD79A",
    "CD79B",
    "CDC73",
    "CDH1",
    "CDH11",
    "CDK4",
    "CDK6",
    "CDK12",
    "CDKN1B",
    "CDKN2A",
    "CDKN2B",
    "CDKN2C",
    "CDX2",
    "CEBPA",
    "CHCHD7",
    "CHD2",
    "CHD4",
    "CHEK2",
    "CHIC2",
    "CIC",
    "CIITA",
    "CLIP1",
    "CLTC",
    "CLTCL1",
    "CNBP",
    "CNOT3",
    "CNTRL",
    "COL1A1",
    "COL2A1",
    "COX6C",
    "CREB1",
    "CREB3L1",
    "CREB3L2",
    "CREBBP",
    "CRLF2",
    "CRTC1",
    "CRTC3",
    "CSF1R",
    "CSF3R",
    "CTCF",
    "CTNNB1",
    "CTNND2",
    "CUX1",
    "CXCR4",
    "CYLD",
    "DAXX",
    "DCAF12L2",
    "DDB2",
    "DDIT3",
    "DDR2",
    "DDX3X",
    "DDX5",
    "DDX6",
    "DDX10",
    "DEK",
    "DGCR8",
    "DICER1",
    "DNM2",
    "DNMT3A",
    "EBF1",
    "ECT2L",
    "EED",
    "EGFR",
    "EIF3E",
    "EIF4A2",
    "ELF4",
    "ELK4",
    "ELL",
    "ELN",
    "EML4",
    "EP300",
    "EPAS1",
    "EPHA3",
    "EPHA7",
    "EPS15",
    "ERBB2",
    "ERBB3",
    "ERBB4",
    "ERC1",
    "ERCC2",
    "ERCC3",
    "ERCC4",
    "ERCC5",
    "ERG",
    "ESR1",
    "ETNK1",
    "ETV1",
    "ETV4",
    "ETV5",
    "ETV6",
    "EWSR1",
    "EXT1",
    "EXT2",
    "EZH2",
    "EZR",
    "FAM131B",
    "FAM46C",
    "FANCA",
    "FANCC",
    "FANCD2",
    "FANCE",
    "FANCF",
    "FANCG",
    "FAS",
    "FBXO11",
    "FBXW7",
    "FCGR2B",
    "FCRL4",
    "FEN1",
    "FES",
    "FEV",
    "FGFR1",
    "FGFR1OP",
    "FGFR2",
    "FGFR3",
    "FGFR4",
    "FH",
    "FHIT",
    "FIP1L1",
    "FLCN",
    "FLI1",
    "FLT3",
    "FLT4",
    "FNBP1",
    "FOXA1",
    "FOXL2",
    "FOXO1",
    "FOXO3",
    "FOXO4",
    "FOXP1",
    "FSTL3",
    "FUBP1",
    "FUS",
    "GAS7",
    "GATA1",
    "GATA2",
    "GATA3",
    "GMPS",
    "GNA11",
    "GNAQ",
    "GNAS",
    "GOLGA5",
    "GOPC",
    "GPC3",
    "GPHN",
    "GRIN2A",
    "H3F3A",
    "H3F3B",
    "HERPUD1",
    "HEY1",
    "HIP1",
    "HIST1H3B",
    "HIST1H4I",
    "HLF",
    "HMGA1",
    "HMGA2",
    "HNF1A",
    "HNRNPA2B1",
    "HOXA9",
    "HOXA11",
    "HOXA13",
    "HOXC11",
    "HOXC13",
    "HOXD11",
    "HOXD13",
    "HRAS",
    "HSP90AA1",
    "HSP90AB1",
    "ID3",
    "IDH1",
    "IDH2",
    "IGF2BP2",
    "IKBKB",
    "IKZF1",
    "IKZF2",
    "IKZF3",
    "IL2",
    "IL6ST",
    "IL7R",
    "IRF4",
    "ITK",
    "JAK1",
    "JAK2",
    "JAK3",
    "JAZF1",
    "JUN",
    "KAT6A",
    "KAT6B",
    "KCNJ5",
    "KDM5A",
    "KDM5C",
    "KDM6A",
    "KDR",
    "KDSR",
    "KEAP1",
    "KIAA1549",
    "KIF5B",
    "KIT",
    "KLF4",
    "KLF6",
    "KMT2A",
    "KMT2C",
    "KMT2D",
    "KRAS",
    "KTN1",
    "LASP1",
    "LCK",
    "LCP1",
    "LEF1",
    "LHFP",
    "LIFR",
    "LMNA",
    "LMO1",
    "LMO2",
    "LPP",
    "LRIG3",
    "LYL1",
    "LZTR1",
    "MAF",
    "MAFB",
    "MALT1",
    "MAML2",
    "MAP2K1",
    "MAP2K2",
    "MAP2K4",
    "MAP3K1",
    "MAP3K13",
    "MAX",
    "MDM2",
    "MDM4",
    "MDS2",
    "MECOM",
    "MED12",
    "MEN1",
    "MET",
    "MITF",
    "MKL1",
    "MLF1",
    "MLH1",
    "MLLT1",
    "MLLT3",
    "MLLT4",
    "MLLT6",
    "MLLT10",
    "MLLT11",
    "MN1",
    "MNX1",
    "MPL",
    "MSH2",
    "MSH6",
    "MSI2",
    "MSN",
    "MTCP1",
    "MUC1",
    "MUTYH",
    "MYB",
    "MYC",
    "MYCL",
    "MYCN",
    "MYD88",
    "MYH9",
    "MYH11",
    "MYO5A",
    "NAB2",
    "NACA",
    "NBN",
    "NCOA1",
    "NCOA2",
    "NCOA4",
    "NCOR1",
    "NCOR2",
    "NDRG1",
    "NF1",
    "NF2",
    "NFATC2",
    "NFE2L2",
    "NFIB",
    "NFKB2",
    "NFKBIE",
    "NIN",
    "NKX2-1",
    "NONO",
    "NOTCH1",
    "NOTCH2",
    "NPM1",
    "NR4A3",
    "NRAS",
    "NSD1",
    "NSD2",
    "NT5C2",
    "NTRK1",
    "NTRK3",
    "NUMA1",
    "NUP98",
    "NUP214",
    "NUTM1",
    "NUTM2A",
    "OLIG2",
    "OMD",
    "P2RY8",
    "PAFAH1B2",
    "PALB2",
    "PATZ1",
    "PAX3",
    "PAX5",
    "PAX7",
    "PAX8",
    "PBRM1",
    "PBX1",
    "PCM1",
    "PDCD1LG2",
    "PDE4DIP",
    "PDGFB",
    "PDGFRA",
    "PDGFRB",
    "PER1",
    "PHF6",
    "PHOX2B",
    "PICALM",
    "PIK3CA",
    "PIK3CB",
    "PIK3R1",
    "PIM1",
    "PLAG1",
    "PLCG1",
    "PML",
    "PMS1",
    "PMS2",
    "POLE",
    "POT1",
    "POU2AF1",
    "POU5F1",
    "PPARG",
    "PPFIBP1",
    "PPP2R1A",
    "PPP6C",
    "PRCC",
    "PRDM1",
    "PRDM16",
    "PRF1",
    "PRKAR1A",
    "PRKCB",
    "PRRX1",
    "PSIP1",
    "PTCH1",
    "PTEN",
    "PTPN11",
    "PTPN13",
    "PTPRB",
    "PTPRC",
    "PTPRK",
    "RABEP1",
    "RAC1",
    "RAD21",
    "RAD51B",
    "RAF1",
    "RALGDS",
    "RANBP2",
    "RAP1GDS1",
    "RARA",
    "RB1",
    "RBM15",
    "RECQL4",
    "REL",
    "RET",
    "RHOA",
    "RHOH",
    "RMI2",
    "RNF43",
    "ROS1",
    "RPL5",
    "RPL10",
    "RPL22",
    "RPN1",
    "RSPO2",
    "RSPO3",
    "RUNX1",
    "RUNX1T1",
    "RUNX2",
    "S100A7",
    "SALL4",
    "SBDS",
    "SDC4",
    "SDHA",
    "SDHAF2",
    "SDHB",
    "SDHC",
    "SDHD",
    "SEPT5",
    "SEPT6",
    "SEPT9",
    "SET",
    "SETBP1",
    "SETD2",
    "SF3B1",
    "SFPQ",
    "SH2B3",
    "SH3GL1",
    "SLC34A2",
    "SLC45A3",
    "SMAD2",
    "SMAD3",
    "SMAD4",
    "SMARCA4",
    "SMARCB1",
    "SMARCD1",
    "SMARCE1",
    "SMO",
    "SND1",
    "SOCS1",
    "SOX2",
    "SOX21",
    "SPECC1",
    "SPEN",
    "SPOP",
    "SRC",
    "SRGAP3",
    "SRSF2",
    "SRSF3",
    "SS18",
    "SS18L1",
    "SSX1",
    "SSX2",
    "SSX4",
    "STAG2",
    "STAT3",
    "STAT5B",
    "STAT6",
    "STIL",
    "STK11",
    "STRN",
    "SUFU",
    "SUZ12",
    "SYK",
    "TAF15",
    "TAL1",
    "TAL2",
    "TBL1XR1",
    "TBX3",
    "TCEA1",
    "TCF3",
    "TCF7L2",
    "TCL1A",
    "TERT",
    "TET1",
    "TET2",
    "TFE3",
    "TFEB",
    "TFG",
    "TFPT",
    "TFRC",
    "THRAP3",
    "TLX1",
    "TLX3",
    "TMPRSS2",
    "TNC",
    "TNFAIP3",
    "TNFRSF14",
    "TNFRSF17",
    "TOP1",
    "TP53",
    "TP63",
    "TPM3",
    "TPM4",
    "TPR",
    "TRAF7",
    "TRIM24",
    "TRIM27",
    "TRIM33",
    "TRIP11",
    "TRRAP",
    "TSC1",
    "TSC2",
    "TSHR",
    "U2AF1",
    "UBR5",
    "USP6",
    "USP8",
    "VHL",
    "VTI1A",
    "WAS",
    "WHSC1",
    "WHSC1L1",
    "WIF1",
    "WRN",
    "WT1",
    "WWTR1",
    "XPA",
    "XPC",
    "XPO1",
    "YWHAE",
    "ZBTB16",
    "ZCCHC8",
    "ZEB1",
    "ZMYM2",
    "ZNF331",
    "ZNF384",
    "ZNF521",
    "ZRSR2",
    "EVI1",
    "MDS1",
    "MLL",
    "C15orf65",
}

_CIG = re.compile(r"(\d+)([MIDNSHP=X])")


def _cigar(cigar):
    p = _CIG.findall(cigar)
    if not p:
        return 0, 0, 0
    lead = int(p[0][0]) if p[0][1] == "S" else 0
    trail = int(p[-1][0]) if p[-1][1] == "S" else 0
    ref = sum(int(n) for n, op in p if op in "MDN=X")
    return lead, trail, ref


def parse_sam(path):
    with open(path) as fh:
        for line in fh:
            if not line.strip() or line.startswith("@"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 11:
                continue
            flag = int(f[1])
            if flag & 0x4 or f[2] == "*":
                continue
            pos, mapq, cigar = int(f[3]), int(f[4]), f[5]
            lead, trail, ref = _cigar(cigar)
            if lead >= trail:
                junction, clip = pos, lead
            else:
                junction, clip = pos + ref, trail
            yield {
                "chrom": f[2],
                "pos": pos,
                "mapq": mapq,
                "strand": "-" if flag & 0x10 else "+",
                "ref_span": ref,
                "clip_len": clip,
                "junction": junction,
                "molbar": f[0].split("_molbar")[0],
            }


def dedupe_umi(reads):
    best = {}
    for r in reads:
        k = r["molbar"]
        if k not in best or r["mapq"] > best[k]["mapq"]:
            best[k] = r
    return list(best.values())


def load_sample(raw):
    """Load configured SAM categories, using the source query-name-prefix rule."""
    out = {"HIGH": [], "TRUE_ALL": []}
    for key in out:
        path = sam_path(raw, key, required=(key == "HIGH"))
        if path is not None:
            out[key] = dedupe_umi(list(parse_sam(path)))
    return out


def call_sites(reads, window=JUNCTION_WINDOW):
    """Chain adjacent same-chromosome junctions within the selected window.

    Support counts distinct query-name prefixes. Return loci by decreasing
    support; centers are integer medians. Alignment strand is ignored.
    """
    by = defaultdict(list)
    for r in reads:
        by[r["chrom"]].append((r["junction"], r["molbar"]))
    sites = []
    for chrom, items in by.items():
        items.sort()
        cl = [items[0]]
        for j, mb in items[1:]:
            if j - cl[-1][0] <= window:
                cl.append((j, mb))
            else:
                sites.append(
                    {
                        "chrom": chrom,
                        "pos": int(np.median([c[0] for c in cl])),
                        "support": len({c[1] for c in cl}),
                    }
                )
                cl = [(j, mb)]
        sites.append(
            {
                "chrom": chrom,
                "pos": int(np.median([c[0] for c in cl])),
                "support": len({c[1] for c in cl}),
            }
        )
    return sorted(sites, key=lambda s: -s["support"])


class GeneModel:
    def __init__(self, path):
        path = str(path)
        self.by_chrom = defaultdict(list)
        op = gzip.open if path.endswith(".gz") else open
        with op(path, "rt") as fh:
            for line in fh:
                c = line.rstrip("\n").split("\t")
                if len(c) < 16:
                    continue
                name, chrom, strand = c[1], c[2], c[3]
                txS, txE, cdsS, cdsE = int(c[4]), int(c[5]), int(c[6]), int(c[7])
                exS = [int(x) for x in c[9].split(",") if x]
                exE = [int(x) for x in c[10].split(",") if x]
                gene = c[12]
                coding = name.startswith("NM_")
                self.by_chrom[chrom].append(
                    {
                        "gene": gene,
                        "strand": strand,
                        "txS": txS,
                        "txE": txE,
                        "cdsS": cdsS,
                        "cdsE": cdsE,
                        "exS": exS,
                        "exE": exE,
                        "coding": coding,
                    }
                )

        if not self.by_chrom:
            raise ValueError(
                "Gene reference contains no extended refGene/genePred rows; a nine-column GTF is not this input format"
            )
        self.starts = {}
        for chrom, lst in self.by_chrom.items():
            lst.sort(key=lambda t: t["txS"])
            self.starts[chrom] = [t["txS"] for t in lst]

    def annotate(self, chrom, pos):
        """Return the source feature, distance, strand, gene-list flag, and label.

        Intragenic distance is strand-aware TSS distance. Intergenic distance
        is signed genomic distance to the nearest transcript boundary.
        """
        lst = self.by_chrom.get(chrom, [])
        containing = [t for t in lst if t["txS"] <= pos < t["txE"]]
        if containing:
            best, best_rank = None, -1
            for t in containing:
                region = self._feature(t, pos)
                rank = {"CDS exon": 4, "UTR exon": 3, "exon": 3, "intron": 2}.get(
                    region, 1
                )
                rank += 0.5 if t["coding"] else 0
                if rank > best_rank:
                    best, best_rank, best_region = t, rank, region
            tss = best["txS"] if best["strand"] == "+" else best["txE"]
            dist = (pos - tss) if best["strand"] == "+" else (tss - pos)
            gene = best["gene"]
            return {
                "gene": gene,
                "region": best_region,
                "dist_tss": dist,
                "strand": best["strand"],
                "oncogene": gene.upper() in CANCER_GENES,
                "label": f"{gene} ({best_region})",
            }

        near, d = self._nearest(chrom, pos)
        if near is None:
            return {
                "gene": "intergenic",
                "region": "intergenic",
                "dist_tss": None,
                "strand": ".",
                "oncogene": False,
                "label": "intergenic",
            }
        gene = near["gene"]
        sign = "+" if d > 0 else "-"
        kb = abs(d) / 1000.0
        return {
            "gene": gene,
            "region": "intergenic",
            "dist_tss": d,
            "strand": near["strand"],
            "oncogene": gene.upper() in CANCER_GENES,
            "label": f"{gene} ({sign}{kb:.0f} kb)",
        }

    def transcript_for(self, chrom, pos):
        """Representative transcript spanning pos (prefer coding, then longest)."""
        cont = [t for t in self.by_chrom.get(chrom, []) if t["txS"] <= pos < t["txE"]]
        if not cont:
            return None
        cont.sort(key=lambda t: (t["coding"], t["txE"] - t["txS"]), reverse=True)
        return cont[0]

    def _feature(self, t, pos):
        in_exon = any(s <= pos < e for s, e in zip(t["exS"], t["exE"]))
        if in_exon:
            if t["cdsS"] <= pos < t["cdsE"] and t["cdsS"] < t["cdsE"]:
                return "CDS exon"
            return "UTR exon"
        return "intron"

    def _nearest(self, chrom, pos):
        lst = self.by_chrom.get(chrom, [])
        if not lst:
            return None, None
        bestd, bestt = None, None

        for t in lst:
            if t["txS"] <= pos <= t["txE"]:
                return t, 0
            d = (t["txS"] - pos) if pos < t["txS"] else (t["txE"] - pos)
            if bestd is None or abs(d) < abs(bestd):
                bestd, bestt = d, t
        return bestt, bestd


if __name__ == "__main__":
    SAMPLE_LABEL, _ = sample_layout()
    gm = GeneModel(REFGENE)
    print(
        f"Loaded gene model: {sum(len(v) for v in gm.by_chrom.values())} transcripts, "
        f"{len(gm.by_chrom)} contigs\n"
    )
    rows = []
    for raw, lab in SAMPLE_LABEL.items():
        d = load_sample(raw)
        sites = call_sites(d["HIGH"], window=100)
        total = sum(s["support"] for s in sites) or 1
        print(f"=== {lab}  ({len(sites)} loci, {total} molecules) — top 6 ===")
        for s in sites[:6]:
            a = gm.annotate(s["chrom"], s["pos"])
            flag = "  *ONCOGENE*" if a["oncogene"] else ""
            print(
                f"  {s['chrom']}:{s['pos']:<10d} {s['support']:>5d} mol "
                f"({100 * s['support'] / total:4.1f}%)  {a['label']}{flag}"
            )
        print()
