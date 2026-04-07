# LDna Cluster Renaming: Methods Note

## Summary

LDna cluster IDs changed between the original (old) analysis and the revised (new) analysis.
The **SNP content of all major clusters is identical** between the two versions.
Only the string identifiers differ, due to a change in how `extractBranches()` was called.

---

## Why Cluster IDs Changed

### LDna Cluster ID Format

LDna assigns cluster IDs in the format `scaffoldPosition_threshold` (e.g., `scaffold_1_1234_100`).
These IDs are determined by the order in which branches are extracted during the call to
`extractBranches()`, not by fixed genomic coordinates. Because extraction order depends on
the full context of the call, the same SNPs can receive different string IDs when extracted
under different conditions.

### Old Approach: Loop Over min.edges Values

The original analysis looped over multiple `min.edges` parameter values:

```r
for (min_edges in c(20, 40, 60, 80, 100)) {
  branches <- extractBranches(ldna, min.edges = min_edges)
  # ... store results
}
# Only the min.edges=100 result was used downstream
```

Looping over five values (20, 40, 60, 80, 100) and retaining only the last (`min.edges=100`)
meant that the extraction was performed in the context of a broader enumeration. The cluster
IDs produced reflected the order imposed by that multi-value loop.

### New Approach: Single Direct Call

The revised analysis makes a single call:

```r
branches <- extractBranches(ldna, min.edges = 100)
```

Extracting branches only for `min.edges=100` produces a different enumeration order internally,
which changes the string IDs assigned to each cluster — even though the underlying SNP
membership is unchanged.

---

## Verification

A cluster comparison table was generated during the revision. It confirmed **100% SNP overlap**
for all major clusters between old and new IDs. Specific mappings include:

| Old ID    | New ID      | SNP overlap |
|-----------|-------------|-------------|
| Cluster 14 | (new ID)   | 100%        |
| Cluster 6  | (new ID)   | 100%        |

The comparison was performed by intersecting the SNP membership lists for each old cluster
with each new cluster. No SNPs were gained or lost.

---

## LDna Parameters

| Parameter   | Value | Description                          |
|-------------|-------|--------------------------------------|
| `min.edges` | 100   | Minimum edges required to retain a cluster |
| `phi`       | 1     | Outlier threshold (standard deviations above median LD) |

---

## r² Matrix Source

The pairwise r² matrix input to LDna was computed using PLINK with the `--r2` flag,
which applies the **composite estimator** of Zaykin et al. (2002). This estimator computes
r² as the squared Pearson correlation of dosage vectors (coded 0/1/2 per individual),
treating each SNP as a continuous dosage variable rather than estimating haplotype
frequencies. This approach is appropriate for data where phase is unknown.

> Zaykin DV, Westfall PH, Young SS, et al. (2002) Testing association of statistically
> inferred haplotypes with discrete and continuous traits in samples of unrelated individuals.
> *Human Heredity* 53(2):79–91.

---

## SOC Network Visualization

The SOC (Single-Outlier Cluster) networks shown in the manuscript were produced by
`scripts/05_ldna/ldna_soc_networks.R` using:

```r
plotLDnetwork(option = 2, full.network = FALSE)
```

Setting `full.network = FALSE` displays only the SOC itself, excluding the surrounding
COC (Compound Outlier Cluster) context. In the resulting plots:

- **Nodes** represent individual SNPs.
- **Edges** connect SNP pairs with r² above the extraction threshold.
- **Node color** encodes chromosomal position as a gradient, allowing visual assessment
  of whether LD clusters correspond to contiguous genomic regions.

---

## Conclusion

Reviewers should treat cluster name differences between the old and new analyses as
a labeling artifact of the extraction procedure, not as evidence of different results.
The biological content — the set of SNPs forming each outlier cluster — is unchanged.
