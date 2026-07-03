## Overview

Minibwa combines bwa-mem seeding with minimap2 chaining and base alignment. It
aims to replace bwa-mem for short-read alignment with better performance at
comparable accuracy. It also works for accurate long read alignment to some
extent. The following table summarizes major components in minibwa and their
relationship with bwa-mem, ropebwt3 and minimap2:

|Component|Source  | Modifications |
|:--------|:-------|:--------------|
|FM-index |bwa-mem | Near identical |
|SMEM     |ropebwt3| Reimplemented and batched |
|SA query |ropebwt3| Batched differently |
|Seeding  |bwa-mem | Different algorithm |
|Chaining |minimap2| Adapted for variable-length seeds |
|Alignment|minimap2| Same ksw2 routines but different gluing code |
|Pairing  |bwa-mem | With pre-alignment filtering |

## Similarity to bwa-mem and minimap2

Similarity to bwa-mem:

 * Almost the same data structure for FM-index. Nonetheless, at the API level,
   minibwa uses half-closed-half-open suffix array (SA) intervals, similar to
   ropebwt3, but bwa-mem uses closed intervals. Minibwa also increases the SA
   sampling rate from 1/32 to 1/16.

Similarity to minimap2:

 * Minibwa chaining is adapted from minimap2, with modifications for
   variable-length seeds in minibwa. Bwa-mem tree-based chaining is no longer
   used.

 * Like minimap2, minibwa uses ksw2 for base alignment. Bwa-mem extension code
   is removed.

## Differences from bwa-mem and minimap2

 * Minibwa constructs FM-index with libsais by default. This library supports
   multi-threading (requiring OpenMP). It is much faster at the cost of more
   memory. The old bwa-mem BWT construction algorithm is still available.

 * Although the memory layout of FM-index is near identical, the index file
   formats are very different now. Minibwa introduces `l2bit` which is similar
   to UCSC's 2-bit format but supports contigs longer than 32-bit integers.

 * Minibwa uses Travis Gagie's algorithm for SMEM finding (function
   `mb_bwt_smem`). This algorithm was first implemented in ropebwt3. Although
   for standard SMEM finding in short reads, it is only slightly faster than
   the original algorithm, the new method is simpler and enables further
   optimization:

   - Minibwa caches the SA intervals of all 10-mers (`mb_bwt_cache`) which
     reduces random memory access. The old SMEM algorithm is incompatible with
     this pattern.

   - Minibwa batches multiple sequences for SMEM finding to hide latency with
     prefetch (`mb_bwt_smem_batch`), speeding up SMEM finding by 2.5 folds.
     Coding agents offerred great help in debugging the complex logic involved.
     It is possible to accelerate the old SMEM algorithm with latency hiding
     but this will be much more complex.

   - Minibwa uses two rounds of the SMEM algorithm for seeding
     (`mb_seed_intv`) and thus benefits from latency hiding throughout the
     seeding step (`mb_seed_intv_batch`). Bwa-mem does not use this strategy
     because with the old algorithm, the second round would be inefficient.

 * When retrieving SA values, minibwa also batches multiple BWT positions
   (`mb_bwt_sa_batch`). Bwa-mem2 uses the same trick.

 * Minibwa slightly speeds up SSE-based alignment originated from minimap2. The
   changes will be merged back to minimap2 later.

 * During mate rescue, minibwa uses a small 7-mer hash table to filter out poor
   alignment without full Smith-Waterman (`mb_ungap`). The algorithm is
   inspired by the Hough transform for line finding.

 * Minibwa reduces alignments in highly repetitive regions because short reads
   in these regions cannot be aligned well anyway. On simulated data, the
   accuracy remains similar.

 * Minibwa records soft-masked (lowercase) reference intervals in the `l2bit`
   file and honors them during mapping (`l2b_mask_overlap`). It drops
   high-occurrence anchors that fall entirely inside a soft-masked interval
   during seeding and caps the mapping quality at 10 when most of an alignment's
   reference span is soft-masked (`mb_cap_mapq_by_mask`). A reference without
   soft-masked bases is unaffected.

 * With `index --human`, minibwa sets the `L2B_F_NO_AMBI_SEED` flag in the
   `l2bit` file and records the intervals of ambiguous (`N`) reference bases. At
   mapping time it then discards any seed anchor overlapping such an interval
   (`l2b_has_ambi_global`), so no seeds fall inside runs of randomized `N`
   bases. This mode cannot be combined with the BWA pac output (`fa2bit -p`),
   which has no room to store the flag.

 * Minibwa natively aligns directional BS-seq reads.

 * With `index --human-resources`, minibwa imports HMF human resource BED/TSV
   files from an extracted directory or hmftools tarball (including nested
   `.bed.gz`) and writes an `<prefix>.hmf.tsv` sidecar with per-resource record
   counts, covered bases, FNV-1a checksums, and genome-build metadata
   (`mb_human_resources_import`). It rejects GRCh37/GRCh38 build mismatches at
   import time, and `map` warns when a loaded index disagrees with the sidecar
   build (`mb_human_resources_warn_mismatch`).
