# Human-Focused Minibwa Alignment Context

This context defines the domain language for this fork of minibwa. The fork's goal is to improve minibwa specifically for human-genome mapping: faster and more accurate read placement against human reference assemblies, with architectural attention on repeat-heavy regions, reference representation, paired-end behavior, methylation support, and practical short-read workloads.

## Language

### Human Genome Focus

**Human Mapping Workload**:
The primary workload for this fork: mapping short reads and accurate long reads from human samples to a human Reference Assembly. Architecture suggestions should treat this workload as the optimization target unless a caller explicitly needs broad species-agnostic behavior.
_Avoid_: generic mapping workload, arbitrary genome workload

**Human Reference Assembly**:
The human genome assembly used as the Reference Genome, such as a primary GRCh38-style assembly or a complete T2T-style assembly. A Human Reference Assembly defines the Contigs, coordinate system, and repeat landscape that dominate mapping behavior.
_Avoid_: reference genome when the human-specific coordinate system matters

**Primary Assembly**:
The main set of assembled human chromosomes used for standard read placement. This is the baseline reference shape this fork should optimize before adding more complex graph, alternate, or pan-genome behavior.
_Avoid_: main genome, canonical genome

**Alternate Contig**:
A reference sequence representing an alternate haplotype or alternate locus rather than the Primary Assembly path. Current minibwa documentation says alternate contigs are not properly supported, so architecture proposals should treat them as an explicit limitation unless they directly address that gap.
_Avoid_: alt when prose needs the limitation to be clear

**Decoy Sequence**:
An extra reference sequence included to absorb reads from contaminants, unlocalized sequence, or hard-to-place genomic sequence. Decoy Sequences can improve practical human mapping but complicate hit ranking and Mapping Quality.
_Avoid_: extra contig

**Repeat-Heavy Region**:
A region of the Human Reference Assembly where many similar placements compete, including segmental duplications, centromeric sequence, simple repeats, and other high-copy sequence. Repeat-Heavy Regions are a first-class human-genome concern because they stress seeding, chaining, pairing, and Mapping Quality.
_Avoid_: repetitive region when the human mapping consequence matters

**Clinically Sensitive Region**:
A human locus where small mapping errors can matter disproportionately to downstream analysis, such as medically relevant genes, immune loci, or pharmacogenomic regions. This term names validation focus, not a separate alignment mode.
_Avoid_: important region, special locus

**Human Accuracy Regression**:
A change that worsens placement, Mapping Quality, pairing, or alignment detail on human truth sets or curated human stress regions, even if generic benchmarks still pass.
_Avoid_: accuracy regression when the human-specific benchmark scope matters

**Human Performance Regression**:
A change that worsens throughput, memory use, index size, or scaling on human-sized references and realistic human read sets.
_Avoid_: performance regression when the human reference size matters

### Reference and Indexing

**Reference Genome**:
The input FASTA sequence collection that Query Reads are mapped against. In this fork, the expected Reference Genome is usually a Human Reference Assembly. A Reference Genome contains one or more **Contigs**.
_Avoid_: genome file, target file, ref

**Contig**:
A named reference sequence with a length and offset in the indexed reference. A **Hit** points to exactly one Contig after coordinates are decoded from the indexed sequence space.
_Avoid_: chromosome, target, sequence name

**L2Bit Reference**:
Minibwa's 2-bit encoded Reference Genome, stored as `.l2b`, including contig metadata, packed bases, ambiguous intervals, and masked intervals. It is the coordinate authority used after the FM-index finds candidate positions.
_Avoid_: pac, 2bit file, reference blob

**FM-index**:
The searchable index over the Reference Genome, stored as `.mbw`, containing the BWT, cumulative base counts, an optional k-mer cache, and sampled suffix array values. One loaded **Index** owns one FM-index.
_Avoid_: BWT file when referring to the whole searchable structure

**BWT**:
The Burrows-Wheeler Transform sequence inside the FM-index. Use BWT only for the transform and rank/extend operations, not for the whole Index.
_Avoid_: index, FM-index

**Suffix Array Interval**:
A half-open interval in FM-index space representing all reference occurrences of a query substring. Seed finding produces Suffix Array Intervals before they are resolved into Anchors.
_Avoid_: SA range, occurrence range, interval when the suffix-array meaning matters

**Index**:
The loaded mapping structure that combines an **L2Bit Reference** and one **FM-index**. In the Human Mapping Workload, the Index must be understood as a large, repeat-dense, human-sized structure where memory layout and cache behavior are architectural concerns.
_Avoid_: database, reference handle, BWT

**Index Prefix**:
The filesystem stem used to locate the `.l2b`, `.mbw`, and optional methylation `.meth.mbw` files for an Index.
_Avoid_: reference path, output prefix

### Reads and Mapping

**Query Read**:
An input FASTA/FASTQ sequence to align against the Reference Genome. In the Human Mapping Workload, Query Reads are usually short-read sequencing reads, accurate long reads, or directional bisulfite reads. A Query Read may be mapped alone or as one segment of a Paired-End Fragment.
_Avoid_: query sequence, read sequence, sequence

**Read Segment**:
One Query Read within a multi-segment input record. In ordinary paired-end mapping, a Paired-End Fragment has exactly two Read Segments.
_Avoid_: mate when the relationship is not specifically paired-end

**Paired-End Fragment**:
Two Read Segments expected to originate from one DNA fragment with a supported orientation and insert size distribution.
_Avoid_: read pair when discussing the biological fragment

**Mini-batch**:
A chunk of Query Reads read from input and mapped together by worker threads. A Mini-batch may be split into **SMEM Sub-batches** for latency-hiding seed discovery.
_Avoid_: batch when the size and worker scheduling matter

**SMEM Sub-batch**:
A smaller group of Query Reads sent through batched SMEM discovery with shared prefetch and memory reuse.
_Avoid_: batch if it could mean the larger Mini-batch

### Seeds, Anchors, and Chains

**SMEM**:
A super-maximal exact match between a Query Read and the Reference Genome. SMEM discovery returns Suffix Array Intervals that may become Seed Intervals.
_Avoid_: seed when the exact-match maximality matters

**Seed Interval**:
A candidate exact-match interval found during seeding, represented in suffix-array space before reference coordinates are materialized.
_Avoid_: seed hit, raw seed

**Anchor**:
A materialized Seed Interval with query position, reference coordinate, strand, and length. Anchors are the input to chaining.
_Avoid_: match, seed when coordinates are already materialized

**Chain**:
An ordered group of Anchors that supports one candidate placement of a Query Read on a Contig. Chaining scores and filters Anchors before base alignment.
_Avoid_: region, cluster

**High-Copy Chain**:
A Chain dominated by Anchors from Repeat-Heavy Regions where many similar placements compete. High-Copy Chains are important because human-genome optimizations often trade speed against sensitivity in these regions.
_Avoid_: repetitive chain, bad chain

**Hit**:
A candidate or finalized mapping of a Query Read to a Contig, including query coordinates, reference coordinates, chaining score, mapping quality, flags, and optional alignment details. A Query Read can produce zero or more Hits.
_Avoid_: alignment when base-level alignment may not have run

**Primary Hit**:
The Hit selected as the main placement for output and mapping quality comparison. A Query Read has at most one Primary Hit per output record.
_Avoid_: best alignment

**Supplementary Hit**:
A Hit representing another segment of the same Query Read, usually from split mapping. It is not an alternative placement of the same aligned bases.
_Avoid_: secondary hit

**Secondary Hit**:
A non-primary alternative placement for the same Query Read. It competes with the Primary Hit for mapping quality.
_Avoid_: supplementary hit

### Alignment and Output

**Base Alignment**:
The Smith-Waterman-style refinement that turns a chained Hit into base-level coordinates and edit operations. Base Alignment may be skipped when only candidate mappings are needed.
_Avoid_: DP when discussing the user-visible mapping stage

**CIGAR**:
The run-length encoded edit path for a Base Alignment. A Hit has CIGAR data only after Base Alignment has populated its extra alignment details.
_Avoid_: cigar string when the stored form is packed operations

**Mapping Quality**:
The confidence score assigned to a Hit after competing Hits, chain scores, pairing evidence, and read mode are considered. In this fork, Mapping Quality should be evaluated against human-specific ambiguity, especially Alternate Contigs, Decoy Sequences, and Repeat-Heavy Regions.
_Avoid_: mapq only in prose

**SAM Output**:
The sequence alignment output format used for short reads, paired-end reads, and workflows that need CIGAR, flags, tags, and headers.
_Avoid_: alignment output when the format matters

**PAF Output**:
The pairwise mapping output format used for long-read or lightweight mapping workflows.
_Avoid_: tab output

**Alignment Tag**:
An optional SAM or PAF annotation derived from the Hit and Base Alignment, such as `cs`, `ds`, `MD`, or `XA`.
_Avoid_: extra field

### Paired-End and Methylation

**Insert Size Distribution**:
The observed distance and orientation model for Paired-End Fragments. Pairing and Mate Rescue use this distribution to judge proper pairs, and it is especially valuable for resolving ambiguous human short-read placements.
_Avoid_: PE stats, insert stats

**Mate Rescue**:
The paired-end recovery step that searches near one confident Read Segment to align its mate when ordinary mapping is weak or missing.
_Avoid_: rescue alignment unless referring only to the Base Alignment inside rescue

**Methylation Mode**:
The directional bisulfite sequencing mode where reads and the Index are interpreted through C-to-T or G-to-A converted sequence space.
_Avoid_: meth mode in prose

**Methylation Type**:
The conversion identity for a Query Read or reference coordinate: none, C-to-T, or G-to-A. Methylation Type affects seeding, coordinate decoding, and sequence conversion.
_Avoid_: strand when the conversion state is meant

**Adaptive Mapping Mode**:
The mode that chooses short-read or long-read internal parameters from Query Read length. It is a parameter-selection concept, not a separate mapper.
_Avoid_: auto mode, adaptive mapper

## Flagged Ambiguities

**Generic Minibwa vs Human-Focused Fork**:
Use **Human Mapping Workload** when discussing this fork's optimization target. Preserve generic minibwa behavior where it is cheap and already supported, but do not let hypothetical non-human workloads dominate architecture decisions.

**Index vs FM-index vs BWT**:
Use **Index** for the loaded `mb_idx_t` concept that owns both reference and search structures. Use **FM-index** for the searchable `.mbw` structure. Use **BWT** only for the transform and its rank/extend operations.

**Reference Genome vs Human Reference Assembly**:
Use **Human Reference Assembly** when the coordinate system, contig set, repeat landscape, or benchmark target matters. Use **Reference Genome** for the generic input accepted by the CLI and public interface.

**Primary Assembly vs Alternate Contig vs Decoy Sequence**:
Use **Primary Assembly** for the main chromosome paths. Use **Alternate Contig** for alternate haplotype or locus representations. Use **Decoy Sequence** for extra sequence intended to absorb otherwise misleading reads.

**Hit vs Base Alignment**:
Use **Hit** for a candidate mapping whether or not Smith-Waterman refinement has run. Use **Base Alignment** only for the refinement step or its resulting edit path.

**Anchor vs Seed Interval**:
Use **Seed Interval** while the candidate is still in suffix-array space. Use **Anchor** once it has query and reference coordinates.

**Secondary Hit vs Supplementary Hit**:
Use **Secondary Hit** for competing placements of the same Query Read. Use **Supplementary Hit** for split segments that together explain different parts of the same Query Read.

**Batch vs Sub-batch**:
Use **Mini-batch** for input and worker scheduling. Use **SMEM Sub-batch** for the smaller seed-discovery group controlled by SMEM batching options.

## Relationships

- One **Human Mapping Workload** uses one **Human Reference Assembly** as its usual Reference Genome.
- One **Human Reference Assembly** contains one **Primary Assembly** and may contain zero or more **Alternate Contigs** or **Decoy Sequences**.
- One **Reference Genome** contains one or more **Contigs**.
- One **Index Prefix** identifies one **L2Bit Reference**, one primary **FM-index**, and optionally one methylation **FM-index**.
- One loaded **Index** owns one **L2Bit Reference** and one active **FM-index**.
- One **Query Read** produces zero or more **Seed Intervals**.
- One **Seed Interval** may materialize into zero or more **Anchors**, depending on occurrence limits and coordinate decoding.
- One **Chain** contains one or more **Anchors**.
- One **Query Read** produces zero or more **Hits**.
- One **Hit** points to exactly one **Contig** after coordinate decoding.
- One **Hit** may have zero or one **CIGAR**.
- One **Paired-End Fragment** contains exactly two **Read Segments**.
- One **Insert Size Distribution** is estimated from many Paired-End Fragments.
- One **Repeat-Heavy Region** can produce many competing **Seed Intervals**, **Anchors**, **Chains**, and **Hits**.
- One **Human Accuracy Regression** or **Human Performance Regression** should be judged against human-sized references and human read sets.

## Example Dialogue

Developer: "This Query Read produced many Suffix Array Intervals in a Repeat-Heavy Region. Should this module return Anchors directly?"

Domain expert: "Only after coordinate decoding. While the candidates are still in FM-index space, call them Seed Intervals. Once they have query position, Contig coordinate, strand, and length, they are Anchors."

Developer: "The Chain module can then accept Anchors and return Hits, but those Hits may not have CIGAR data yet. For this fork, should we keep extra candidates in repeat-heavy human regions?"

Domain expert: "A Hit is a candidate placement. Base Alignment adds CIGAR and alignment tags later, unless the user requested mapping without Base Alignment. The human-specific question is whether keeping those candidates improves Mapping Quality or just creates a Human Performance Regression."

Developer: "For paired-end input, should pairing work on Query Reads or Paired-End Fragments?"

Domain expert: "Pairing reasons about Paired-End Fragments. Each fragment has two Read Segments, and the Insert Size Distribution tells us whether their Hits form a proper pair or need Mate Rescue."

Developer: "Should an architecture review assume Alternate Contigs are supported?"

Domain expert: "No. Current minibwa says Alternate Contigs are not properly supported. Treat that as an explicit limitation, or propose a module that deepens the reference-coordinate model enough to address it."
