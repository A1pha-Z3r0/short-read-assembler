# README

**What’s a “contig” ? (and why we care)**

> “A contig is a set of DNA segments or sequences that overlap in a way that provides a contiguous representation of a genomic region.” 


In practice, de novo assembly turns many short reads into longer   contigs. Short reads are cheap and abundant, but overlap ambiguities (repeats, errors) make assembly tricky, especially when reads are very short:

> “…assembly of very short reads is still a challenging issue.” 

This repo implements a naïve but principled de Bruijn graph (DBG) approach that:

* Builds a weighted directed graph of k-mers,

* Finds Eulerian contigs when the graph structure supports it,

* Otherwise emits unitigs (maximal non-branching paths) and optionally chains them only when the next step is uniquely determined,

* Computes and plots per-base coverage from evidence(reads or k-mers), and **crucially breaks** when unsure to avoid mis-joins.


## Why this design?

* **Graph first, greediness last:** The graph summarizes evidence across all reads, so decisions (Eulerian/unitig) reflect global consistency rather than local, order-dependent merges.

* **Weighted edges encode multiplicity (coverage/repeats):** High-weight edges often indicate repeats; our logic naturally breaks at ambiguous branches.

* **Safety over length:** If the graph presents multiple options, we stop.


## Setup & Dependencies

**1. Create and activate virtual environment**
```bash
python -m venv venv
source venv/bin/activate
```

**2. Install dependencies**
```bash
pip install -r requirements.txt # on Windows: .venv\Scripts\activate
```

Libraries used: **biopython**, **networkx**, **matplotlib**.


## How to run app.py

**1. Find the direcotry of app.py and run this command** 
```bash
python app.py --input xyz.txt --output xyz.fasta --k 12
```

The two accepted formats for short reads inputs are:

**1. ```.txt``` files**

**2. ```.fasta``` files**

## Quick tips

* Choose k ≤ min(read length). If k is too large, you’ll get no k-mers.

* With very small k (e.g., 3–4), expect many branches → short unitigs.

* Duplicates are not errors—they reflect coverage and repeats.

* No reverse complements are used in this pipeline.


# References:
* Ji, Y., Shi, Y., Ding, G., & Li, Y. (2011). A new strategy for better genome assembly from very short reads. BMC Bioinformatics, 12, 493. https://doi.org/10.1186/1471-2105-12-493. 


* Zerbino, D. R., & Birney, E. (2008). Velvet: Algorithms for de novo short read assembly using de Bruijn graphs. Genome Research, 18(5), 821–829. https://doi.org/10.1101/gr.074492.107. 


* National Human Genome Research Institute. (2025, September 28). Contig. Talking Glossary of Genomic and Genetic Terms. https://www.genome.gov/genetics-glossary/Contig.

