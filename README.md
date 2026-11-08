# Progressive Multiple Sequence Alignment from Scratch

A Python implementation of progressive multiple sequence alignment (MSA) using the Needleman–Wunsch algorithm, BLOSUM62 substitution scoring, UPGMA guide tree construction, and profile–profile alignment.

This project demonstrates how multiple protein sequences can be aligned progressively without relying on pre-built sequence alignment libraries.

## Features

- **Pairwise Global Alignment:** Implements Needleman–Wunsch using dynamic programming and traceback.
- **BLOSUM62 Scoring:** Transposes the provided substitution matrix and increases diagonal scores by 2.
- **Sequence Similarity:** Calculates the proportion of identical amino acids in each pairwise alignment.
- **UPGMA Guide Tree:** Builds a hierarchical guide tree based on average sequence similarity.
- **Progressive Alignment:** Aligns sequence groups using profile–profile dynamic programming.
- **Conserved Positions:** Identifies columns containing identical amino acids across all sequences.
- **Gap Penalty Comparison:** Supports different gap penalties to examine their effects on alignment results.

## How It Works

The implementation follows three main stages:

**1. Pairwise Sequence Alignment**

Every pair of protein sequences is globally aligned using Needleman–Wunsch. Alignment scores are calculated using the modified BLOSUM62 matrix and a linear gap penalty.

Sequence similarity is calculated as:

\[
\text{Similarity} = \frac{\text{Identical aligned positions}}{\text{Alignment length}}
\]

**2. UPGMA Guide Tree Construction**

A similarity matrix is constructed from pairwise alignments. At each iteration, the two most similar clusters are merged, and cluster similarities are recalculated using the average similarity between their members.

**3. Progressive Multiple Sequence Alignment**

The guide tree determines the order in which sequences are aligned. Existing alignments are treated as profiles and combined using dynamic programming. Gaps introduced into a group are preserved during subsequent merging steps.

The final alignment is displayed in blocks of 50 columns, with fully conserved positions marked by `*`.

## Project Structure

```text
progressive-multiple-sequence-alignment/
├── main.py
├── Input.txt
├── Blosum62.txt
└── README.md
```

## How to Run

**Requirements:** Python 3. No third-party Python packages are required.

Clone the repository:

```bash
git clone https://github.com/emirgencler/progressive-multiple-sequence-alignment.git
cd progressive-multiple-sequence-alignment
```

Run the program:

```bash
python main.py
```

When prompted, enter a negative gap penalty, such as `-8` or `-4`.

Ensure that `Input.txt` and `Blosum62.txt` are in the working directory.

## Example Results

The provided dataset contains five protein sequences: HUMAN, RABIT, FUKDA, PIG, and CHICK.

| Result | Gap = -8 | Gap = -4 |
|---|---|---|
| Most similar pair | HUMAN–RABIT | HUMAN–RABIT |
| Pairwise similarity | 0.92857 | 0.92857 |
| Final alignment length | 156 | 157 |
| Fully conserved columns | 60 | 60 |
| Second cluster addition | PIG | CHICK |

**Guide tree (gap = -8):**

```text
((((HUMAN,RABIT),PIG),CHICK),FUKDA)
```

**Guide tree (gap = -4):**

```text
((((HUMAN,RABIT),CHICK),PIG),FUKDA)
```

The results show that changing the gap penalty can affect both the guide tree topology and the final alignment length.

## Limitations

This is an educational implementation rather than a production-level alignment tool.

It uses linear gap penalties and progressive alignment without iterative refinement. As a result, gaps introduced in earlier stages may affect subsequent alignments.

## Technologies

- Python 3
- Dynamic Programming
- Needleman–Wunsch
- BLOSUM62
- UPGMA
- Progressive Multiple Sequence Alignment

## Author

Emir Gençler  
Computer Engineering  
Academic Project — Fall 2026
