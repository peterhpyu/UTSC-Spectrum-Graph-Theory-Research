> ⚠️ **Important Note on Confidentiality**
>
> This project is a academic research project, conducted in collaboration with Professor Michael Scott Cavers for a future co-authored paper.
>
> As this research contains **unpublished original findings** (e.g., the solutions to Problem 1.4 and new theoretical proofs), the core computational scripts and the `results.txt` dataset are kept in a **private repository** to prevent being "scooped" prior to publication on arXiv.
>
> This `README.md` document serves as a detailed overview of my personal contributions, methodologies, and the capabilities I developed during this research project.

---

# UTSC MATD92 Research: Spectral Graph Theory

## Overview
This project investigates open problems in spectral graph theory, specifically focusing on the classification and construction of bipartite directed graphs that possess exactly three distinct eigenvalues.

This research provides foundational mathematical insights for modern AI, particularly in **Large-scale Network Analysis, Graph Partitioning** and **Algebraic Graph Algorithms**.

* **Supervisor:** Professor Michael Cavers
* **Status:** Co-Author on in-progress academic paper ("Eigenvalues of digraphs")

Note on Public Code: ‘The spectral_graph_enumerator.py’ file in this public repository is a Reference example script from the Logic Verification of my research. It demonstrates the basic computational approach.
The main, unpublished computational scripts, which contain the advanced pruning algorithms and generated the final results (e.g., the solutions to Problem 1.4), remain in a private repository to protect the integrity of the future co-authored paper.

---

## My Contributions
As a co-author and the primary computational researcher, I completed the full research cycle: from theoretical analysis and experimental design to data generation and discovery.

### 1. ⚙️ Theoretical Analysis & Formalization
I contributed directly to the theoretical foundation of the paper, including:

* **Deriving the "Cloning Law":** I formulated the algebraic law ($k_{\text{new}} = k_{\text{core}} + k_{\delta}$) that precisely governs how the graph's $k$-value (number of 2-cycles) and spectrum evolve during a "cloning" (vertex duplication) operation.
* **Proposing a New Rank Theorem:** Based on computational evidence from my `n=7` (0 solutions) and `n=8` (found solutions) searches, I proposed and independently authored the formal proof for a new theoretical upper bound: $\text{rank}(BC) \le \min(n_1, n_2) - 1$. This argument was reviewed and confirmed by my supervisor.
* **Authoring Formal Propositions:** I authored the formal algebraic proof for Proposition 1.1 (“A General Construction for Rank-1 Product Matrices”), which provides a general recipe for constructing these specific graphs.

### 2. ⚡️ Computational Framework Development
I independently designed and developed the entire computational framework in Python/NumPy to exhaustively search for these rare graph structures.

* **High-Performance Pruning:** The search space for `n=8` exceeded 5.3 billion combinations. To make this computationally feasible on a local machine (requiring 50+ hours), I designed and implemented a multi-layer algebraic pruning system, including:
    * **Rank-based Pruning (Prune #1):** Skipping entire `(n1, n2)` splits (e.g., `n=7, n1=2, n2=5`) where $\min(n_1, n_2) < 3$, as they could never satisfy the $\text{rank} \ge 3$ condition.
    * **Factor-Rank Pruning (Prune #2 & #3):** Checking `rank(B) < 3` or `rank(C) < 3` before the main computation to skip millions of inner loops.
    * **Redundancy Pruning (Prune #1.5 & #2.5):** Using fast boolean checks (`has_identical_rows`) to quickly filter out matrices that would fail the more expensive `matrix_rank` calculations.
    * **Connectivity Pruning:** Using `has_zero_row_or_col` to pre-filter non-strongly connected graphs before the costly `matrix_power` operation.

### 3. 💡 Computational Discovery
My framework successfully led to a major breakthrough in the project:

* **Solved Problem 1.4:** My `n=8` script discovered the first-known computational examples of three-eigenvalue graphs where $\text{rank}(BC)=3$.
* **Significance:** This solved an open question posed by my supervisor and falsified the implicit hypothesis that the rank was limited to 2.
* **Unified the Theory:** This discovery, combined with my null-result at `n=7`, provided the **crucial evidence** needed to confirm my new theorem: $\text{rank}(BC) \le \min(n_1, n_2) - 1$.

---

## Author
Developed by Haopeng Yu  
📧 peter.hp.yu@gmail.com
🌐 GitHub: https://github.com/peterhpyu
