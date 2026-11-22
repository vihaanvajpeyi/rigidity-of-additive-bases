# Rigidity of additive bases with at most seven representations, and an improvement of Sandor's inequality -- Code and certificates


Code and data for the paper

> V. Vajpeyi, *Rigidity of additive bases with at most seven representations, and an improvement of Sándor's inequality*, .

The paper studies sets $A$ in a finite abelian group with $1 \le R_A(g) \le 7$ for every $g$, where $R_A(g)$ is the number of ordered pairs $(a,b) \in A^2$ with $a+b=g$. It also proves a strict improvement of Sándor's inequality for asymptotic bases of the integers. This repository contains everything needed to reproduce the computations in the paper: the search programs, the sets they found, and scripts that check the paper's identities, inequalities, and numerical constants.



## Checks

These run in a few seconds each and reproduce the main computational claims.

```bash
python3 certificates.py      # the six certificates of Appendix A
python3 verify.py            # identity (1.1) and its consequences on 36 sets
./exh_classify 37 4          # expect: count=40 size10:40
```

## Data

| File | Contents |
|---|---|
| `rm_U6.txt`, `found_103.txt`, `found_105.txt`, `s3a.txt`, `s3b.txt`, `found_s5.txt` | Sets $A \subseteq \mathbb{Z}_m$ with $1 \le R_A \le 6$ (or $\le 7$) found by the searches. This includes the certificates for $m = 101, 103, 104, 105, 107, 108$ (Appendix A). |
| `class_2_41.txt`, `class_42_48.txt` | Output of the exhaustive classification for small $m$ (Remark 2.5 and Section 8, item 1). |

## Scripts that check the paper

| Script | What it checks | Paper reference |
|---|---|---|
| `certificates.py` | Recomputes $R_A$ for the six certificates and the set in $\mathbb{Z}_{43}$ with $N_4 = 0$. Checks both sides of identity (1.1). | Appendix A; Section 8, items 2 and 4 |
| `verify.py` | Identity (1.1), Corollaries 1.3 and 1.4, and Proposition 4.1 on all 36 sets with $1 \le R_A \le 6$. | Theorem 1.2; Section 8, item 3 |
| `stats.py` | The same checks, plus the Singer one-exception table and the range of $N_4/m$. Imports `sidon.py`. | Table 2; Section 5; Section 8 |
| `sidon.py` | Singer, Bose, and Ruzsa Sidon sets in cyclic groups. The Singer set uses the convention stated in Section 5. | Proposition 5.1 |
| `randcheck.py` | Random checks of identity (1.1) and inequalities (1.2), (1.3), and (3.1) in cyclic and non-cyclic groups. | Section 8, item 3 |
| `constants.py` | Table 1, the constants $c_k$, the numbers in Corollaries 7.5 and 7.6, and the claim that the integer bound improves exactly when $K$ is a perfect square (checked for $K \le 20000$). | Sections 2 and 7 |
| `cov_identity_check.py` | The covariance identity (Lemma 7.1) and Corollary 7.2, to machine precision, for $F$ and for the cut-off function $\varphi F$. | Section 7 |
| `ident.py` | Random-set experiment: $V/M$ and the fourth moment of $G$. | Section 8, item 5 |
| `phi.py` | The same experiment for $\Phi = \varphi F$. Usage: `python3 phi.py SEED N`, e.g. `python3 phi.py 1 2000`. | Section 8, item 5 |
| `kappa_cutoff.py` | $\kappa$ for the cut-off function $\varphi F$ (close to 1). | Remark 7.7; Section 8, item 5 |
| `fourth_moment_uncut.py` | $\int\lvert G_F\rvert^4 / (\int\lvert G_F\rvert^2)^2$ for the uncut function (close to 3, for $250 \le N \le 16000$). | Remark 7.7 |
| `poisson_heuristic.py` | The crude Poisson count in "Why small moduli are not representative". | Section 8 |

## Search programs

| Program | Usage | Purpose |
|---|---|---|
| `exh_classify.c` | `./exh_classify m B [first]` | Exhaustive search over all $A \subseteq \mathbb{Z}_m$ with $0 \in A$ and $1 \le R_A \le B$. With `first`, stops at the first set found. |
| `exh.c` | `./exh m k U [stop_first]` | The original exhaustive search, for sets of size $k$ with $1 \le R_A \le U$. |
| `sa.c` | `./sa m k U seconds [seed]` | Simulated annealing for sets of size $k$ with $1 \le R_A \le U$. Used to find the certificates. |
| `sa_n4.c` | `./sa_n4 m k U seconds [seed]` | Annealing variant that searches for sets with $N_4 = 0$. |

To reproduce the classification in Remark 2.5 (for $m \le 47$, some $A \subseteq \mathbb{Z}_m$ has $1 \le R_A \le 5$ exactly when $m \le 25$ or $m \in \{27, 28, 35, 37, 39\}$):

```bash
for m in $(seq 2 47); do ./exh_classify $m 5 first; done
```

This takes a few minutes in total. Translating a set does not change whether it qualifies, so assuming $0 \in A$ loses nothing.

## Expected output

| Command | Expected |
|---|---|
| `python3 certificates.py` | All six sets have $\min R_A = 1$ and $\max R_A = 6$. For $m = 103$: $(N_1,\dots,N_6) = (9, 37, 5, 22, 5, 25)$, $t = 19$, $\mathcal{E}(A) = 48$, and both sides of (1.1) equal 136. |
| `python3 verify.py` | `all checks passed on 36 certificates` |
| `./exh_classify 37 4` | `count=40 size10:40` |
| `python3 constants.py` | $P_2 \le 0.02277$, $P_1 \ge 0.5552$, $P_3 \ge 0.2938$; $c_4 = 7 - 4\sqrt{3} = 0.07179\ldots$ |
| `python3 cov_identity_check.py` | Relative errors below $10^{-12}$ |

The annealing searches are randomized. They may need several runs or a longer time limit to find a set for a given $m$, and failing to find one says nothing about whether it exists.
