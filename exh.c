// Exhaustive search: A subset Z_m, 0 in A, |A| = k, R_A(n) <= U for all n (ordered count), and R_A(n) >= 1 for all n.
// Prints number of such sets and the first one found.
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
static int m, k, U;
static int A[128], r[256];
static unsigned long long cnt = 0, nodes = 0;
static int first[128], have = 0;
static int stop_first = 0;
static void rec(int depth, int start, int covered) {
    nodes++;
    if (depth == k) {
        if (covered == m) { if (!have) { memcpy(first, A, sizeof(int)*k); have = 1; } cnt++; }
        return;
    }
    if (stop_first && have) return;
    int rem = k - depth;
    // each new element x creates at most (depth) + 1 new covered residues, later ones more: bound on coverage
    // max newly covered by remaining rem elements: sum_{j=0}^{rem-1} (depth + j + 1)
    int maxnew = rem*depth + rem*(rem+1)/2;
    if (covered + maxnew < m) return;
    for (int x = start; x <= m - rem; x++) {
        int ok = 1, i, newc = 0;
        int d = (2*x) % m;
        for (i = 0; i < depth; i++) { int s = A[i] + x; if (s >= m) s -= m; if (r[s] == 0) newc++; r[s] += 2; if (r[s] > U) ok = 0; }
        if (r[d] == 0) newc++; r[d] += 1; if (r[d] > U) ok = 0;
        if (ok) { A[depth] = x; rec(depth+1, x+1, covered + newc); }
        r[d] -= 1;
        for (i = 0; i < depth; i++) { int s = A[i] + x; if (s >= m) s -= m; r[s] -= 2; }
        if (stop_first && have) return;
    }
}
int main(int argc, char **argv) {
    m = atoi(argv[1]); k = atoi(argv[2]); U = atoi(argv[3]); stop_first = argc > 4 ? atoi(argv[4]) : 0;
    memset(r, 0, sizeof r);
    A[0] = 0; r[0] = 1;
    rec(1, 1, 1);
    printf("m=%d k=%d U=%d count=%llu nodes=%llu", m, k, U, cnt, nodes);
    if (have) { printf(" first={"); for (int i = 0; i < k; i++) printf("%d%s", first[i], i+1<k?",":""); printf("}"); }
    printf("\n");
    return 0;
}
