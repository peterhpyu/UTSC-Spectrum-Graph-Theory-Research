import numpy as np

def is_strongly_connected(A):
    n = A.shape[0]
    M = np.linalg.matrix_power(A + np.eye(n, dtype=int), n)
    return np.all(M > 0)

def has_zero_row_or_col(M):
    """排除有全零行/列的情况（一定不强连通）"""
    return np.any(M.sum(1) == 0) or np.any(M.sum(0) == 0)

def build_matrix_from_bits(code, r, c):
    """把整数 code 的二进制展开成 r×c 的 0/1 矩阵"""
    bits = [(code >> k) & 1 for k in range(r * c)]
    return np.array(bits, dtype=int).reshape(r, c)

def two_cycle_count(B, C):
    """k = # {(i,j): B_ij=1 且 C_ji=1}"""
    return int(np.sum((B == 1) & (C.T == 1)))

def has_three_eigs_form(A, k, tol=1e-6):
    """判断 A 是否恰好有三种不同实特征值 ≈ {−√k, 0, √k}"""
    eigs = np.linalg.eigvals(A.astype(float))
    # 若存在显著虚部，排除
    if np.max(np.abs(eigs.imag)) > 1e-8:
        return False, eigs
    vals = np.sort(eigs.real)
    # 去重（合并相近的数）
    distinct = []
    for v in vals:
        if not distinct or abs(v - distinct[-1]) > tol:
            distinct.append(float(v))
    if len(distinct) != 3:
        return False, distinct
    t = np.sqrt(k)
    target = np.array([-t, 0.0, t])
    return np.allclose(np.sort(distinct), target, atol=tol, rtol=0.0), distinct

def canonical_key(B: np.ndarray, C: np.ndarray):
    """
    生成 (B,C) 的“转置等价类”规范键：
    在 (B,C) 与 (C^T, B^T) 两个表示中取字典序较小者。
    """
    B8 = B.astype(np.uint8, copy=False)
    C8 = C.astype(np.uint8, copy=False)
    k1 = (B8.shape, C8.shape, B8.tobytes(), C8.tobytes())
    k2 = (C8.T.shape, B8.T.shape, C8.T.tobytes(), B8.T.tobytes())
    return k1 if k1 <= k2 else k2

# ----- main -----
if __name__ == "__main__":
    output_file = "results.txt"  # output file name

    with open(output_file, "w", encoding="utf-8") as f:
        grand_total = 0     # num of all assembled A checked
        grand_good  = 0     # num of matrices satisfying Problem 1.1 (after dedup)
        dedup_skipped = 0   # how many were removed as transpose-duplicates

        seen_keys = set()   # <<< 新增：记录已输出的转置等价类代表

        for n in range(2, 7):  # n = 2,3,4,5,6
            good_this_n = 0     # count for current n
            f.write(f"\n=== n = {n} ===\n")
            print(f"\n=== n = {n} ===")

            for n1 in range(1, n//2 + 1):    # 只取一半划分，避免 (n1,n2)/(n2,n1) 跨划分重复
                n2 = n - n1
                kept_this_split = 0   # count for current (n1, n2)

                # enumerate B ∈ {0,1}^{n1×n2}
                for codeB in range(1 << (n1 * n2)):
                    B = build_matrix_from_bits(codeB, n1, n2)

                    # exclude B = J (trivial case)
                    if np.all(B == 1):
                        continue
                    # quick prune: not strongly connected anyway
                    if has_zero_row_or_col(B):
                        continue

                    # enumerate C ∈ {0,1}^{n2×n1}
                    for codeC in range(1 << (n2 * n1)):
                        C = build_matrix_from_bits(codeC, n2, n1)

                        # exclude C = J (trivial case)
                        if np.all(C == 1):
                            continue

                        if has_zero_row_or_col(C):
                            continue

                        # Assemble A from blocks B and C (not product)
                        A = np.block([
                            [np.zeros((n1, n1), dtype=int), B],
                            [C, np.zeros((n2, n2), dtype=int)]
                        ])

                        grand_total += 1

                        # strong connectivity test
                        if not is_strongly_connected(A):
                            continue

                        # compute k (2-cycles). If k=0, cannot have 3 distinct eigvals
                        k = two_cycle_count(B, C)
                        if k == 0:
                            continue

                        # check eigenvalue pattern {-sqrt(k), 0, sqrt(k)}
                        ok, eigs = has_three_eigs_form(A, k)
                        if not ok:
                            continue

                        # <<< 新增
                        key = canonical_key(B, C)
                        if key in seen_keys:
                            dedup_skipped += 1
                            continue
                        seen_keys.add(key)

                        kept_this_split += 1
                        good_this_n += 1
                        grand_good += 1

                        f.write(f"\nn={n}, n1={n1}, n2={n2}, k={k}, sqrt(k)={np.sqrt(k):.4f}\n")
                        f.write("B=\n" + str(B) + "\n")
                        f.write("C=\n" + str(C) + "\n")
                        f.write("A=\n" + str(A) + "\n")
                        f.write(f"Eigenvalues ≈ {eigs}\n")

                print(f"[n={n}, n1=?, n2=?] kept(after dedup)={kept_this_split}")
                f.write(f"[n={n}] split kept(after dedup)={kept_this_split}\n")

            print(f"==> Total kept for n={n}: {good_this_n}")
            f.write(f"==> Total kept for n={n}: {good_this_n}\n")

        # summary (still inside `with` to write safely)
        print("\n=== Summary ===")
        print(f"Scanned candidates (assembled A): {grand_total}")
        print(f"Good matrices total (after dedup): {grand_good}")
        print(f"Transpose-duplicates skipped: {dedup_skipped}")
        f.write("\n=== Summary ===\n")
        f.write(f"Scanned candidates (assembled A): {grand_total}\n")
        f.write(f"Good matrices total (after dedup): {grand_good}\n")
        f.write(f"Transpose-duplicates skipped: {dedup_skipped}\n")

    print(f"\n✅ All results have been saved to '{output_file}'")
