//! SIMD-accelerated primitives for sequence analysis using pulp.
//!
//! All functions have scalar fallbacks and produce bit-identical results
//! regardless of the available ISA (SSE2 / AVX2 / AVX-512 / scalar).

use pulp::Arch;

/// Count occurrences of A, C, G, T, N in `seq` and find the minimum byte value
/// in `qual`.
///
/// Returns `(a, c, g, t, n, quality_min)`.
#[inline]
pub fn count_bases_and_min_qual(seq: &[u8], qual: &[u8]) -> (u64, u64, u64, u64, u64, u8) {
    let arch = Arch::new();
    arch.dispatch(CountBasesAndMinQual { seq, qual })
}

struct CountBasesAndMinQual<'a> {
    seq: &'a [u8],
    qual: &'a [u8],
}

impl pulp::WithSimd for CountBasesAndMinQual<'_> {
    type Output = (u64, u64, u64, u64, u64, u8);

    #[inline(always)]
    fn with_simd<S: pulp::Simd>(self, simd: S) -> Self::Output {
        count_bases_and_min_qual_simd(simd, self.seq, self.qual)
    }
}

#[inline(always)]
fn count_bases_and_min_qual_simd<S: pulp::Simd>(
    simd: S,
    seq: &[u8],
    qual: &[u8],
) -> (u64, u64, u64, u64, u64, u8) {
    let mut a_count = 0u64;
    let mut c_count = 0u64;
    let mut g_count = 0u64;
    let mut t_count = 0u64;
    let mut n_count = 0u64;

    let (seq_head, seq_tail) = S::as_simd_u8s(seq);

    let splat_a = simd.splat_u8s(b'A');
    let splat_c = simd.splat_u8s(b'C');
    let splat_g = simd.splat_u8s(b'G');
    let splat_t = simd.splat_u8s(b'T');
    let splat_n = simd.splat_u8s(b'N');
    let zero = simd.splat_u8s(0);

    // Accumulators per lane (u8, so flush every 255 iterations to avoid overflow)
    let mut a_acc = simd.splat_u8s(0);
    let mut c_acc = simd.splat_u8s(0);
    let mut g_acc = simd.splat_u8s(0);
    let mut t_acc = simd.splat_u8s(0);
    let mut n_acc = simd.splat_u8s(0);

    let mut chunk_count = 0u32;

    for &chunk in seq_head {
        // equal_u8s returns m8s (mask type: 0xFF match, 0x00 no match)
        // Convert mask to u8s, then sub from zero to get 1 for match
        let a_mask = simd.transmute_u8s_m8s(simd.equal_u8s(chunk, splat_a));
        let c_mask = simd.transmute_u8s_m8s(simd.equal_u8s(chunk, splat_c));
        let g_mask = simd.transmute_u8s_m8s(simd.equal_u8s(chunk, splat_g));
        let t_mask = simd.transmute_u8s_m8s(simd.equal_u8s(chunk, splat_t));
        let n_mask = simd.transmute_u8s_m8s(simd.equal_u8s(chunk, splat_n));

        // mask is 0xFF for match; sub_u8s wraps: 0 - 0xFF = 1
        a_acc = simd.sub_u8s(a_acc, a_mask);
        c_acc = simd.sub_u8s(c_acc, c_mask);
        g_acc = simd.sub_u8s(g_acc, g_mask);
        t_acc = simd.sub_u8s(t_acc, t_mask);
        n_acc = simd.sub_u8s(n_acc, n_mask);

        chunk_count += 1;
        if chunk_count == 255 {
            a_count += horizontal_sum_u8s::<S>(a_acc);
            c_count += horizontal_sum_u8s::<S>(c_acc);
            g_count += horizontal_sum_u8s::<S>(g_acc);
            t_count += horizontal_sum_u8s::<S>(t_acc);
            n_count += horizontal_sum_u8s::<S>(n_acc);
            a_acc = zero;
            c_acc = zero;
            g_acc = zero;
            t_acc = zero;
            n_acc = zero;
            chunk_count = 0;
        }
    }

    // Flush remaining
    a_count += horizontal_sum_u8s::<S>(a_acc);
    c_count += horizontal_sum_u8s::<S>(c_acc);
    g_count += horizontal_sum_u8s::<S>(g_acc);
    t_count += horizontal_sum_u8s::<S>(t_acc);
    n_count += horizontal_sum_u8s::<S>(n_acc);

    // Scalar tail
    for &b in seq_tail {
        match b {
            b'A' => a_count += 1,
            b'C' => c_count += 1,
            b'G' => g_count += 1,
            b'T' => t_count += 1,
            b'N' => n_count += 1,
            _ => {}
        }
    }

    // Quality min — SIMD min reduction
    let mut min_q = 255u8;
    let (qual_head, qual_tail) = S::as_simd_u8s(qual);
    let mut min_acc = simd.splat_u8s(255);
    for &chunk in qual_head {
        min_acc = simd.min_u8s(min_acc, chunk);
    }
    // Reduce min across SIMD register
    let min_bytes = bytemuck::cast_slice::<S::u8s, u8>(core::slice::from_ref(&min_acc));
    for &b in min_bytes {
        if b < min_q {
            min_q = b;
        }
    }
    for &b in qual_tail {
        if b < min_q {
            min_q = b;
        }
    }

    (a_count, c_count, g_count, t_count, n_count, min_q)
}

/// Horizontal sum of a u8 SIMD register -> u64.
#[inline(always)]
fn horizontal_sum_u8s<S: pulp::Simd>(v: S::u8s) -> u64 {
    let bytes = bytemuck::cast_slice::<S::u8s, u8>(core::slice::from_ref(&v));
    bytes.iter().map(|&b| b as u64).sum()
}

/// Sum quality bytes and find the minimum. Returns `(sum, min)`.
#[inline]
pub fn sum_and_min_quality(qual: &[u8]) -> (u64, u8) {
    let arch = Arch::new();
    arch.dispatch(SumAndMinQuality { qual })
}

struct SumAndMinQuality<'a> {
    qual: &'a [u8],
}

impl pulp::WithSimd for SumAndMinQuality<'_> {
    type Output = (u64, u8);

    #[inline(always)]
    fn with_simd<S: pulp::Simd>(self, simd: S) -> Self::Output {
        sum_and_min_quality_simd(simd, self.qual)
    }
}

#[inline(always)]
fn sum_and_min_quality_simd<S: pulp::Simd>(simd: S, qual: &[u8]) -> (u64, u8) {
    let mut total_sum = 0u64;
    let mut min_val = 255u8;

    let (head, tail) = S::as_simd_u8s(qual);

    let mut min_acc = simd.splat_u8s(255);

    for &chunk in head {
        min_acc = simd.min_u8s(min_acc, chunk);
        // Sum: cast to bytes and accumulate
        let bytes = bytemuck::cast_slice::<S::u8s, u8>(core::slice::from_ref(&chunk));
        for &b in bytes {
            total_sum += b as u64;
        }
    }

    // Reduce min
    let min_bytes = bytemuck::cast_slice::<S::u8s, u8>(core::slice::from_ref(&min_acc));
    for &b in min_bytes {
        if b < min_val {
            min_val = b;
        }
    }

    // Scalar tail
    for &b in tail {
        total_sum += b as u64;
        if b < min_val {
            min_val = b;
        }
    }

    (total_sum, min_val)
}

/// Count N bases per position, adding to pre-allocated `n_counts` and `not_n_counts`.
///
/// Branchless: increment both arrays using arithmetic on the comparison result.
#[inline]
pub fn count_n_per_position(seq: &[u8], n_counts: &mut [u64], not_n_counts: &mut [u64]) {
    debug_assert!(n_counts.len() >= seq.len());
    debug_assert!(not_n_counts.len() >= seq.len());

    for (i, &b) in seq.iter().enumerate() {
        let is_n = (b == b'N') as u64;
        n_counts[i] += is_n;
        not_n_counts[i] += 1 - is_n;
    }
}

/// Count A, C, G, T per position into pre-allocated vectors.
///
/// Uses a 256-entry LUT to convert base → index (0..4) without branching.
/// Index 4 = "other" (skipped).  This eliminates costly branch mispredictions
/// in the inner loop.
#[inline]
pub fn count_bases_per_position(
    seq: &[u8],
    a_counts: &mut [u64],
    c_counts: &mut [u64],
    g_counts: &mut [u64],
    t_counts: &mut [u64],
) {
    debug_assert!(a_counts.len() >= seq.len());
    debug_assert!(c_counts.len() >= seq.len());
    debug_assert!(g_counts.len() >= seq.len());
    debug_assert!(t_counts.len() >= seq.len());

    // LUT: A=0, C=1, G=2, T=3, everything else=4
    const LUT: [u8; 256] = {
        let mut t = [4u8; 256];
        t[b'A' as usize] = 0;
        t[b'C' as usize] = 1;
        t[b'G' as usize] = 2;
        t[b'T' as usize] = 3;
        t
    };

    // Array of 4 mutable slice pointers for indexed access
    let counts: [&mut [u64]; 4] = [a_counts, c_counts, g_counts, t_counts];

    // Unrolled: process 4 bases at a time to amortize loop overhead
    let chunks = seq.len() / 4;
    let remainder = seq.len() % 4;

    for c in 0..chunks {
        let i = c * 4;
        let b0 = LUT[seq[i] as usize];
        let b1 = LUT[seq[i + 1] as usize];
        let b2 = LUT[seq[i + 2] as usize];
        let b3 = LUT[seq[i + 3] as usize];

        // Branchless: LUT value < 4 means valid base
        if b0 < 4 { counts[b0 as usize][i] += 1; }
        if b1 < 4 { counts[b1 as usize][i + 1] += 1; }
        if b2 < 4 { counts[b2 as usize][i + 2] += 1; }
        if b3 < 4 { counts[b3 as usize][i + 3] += 1; }
    }
    for i in (chunks * 4)..seq.len() {
        let idx = LUT[seq[i] as usize];
        if idx < 4 {
            counts[idx as usize][i] += 1;
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_count_bases_and_min_qual() {
        let seq = b"AACGGTTNNN";
        let qual = b"IIIIIIIBBB"; // I=73, B=66
        let (a, c, g, t, n, min_q) = count_bases_and_min_qual(seq, qual);
        assert_eq!(a, 2);
        assert_eq!(c, 1);
        assert_eq!(g, 2);
        assert_eq!(t, 2);
        assert_eq!(n, 3);
        assert_eq!(min_q, b'B');
    }

    #[test]
    fn test_sum_and_min_quality() {
        let qual = b"III"; // I = 73
        let (sum, min) = sum_and_min_quality(qual);
        assert_eq!(sum, 219);
        assert_eq!(min, 73);
    }

    #[test]
    fn test_count_bases_and_min_qual_empty() {
        let (a, c, g, t, n, min_q) = count_bases_and_min_qual(b"", b"");
        assert_eq!((a, c, g, t, n), (0, 0, 0, 0, 0));
        assert_eq!(min_q, 255);
    }

    #[test]
    fn test_count_bases_and_min_qual_long() {
        // Test with > 255 * SIMD_WIDTH bytes to exercise the flush logic
        let seq: Vec<u8> = (0..10000).map(|i| match i % 5 {
            0 => b'A',
            1 => b'C',
            2 => b'G',
            3 => b'T',
            _ => b'N',
        }).collect();
        let qual: Vec<u8> = (0..10000).map(|i| (33 + (i % 41) as u8)).collect();
        let (a, c, g, t, n, min_q) = count_bases_and_min_qual(&seq, &qual);
        assert_eq!(a, 2000);
        assert_eq!(c, 2000);
        assert_eq!(g, 2000);
        assert_eq!(t, 2000);
        assert_eq!(n, 2000);
        assert_eq!(min_q, 33);
    }

    #[test]
    fn test_count_n_per_position() {
        let seq = b"ANNG";
        let mut n_counts = vec![0u64; 4];
        let mut not_n_counts = vec![0u64; 4];
        count_n_per_position(seq, &mut n_counts, &mut not_n_counts);
        assert_eq!(n_counts, vec![0, 1, 1, 0]);
        assert_eq!(not_n_counts, vec![1, 0, 0, 1]);
    }

    #[test]
    fn test_count_bases_per_position() {
        let seq = b"ACGT";
        let mut a = vec![0u64; 4];
        let mut c = vec![0u64; 4];
        let mut g = vec![0u64; 4];
        let mut t = vec![0u64; 4];
        count_bases_per_position(seq, &mut a, &mut c, &mut g, &mut t);
        assert_eq!(a, vec![1, 0, 0, 0]);
        assert_eq!(c, vec![0, 1, 0, 0]);
        assert_eq!(g, vec![0, 0, 1, 0]);
        assert_eq!(t, vec![0, 0, 0, 1]);
    }

    #[test]
    fn test_sum_and_min_quality_long() {
        let qual: Vec<u8> = (0..5000).map(|i| (33 + (i % 41) as u8)).collect();
        let expected_sum: u64 = qual.iter().map(|&b| b as u64).sum();
        let expected_min: u8 = *qual.iter().min().unwrap();
        let (sum, min) = sum_and_min_quality(&qual);
        assert_eq!(sum, expected_sum);
        assert_eq!(min, expected_min);
    }
}

/// SIMD-accelerated per-position base accumulation using u8 accumulators.
///
/// Instead of branching per byte, uses 4× SIMD cmpeq to produce 0xFF masks,
/// then subtracts from u8 accumulators (0 - 0xFF wraps to +1).
/// Caller must flush accumulators to u64 counts every ≤255 calls.
#[inline]
pub fn accumulate_bases_simd(
    seq: &[u8],
    a_acc: &mut [u8],
    c_acc: &mut [u8],
    g_acc: &mut [u8],
    t_acc: &mut [u8],
) {
    debug_assert!(a_acc.len() >= seq.len());
    debug_assert!(c_acc.len() >= seq.len());
    debug_assert!(g_acc.len() >= seq.len());
    debug_assert!(t_acc.len() >= seq.len());

    struct AccOp<'a> {
        seq: &'a [u8],
        a_acc: &'a mut [u8],
        c_acc: &'a mut [u8],
        g_acc: &'a mut [u8],
        t_acc: &'a mut [u8],
    }

    impl pulp::WithSimd for AccOp<'_> {
        type Output = ();

        #[inline(always)]
        fn with_simd<S: pulp::Simd>(self, simd: S) -> Self::Output {
            let len = self.seq.len();

            let splat_a = simd.splat_u8s(b'A');
            let splat_c = simd.splat_u8s(b'C');
            let splat_g = simd.splat_u8s(b'G');
            let splat_t = simd.splat_u8s(b'T');

            let (seq_head, seq_tail) = S::as_simd_u8s(self.seq);
            let (a_head, a_tail) = S::as_mut_simd_u8s(&mut self.a_acc[..len]);
            let (c_head, c_tail) = S::as_mut_simd_u8s(&mut self.c_acc[..len]);
            let (g_head, g_tail) = S::as_mut_simd_u8s(&mut self.g_acc[..len]);
            let (t_head, t_tail) = S::as_mut_simd_u8s(&mut self.t_acc[..len]);

            // SIMD main loop: process full SIMD-width chunks
            for i in 0..seq_head.len() {
                let chunk = seq_head[i];

                // cmpeq → 0xFF where match, 0x00 elsewhere
                let ma = simd.transmute_u8s_m8s(simd.equal_u8s(chunk, splat_a));
                let mc = simd.transmute_u8s_m8s(simd.equal_u8s(chunk, splat_c));
                let mg = simd.transmute_u8s_m8s(simd.equal_u8s(chunk, splat_g));
                let mt = simd.transmute_u8s_m8s(simd.equal_u8s(chunk, splat_t));

                // sub: acc - 0xFF = acc + 1 (wrapping)
                a_head[i] = simd.sub_u8s(a_head[i], ma);
                c_head[i] = simd.sub_u8s(c_head[i], mc);
                g_head[i] = simd.sub_u8s(g_head[i], mg);
                t_head[i] = simd.sub_u8s(t_head[i], mt);
            }

            // Scalar tail
            for j in 0..seq_tail.len() {
                if seq_tail[j] == b'A' { a_tail[j] = a_tail[j].wrapping_add(1); }
                if seq_tail[j] == b'C' { c_tail[j] = c_tail[j].wrapping_add(1); }
                if seq_tail[j] == b'G' { g_tail[j] = g_tail[j].wrapping_add(1); }
                if seq_tail[j] == b'T' { t_tail[j] = t_tail[j].wrapping_add(1); }
            }
        }
    }

    Arch::new().dispatch(AccOp {
        seq,
        a_acc,
        c_acc,
        g_acc,
        t_acc,
    });
}

/// Flush u8 accumulators into u64 count arrays and reset accumulators to zero.
#[inline]
pub fn flush_base_accumulators(
    a_acc: &mut [u8],
    c_acc: &mut [u8],
    g_acc: &mut [u8],
    t_acc: &mut [u8],
    a_counts: &mut [u64],
    c_counts: &mut [u64],
    g_counts: &mut [u64],
    t_counts: &mut [u64],
    len: usize,
) {
    for i in 0..len {
        a_counts[i] += a_acc[i] as u64;
        c_counts[i] += c_acc[i] as u64;
        g_counts[i] += g_acc[i] as u64;
        t_counts[i] += t_acc[i] as u64;
        a_acc[i] = 0;
        c_acc[i] = 0;
        g_acc[i] = 0;
        t_acc[i] = 0;
    }
}

/// SIMD-accelerated per-position N/not-N accumulation using u8 accumulators.
#[inline]
pub fn accumulate_n_simd(
    seq: &[u8],
    n_acc: &mut [u8],
    not_n_acc: &mut [u8],
) {
    debug_assert!(n_acc.len() >= seq.len());
    debug_assert!(not_n_acc.len() >= seq.len());

    struct NOp<'a> {
        seq: &'a [u8],
        n_acc: &'a mut [u8],
        not_n_acc: &'a mut [u8],
    }

    impl pulp::WithSimd for NOp<'_> {
        type Output = ();

        #[inline(always)]
        fn with_simd<S: pulp::Simd>(self, simd: S) -> Self::Output {
            let len = self.seq.len();
            let splat_n = simd.splat_u8s(b'N');

            let (seq_head, seq_tail) = S::as_simd_u8s(self.seq);
            let (n_head, n_tail) = S::as_mut_simd_u8s(&mut self.n_acc[..len]);
            let (nn_head, nn_tail) = S::as_mut_simd_u8s(&mut self.not_n_acc[..len]);

            let all_ones = simd.splat_u8s(0xFF);
            for i in 0..seq_head.len() {
                let chunk = seq_head[i];
                // cmpeq → 0xFF where N, 0x00 elsewhere
                let is_n_u8 = simd.transmute_u8s_m8s(simd.equal_u8s(chunk, splat_n));
                // NOT via XOR with all 1s
                let not_n_u8 = simd.xor_u8s(is_n_u8, all_ones);

                // acc - 0xFF wraps to acc + 1
                n_head[i] = simd.sub_u8s(n_head[i], is_n_u8);
                nn_head[i] = simd.sub_u8s(nn_head[i], not_n_u8);
            }

            for j in 0..seq_tail.len() {
                if seq_tail[j] == b'N' {
                    n_tail[j] = n_tail[j].wrapping_add(1);
                } else {
                    nn_tail[j] = nn_tail[j].wrapping_add(1);
                }
            }
        }
    }

    Arch::new().dispatch(NOp { seq, n_acc, not_n_acc });
}

/// Flush N/not-N u8 accumulators into u64 counts.
#[inline]
pub fn flush_n_accumulators(
    n_acc: &mut [u8],
    not_n_acc: &mut [u8],
    n_counts: &mut [u64],
    not_n_counts: &mut [u64],
    len: usize,
) {
    for i in 0..len {
        n_counts[i] += n_acc[i] as u64;
        not_n_counts[i] += not_n_acc[i] as u64;
        n_acc[i] = 0;
        not_n_acc[i] = 0;
    }
}
