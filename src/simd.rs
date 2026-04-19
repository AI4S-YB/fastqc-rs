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
#[inline]
pub fn count_n_per_position(seq: &[u8], n_counts: &mut [u64], not_n_counts: &mut [u64]) {
    debug_assert!(n_counts.len() >= seq.len());
    debug_assert!(not_n_counts.len() >= seq.len());

    for (i, &b) in seq.iter().enumerate() {
        if b == b'N' {
            n_counts[i] += 1;
        } else {
            not_n_counts[i] += 1;
        }
    }
}

/// Count A, C, G, T per position into pre-allocated vectors.
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

    for (i, &b) in seq.iter().enumerate() {
        match b {
            b'A' => a_counts[i] += 1,
            b'C' => c_counts[i] += 1,
            b'G' => g_counts[i] += 1,
            b'T' => t_counts[i] += 1,
            _ => {}
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
