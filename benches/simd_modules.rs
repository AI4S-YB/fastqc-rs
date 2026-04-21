use criterion::{black_box, criterion_group, criterion_main, Criterion};
use fastqc_rs::simd;

// ============ Scalar baselines ============

#[inline(never)]
fn scalar_count_bases_and_min_qual(seq: &[u8], qual: &[u8]) -> (u64, u64, u64, u64, u64, u8) {
    let mut a = 0u64;
    let mut c = 0u64;
    let mut g = 0u64;
    let mut t = 0u64;
    let mut n = 0u64;
    for &b in seq {
        match b {
            b'A' => a += 1,
            b'C' => c += 1,
            b'G' => g += 1,
            b'T' => t += 1,
            b'N' => n += 1,
            _ => {}
        }
    }
    let mut min_q = 255u8;
    for &b in qual {
        if b < min_q {
            min_q = b;
        }
    }
    (a, c, g, t, n, min_q)
}

#[inline(never)]
fn scalar_sum_and_min_quality(qual: &[u8]) -> (u64, u8) {
    let mut sum = 0u64;
    let mut min = 255u8;
    for &b in qual {
        sum += b as u64;
        if b < min {
            min = b;
        }
    }
    (sum, min)
}

// ============ Benchmarks ============

fn bench_count_bases_150bp(c: &mut Criterion) {
    let seq: Vec<u8> = (0..150).map(|i| match i % 5 {
        0 => b'A', 1 => b'C', 2 => b'G', 3 => b'T', _ => b'N',
    }).collect();
    let qual: Vec<u8> = (0..150).map(|i| 33 + (i % 41) as u8).collect();

    let mut group = c.benchmark_group("count_bases_150bp");
    group.bench_function("scalar", |b| {
        b.iter(|| scalar_count_bases_and_min_qual(black_box(&seq), black_box(&qual)))
    });
    group.bench_function("simd", |b| {
        b.iter(|| simd::count_bases_and_min_qual(black_box(&seq), black_box(&qual)))
    });
    group.finish();
}

fn bench_count_bases_10kb(c: &mut Criterion) {
    let seq: Vec<u8> = (0..10000).map(|i| match i % 5 {
        0 => b'A', 1 => b'C', 2 => b'G', 3 => b'T', _ => b'N',
    }).collect();
    let qual: Vec<u8> = (0..10000).map(|i| 33 + (i % 41) as u8).collect();

    let mut group = c.benchmark_group("count_bases_10kb");
    group.bench_function("scalar", |b| {
        b.iter(|| scalar_count_bases_and_min_qual(black_box(&seq), black_box(&qual)))
    });
    group.bench_function("simd", |b| {
        b.iter(|| simd::count_bases_and_min_qual(black_box(&seq), black_box(&qual)))
    });
    group.finish();
}

fn bench_sum_min_quality_150bp(c: &mut Criterion) {
    let qual: Vec<u8> = (0..150).map(|i| 33 + (i % 41) as u8).collect();

    let mut group = c.benchmark_group("sum_min_quality_150bp");
    group.bench_function("scalar", |b| {
        b.iter(|| scalar_sum_and_min_quality(black_box(&qual)))
    });
    group.bench_function("simd", |b| {
        b.iter(|| simd::sum_and_min_quality(black_box(&qual)))
    });
    group.finish();
}

fn bench_sum_min_quality_10kb(c: &mut Criterion) {
    let qual: Vec<u8> = (0..10000).map(|i| 33 + (i % 41) as u8).collect();

    let mut group = c.benchmark_group("sum_min_quality_10kb");
    group.bench_function("scalar", |b| {
        b.iter(|| scalar_sum_and_min_quality(black_box(&qual)))
    });
    group.bench_function("simd", |b| {
        b.iter(|| simd::sum_and_min_quality(black_box(&qual)))
    });
    group.finish();
}

criterion_group!(benches,
    bench_count_bases_150bp,
    bench_count_bases_10kb,
    bench_sum_min_quality_150bp,
    bench_sum_min_quality_10kb,
);
criterion_main!(benches);
