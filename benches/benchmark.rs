use criterion::{black_box, criterion_group, criterion_main, Criterion};
use ark_std::str::FromStr;
use banderwagon::{Element, msm_gotti::MSMPrecompWnafGotti, Fr};
fn benchmark_precompute_mul(c: &mut Criterion) {
    // Create a vector of 256 elements, each being a multiple of the prime subgroup generator
    let basis_num = 1;
    let mut basic_crs = Vec::with_capacity(basis_num);
    for i in 0..basis_num {
        basic_crs.push(Element::prime_subgroup_generator() * Fr::from((i + 1) as u64));
    }

    let scalars = vec![
        Fr::from_str("13108968793781547619861935127046491459309155893440570251786403306729687672800").unwrap()
    ];

    let precompute = MSMPrecompWnafGotti::new(&basic_crs, 2, 8);

    c.bench_function("prempute_mul", |b| {
        b.iter(|| {
            let result = precompute.mul(black_box(&scalars));
            black_box(result);
        })
    });
}

criterion_group!(benches, benchmark_precompute_mul);
criterion_main!(benches);