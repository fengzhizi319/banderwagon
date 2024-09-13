use ark_ed_on_bls12_381_bandersnatch::{EdwardsProjective, Fr};
use ark_ff::Zero;
use rayon::prelude::*;
use ark_ec::Group;
use ark_ff::{BigInteger, PrimeField};
use ark_std::vec::Vec;
use crate::Element;
use crate::salt_msm::WnafGottiContext;
#[derive(Clone, Debug)]
pub struct MSMPrecompWnafGotti {
    t:         usize,
    b:        usize,
    tables: Vec<Vec<EdwardsProjective>>,
}

impl MSMPrecompWnafGotti {
    pub fn table<G: Group>(mut base: G,t: usize, b: usize) -> Vec<G> {

        let fr_bits =253;
        //b为窗口bit长度，t为预处理的倍数，window_size为每个窗口内的值的所有可能的值
        let window_size = 1 << b;
        //用t把标量分为t个部分，如果t=2，那么标量分为2部分，每部分的位数为fr_bits/2，kP=2A1+A0，
        // A1为高位，A0为低位，A1，A0即为划分出来的部分。A1、A0的所有可能的个数为points_per_column，
        // 如t=2，b=5，那么points_per_column=127
        let points_per_column = (fr_bits + t - 1) / t as usize;
        //窗口的个数
        let window_count = (points_per_column + b - 1) / b; // 等价于 ceil(256 / self.window_size)
        let table_size = window_count * (1 << (b - 1));
        let mut table = Vec::with_capacity(table_size);
        let threshold = 1 << (b - 1);
        println!("threshold: {:?}", threshold);
        for _ in 0..window_count {
            // 初始化 current_base 为 base
            let mut current_base = base;
            // 遍历窗口中的元素,table中为G,2G,3G,...,2^(window_size-1)G
            table.push(current_base);
            for _ in 1..threshold {
                current_base += &base;
                current_base += &base; // Skip one element
                // 将 current_base 添加到表中
                table.push(current_base);
            }
            //2^(window_size-1)G+2^(window_size-1)G=2^(window_size)G
            base = current_base + current_base;
        }
        table
    }
    //b为窗口大小，t为预处理的倍数
    //pub  fn new(basis: &[Element], t: usize, b: usize) -> MSMPrecompWnafGotti {
    pub  fn new1(basis: &[Element], t: usize, window_size: usize) {
        //每个标量的位数
        let fr_bits =253;
        //b为窗口大小，t为预处理的倍数，window_size为每个窗口内的值的所有可能的值
        let window_size = 1 << window_size;
        //用t把标量分为t个部分，如果t=2，那么标量分为2部分，每部分的位数为fr_bits/2，kP=2A1+A0，
        // A1为高位，A0为低位，A1，A0即为划分出来的部分。A1、A0的所有可能的个数为points_per_column，
        // 如t=2，b=5，那么points_per_column=127
        let points_per_column = (fr_bits + t - 1) / t as usize;
        //窗口的个数
        let window_size = (points_per_column * basis.len() + window_size as usize - 1) / window_size as usize;
        let window_count = (256 + window_size - 1) / window_size; // 等价于 ceil(256 / self.window_size)
        let table_size = window_count * (1 << (window_size - 1));
        // let mut table = Vec::with_capacity(table_size);
        // let threshold = 1 << (window_size - 1);
        // println!("threshold: {:?}", threshold);

        // let mut table_basis = vec![bandersnatch::PointExtended::default(); points_per_column * basis.len()];
        //
        // let mut basis_extended = vec![bandersnatch::PointExtended::default(); basis.len()];
        // for (i, elem) in basis.iter().enumerate() {
        //     basis_extended[i] = bandersnatch::PointExtended::from_proj(&elem.inner);
        // }
        // let mut idx = 0;
        // for hi in 0..basis_extended.len() {
        //     table_basis[idx] = basis_extended[hi];
        //     idx += 1;
        //     for _ in 1..points_per_column {
        //         table_basis[idx] = table_basis[idx - 1];
        //         for _ in 0..t {
        //             table_basis[idx].double(&table_basis[idx]);
        //         }
        //         idx += 1;
        //     }
        // }
        // let mut nn_table = vec![bandersnatch::PointExtendedNormalized::default(); window_size * num_windows];
        // for w in 0..num_windows {
        //     let start = w * b as usize;
        //     let mut end = (w + 1) * b as usize;
        //     if end > table_basis.len() {
        //         end = table_basis.len();
        //     }
        //     let window_basis = &table_basis[start..end];
        //
        //     let mut table = vec![bandersnatch::PointExtended::default(); window_size];
        //     fill_window(window_basis, &mut table);
        //     let table_normalized = batch_to_extended_point_normalized(&table);
        //     for i in 0..window_size {
        //         nn_table[w * window_size + i] = table_normalized[i].clone();
        //     }
        // }

    }
    pub fn new(basis: &[Element], t: usize, b: usize)-> MSMPrecompWnafGotti {
        let wnaf_gotti_context = WnafGottiContext::new(t,b);

        // Parallel generation of precompute tables
        // 并行生成预计算表
        MSMPrecompWnafGotti {
            tables: basis.par_iter().map(|base|{
                wnaf_gotti_context.table(base.0)
            }).collect(),
            t,
            b
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{multi_scalar_mul, Element};
    use crate::{try_reduce_to_element, msm_gotti::MSMPrecompWnafGotti, Fr};

    #[test]
    fn testmain(){
        let s = "Hello world!".to_string();
        let s1 = s;
        // println!("s: {:?}", s); // 此行打开编译将报错
        println!("s1: {:?}", s1);
    }
    #[test]

    /// This test checks the correctness of the multi-scalar multiplication (MSM) implementation.
    ///
    /// 该测试检查多标量乘法（MSM）实现的正确性。
    ///
    /// # Explanation
    /// The test creates a vector of 256 elements, each being a multiple of the prime subgroup generator.
    /// It also creates a vector of 256 scalars, each being the negative of the corresponding index + 1.
    /// The test then performs MSM using both the standard and precomputed WNAF methods, and compares the results.
    ///
    /// # 说明
    /// 该测试创建一个包含 256 个元素的向量，每个元素都是素数子群生成元的倍数。
    /// 它还创建了一个包含 256 个标量的向量，每个标量都是相应索引 + 1 的负数。
    /// 然后，测试使用标准和预计算的 WNAF 方法执行 MSM，并比较结果。
    fn correctness_smoke_test() {
        // Create a vector of 256 elements, each being a multiple of the prime subgroup generator
        // 创建一个包含 256 个元素的向量，每个元素都是素数子群生成元的倍数
        let mut basic_g = Vec::with_capacity(256);
        for i in 0..256 {
            basic_g.push(Element::prime_subgroup_generator() * Fr::from((i + 1) as u64));
        }

        // Create a vector of 256 scalars, each being the negative of the corresponding index + 1
        // 创建一个包含 256 个标量的向量，每个标量都是相应索引 + 1 的负数
        let mut scalars = vec![];
        for i in 0..256 {
            scalars.push(-Fr::from(i + 1));
        }

        // Perform MSM using the standard method
        // 使用标准方法执行 MSM
        let result = multi_scalar_mul(&basic_g, &scalars);

        // Perform MSM using the precomputed WNAF method
        // 使用预计算的 WNAF 方法执行 MSM
        use std::time::Instant;

        let start = Instant::now();

        // Perform the operation
        let precomp = MSMPrecompWnafGotti::new(&basic_g, 2,4);

        // let duration = start.elapsed();
        // println!("Time elapsed in MSMPrecompWnaf_new() is: {:?}", duration);
        // let got_result = precomp.mul(&scalars);
        // let got_par_result = precomp.mul_par(&scalars);
        // 
        // // Compare the results
        // // 比较结果
        // assert_eq!(result, got_result);
        // assert_eq!(result, got_par_result);
    }
    #[test]
    fn correctness_precompute_one_table_test() {
        // Create a vector of 256 elements, each being a multiple of the prime subgroup generator
        // 创建一个包含 256 个元素的向量，每个元素都是素数子群生成元的倍数
        let mut basic_g = Vec::with_capacity(1);

        basic_g.push(Element::prime_subgroup_generator() * Fr::from(1));

        // Perform the operation
        let precompute = MSMPrecompWnafGotti::table(basic_g[0].0, 2,5);

    }
}
