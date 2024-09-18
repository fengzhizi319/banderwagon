use crate::salt_msm::WnafContext;
use ark_ed_on_bls12_381_bandersnatch::{EdwardsProjective, Fr};
use ark_ff::Zero;
use rayon::prelude::*;

use crate::Element;
#[derive(Clone, Debug)]
pub struct MSMPrecompWnaf {
    window_size: usize,
    pub tables: Vec<Vec<EdwardsProjective>>,
}

impl MSMPrecompWnaf {
    /// Creates a new MSMPrecompWnaf instance with precomputed tables.
    ///
    /// 创建一个带有预计算表的 MSMPrecompWnaf 实例。
    ///
    /// # Parameters
    /// - `bases`: A slice of `Element` representing the base points.
    /// - `window_size`: The window size for the WNAF context.
    ///
    /// # Returns
    /// A new `MSMPrecompWnaf` instance.
    ///
    /// # 参数
    /// - `bases`: 表示基点的 `Element` 切片。
    /// - `window_size`: WNAF 上下文的窗口大小。
    ///
    /// # 返回值
    /// 一个新的 `MSMPrecompWnaf` 实例。
    pub fn new(bases: &[Element], window_size: usize) -> MSMPrecompWnaf {
        let wnaf_context = WnafContext::new(window_size);

        // Parallel generation of precompute tables
        // 并行生成预计算表
        MSMPrecompWnaf {
            tables: bases.par_iter().map(|base|{
                wnaf_context.table(base.0)
            }).collect(),
            window_size,
        }
    }

    /// Multiplies a scalar with the precomputed table at the given index.
    ///
    /// 用给定索引处的预计算表乘以标量。
    ///
    /// # Parameters
    /// - `scalar`: The scalar to multiply.
    /// - `index`: The index of the precomputed table.
    ///
    /// # Returns
    /// An `Element` resulting from the multiplication.
    ///
    /// # 参数
    /// - `scalar`: 要乘的标量。
    /// - `index`: 预计算表的索引。
    ///
    /// # 返回值
    /// 乘法结果的 `Element`。
    pub fn mul_index(&self, scalar: Fr, index: usize) -> Element {
        let wnaf_context = WnafContext::new(self.window_size);
        Element(
            wnaf_context
                .mul_with_table(&self.tables[index], &scalar)
                .unwrap(),
        )
    }

    /// Performs multi-scalar multiplication (MSM) using the precomputed tables.
    ///
    /// 使用预计算表执行多标量乘法（MSM）。
    ///
    /// # Parameters
    /// - `scalars`: A slice of scalars to multiply.
    ///
    /// # Returns
    /// An `Element` resulting from the MSM.
    ///
    /// # 参数
    /// - `scalars`: 要乘的标量切片。
    ///
    /// # 返回值
    /// MSM 结果的 `Element`。
    pub fn mul(&self, scalars: &[Fr]) -> Element {
        let wnaf_context = WnafContext::new(self.window_size);
        let result: EdwardsProjective = scalars
            .iter()
            .zip(self.tables.iter())
            .filter(|(scalar, _)| !scalar.is_zero())
            .map(|(scalar, table)| wnaf_context.mul_with_table(table, scalar).unwrap())
            .sum();

        Element(result)
    }

    /// Performs parallel multi-scalar multiplication (MSM) using the precomputed tables.
    ///
    /// 使用预计算表执行并行多标量乘法（MSM）。
    ///
    /// # Parameters
    /// - `scalars`: A slice of scalars to multiply.
    ///
    /// # Returns
    /// An `Element` resulting from the parallel MSM.
    ///
    /// # 参数
    /// - `scalars`: 要乘的标量切片。
    ///
    /// # 返回值
    /// 并行 MSM 结果的 `Element`。
    pub fn mul_par(&self, scalars: &[Fr]) -> Element {
        // 创建一个新的 WNAF 上下文
        let wnaf_context = WnafContext::new(self.window_size);

        // 并行执行多标量乘法（MSM）
        let result: EdwardsProjective = scalars
            .par_iter() // 并行迭代标量
            .zip(self.tables.par_iter()) // 并行迭代预计算表
            .filter(|(scalar, _)| !scalar.is_zero()) // 过滤掉零标量
            .map(|(scalar, table)| wnaf_context.mul_with_table(table, scalar).unwrap()) // 使用 WNAF 表进行标量乘法
            .sum(); // 求和得到结果

        // 返回 MSM 结果的 `Element`
        Element(result)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{multi_scalar_mul, Element};

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
        let mut crs = Vec::with_capacity(256);
        for i in 0..256 {
            crs.push(Element::prime_subgroup_generator() * Fr::from((i + 1) as u64));
        }

        // Create a vector of 256 scalars, each being the negative of the corresponding index + 1
        // 创建一个包含 256 个标量的向量，每个标量都是相应索引 + 1 的负数
        let mut scalars = vec![];
        for i in 0..256 {
            scalars.push(-Fr::from(i + 1));
        }

        // Perform MSM using the standard method
        // 使用标准方法执行 MSM
        let result = multi_scalar_mul(&crs, &scalars);

        // Perform MSM using the precomputed WNAF method
        // 使用预计算的 WNAF 方法执行 MSM
        use std::time::Instant;

        let start = Instant::now();

        // Perform the operation
        let precomp = MSMPrecompWnaf::new(&crs, 12);

        let duration = start.elapsed();
        println!("Time elapsed in MSMPrecompWnaf_new() is: {:?}", duration);
        let got_result = precomp.mul(&scalars);
        let got_par_result = precomp.mul_par(&scalars);

        // Compare the results
        // 比较结果
        assert_eq!(result, got_result);
        assert_eq!(result, got_par_result);
    }
}
