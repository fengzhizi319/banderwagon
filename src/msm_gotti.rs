use ark_ed_on_bls12_381_bandersnatch::{EdwardsProjective, Fr};
use ark_ff::Zero;
use rayon::prelude::*;
use ark_ec::CurveGroup;
//use ark_ff::{BigInteger, PrimeField};
use ark_std::vec::Vec;
use crate::Element;
use crate::salt_msm::WnafGottiContext;
#[derive(Clone, Debug)]
pub struct MSMPrecompWnafGotti {
    tables: Vec<Vec<EdwardsProjective>>,
    t:         usize,
    b:        usize,
}

impl MSMPrecompWnafGotti {
    pub fn fill_window<G: CurveGroup>(basis: &mut [G], table: &mut [G]) {
        // 如果 basis 数组为空
        if basis.is_empty() {
            // 将 table 数组填充为零元素
            for i in 0..table.len() {
                table[i] = G::zero();
            }
            return;
        }
        let l = table.len();
        // 递归调用 fill_window，传入去掉第一个元素的 basis 数组和 table 数组的前半部分
        Self::fill_window(&mut basis[1..], &mut table[.. l/ 2]);
        // 遍历 table 数组的前半部分
        for i in 0..l / 2 {
            // 将 table 数组后半部分的对应元素设置为当前元素与 basis 数组第一个元素的和
            table[l / 2 + i] = table[i] + basis[0];
        }
    }

    pub fn table<G: CurveGroup>(mut bases: Vec<G>, t: usize, b: usize) -> Vec<Vec<G>> {
        let fr_bits = 253;
        let window_size = 1 << b;
        let points_per_column = (fr_bits + t - 1) / t as usize;
        let window_count = (points_per_column + b - 1) / b;

        let mut final_nn_table = Vec::new();

        for base in bases.iter_mut() {
            let mut table_basis = Vec::with_capacity(points_per_column);
            let mut dbl = base.clone();
            table_basis.push(base.clone());

            for _i in 1..points_per_column {
                for _j in 0..t {
                    dbl = dbl.double();
                }
                table_basis.push(dbl.clone());
            }

            let mut nn_table = vec![G::zero(); window_count * window_size];
            for i in 0..window_count {
                let w = i;
                let start = w * b as usize;
                let mut end = (w + 1) * b as usize;
                if end > table_basis.len() {
                    end = table_basis.len();
                }
                let window_basis: &mut [G] = &mut table_basis[start..end];
                let mut table = vec![G::zero(); window_size];
                Self::fill_window(window_basis, &mut *table);

                for (j, table_element) in table.into_iter().enumerate() {
                    nn_table[w * window_size + j] = table_element;
                }
            }

            final_nn_table.push(nn_table);
        }

        final_nn_table
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
    pub fn mul(&self, scalars: &[Fr]) -> Element {
        let wnaf_gotti_context = WnafGottiContext::new(self.t,self.b);
        let result: EdwardsProjective = scalars
            .iter()
            .zip(self.tables.iter())
            .filter(|(scalar, _)| !scalar.is_zero())
            .map(|(scalar, table)| wnaf_gotti_context.mul_with_table_gotti(table, scalar).unwrap())
            .sum();

        Element(result)
    }
}



#[cfg(test)]
mod tests {
    use std::str::FromStr;
    use ark_ec::CurveGroup;
    //use ark_ff::PrimeField;
    use super::*;
    use crate::{multi_scalar_mul, Element};
    use crate::{msm_gotti::MSMPrecompWnafGotti, Fr};

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
        let mut basic_crs = Vec::with_capacity(256);
        for i in 0..256 {
            basic_crs.push(Element::prime_subgroup_generator() * Fr::from((i + 1) as u64));
        }
        // Create a vector of 256 scalars, each being the negative of the corresponding index + 1
        // 创建一个包含 256 个标量的向量，每个标量都是相应索引 + 1 的负数
        let mut scalars = vec![];
        for i in 0..256 {
            scalars.push(-Fr::from(i + 1));
        }
        // Perform MSM using the standard method
        // 使用标准方法执行 MSM
        let _result = multi_scalar_mul(&basic_crs, &scalars);
        // Perform MSM using the precomputed WNAF method
        // 使用预计算的 WNAF 方法执行 MSM
        //use std::time::Instant;
        //let start = Instant::now();

        // Perform the operation
        let _precomp = MSMPrecompWnafGotti::new(&basic_crs, 2,4);

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
        let mut basic_crs = Vec::with_capacity(1);

        basic_crs.push(Element::prime_subgroup_generator() * Fr::from(1));

        // Perform the operation
        let precompute=MSMPrecompWnafGotti::new(&basic_crs, 2,8);
        let x= Element::batch_extended_point_normalized(&*precompute.tables[0].clone());
        println!("precompute: {:?}", x.len());

    }
    #[test]
    fn correctness_precompute_one_table_test1() {
        // Create a vector of 256 elements, each being a multiple of the prime subgroup generator
        // 创建一个包含 256 个元素的向量，每个元素都是素数子群生成元的倍数
        let mut basic_crs = Vec::with_capacity(2);
        for i in 0..2 {
            basic_crs.push(Element::prime_subgroup_generator() * Fr::from((i + 1) as u64));
        }

        // Perform the operation
        let basic: Vec<_> = basic_crs.par_iter().map(|base| base.0).collect();
        let precompute_table = MSMPrecompWnafGotti::table(basic, 2, 8);
        //basic_crs.0
        let precompute=MSMPrecompWnafGotti::new(&basic_crs, 2,8);
        println!("precompute: {:?}", precompute_table.len());
        let x= Element::batch_extended_point_normalized(&*precompute.tables[0].clone());
        println!("precompute: {:?}", x.len());
        assert_eq!(precompute_table, precompute.tables);
        //let mut scalars = vec![];

    }
    #[test]
    fn correctness_msm_mul_one_test() {
        // Create a vector of 256 elements, each being a multiple of the prime subgroup generator
        // 创建一个包含 256 个元素的向量，每个元素都是素数子群生成元的倍数

        let basis_num = 1;
        let mut basic_crs = Vec::with_capacity(basis_num);
        for i in 0..basis_num {
            basic_crs.push(Element::prime_subgroup_generator() * Fr::from((i + 1) as u64));
        }
        let mut scalars = vec![];
        // for i in 0..basis_num {
        //     scalars.push(-Fr::from(i + 1));
        // }
        let scalars1=Fr::from_str("13108968793781547619861935127046491459309155893440570251786403306729687672800").unwrap();
        //q-1
        scalars.push(Fr::from_str("13108968793781547619861935127046491459309155893440570251786403306729687672800").unwrap());
        println!("scalars1: {:?}", scalars1.to_string());

        let precompute=MSMPrecompWnafGotti::new(&basic_crs, 2,8);
        use std::time::Instant;
        let start = Instant::now();
        let got_result = precompute.mul(&scalars);
        let duration = start.elapsed();
        println!("Time elapsed in mul is: {:?}", duration);

        let affine_result= got_result.0.into_affine();
        let string_x="33549696307925229982445904590536874618633472405590028303463218160177641247209";
        let string_y="19188667384257783945677642223292697773471335439753913231509108946878080696678";
        let x= affine_result.x.to_string();
        let y= affine_result.y.to_string();
        assert_eq!(string_x, x);
        assert_eq!(string_y, y);
        println!("got_result: {:?}", affine_result);
        //let mut scalars = vec![];

    }
}
