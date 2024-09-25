use ark_ed_on_bls12_381_bandersnatch::{EdwardsProjective, Fr};
use ark_ff::Zero;
use rayon::prelude::*;
use ark_ec::PrimeGroup;
//use ark_ff::{BigInteger, PrimeField};
use ark_std::vec::Vec;
use crate::Element;
use crate::scalar_multi::WnafGottiContext;
#[derive(Clone, Debug)]
pub struct MSMPrecompWnafGotti {
    pub tables: Vec<Vec<EdwardsProjective>>,
    t:         usize,
    b:        usize,
}

impl MSMPrecompWnafGotti {
    pub fn fill_window<G: PrimeGroup>(basis: &mut [G], table: &mut [G]) {
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
            .map(|(scalar, table)| wnaf_gotti_context.mul_with_table(table, scalar).unwrap())
            .sum();

        Element(result)
    }

    pub fn gotti_window(&self, scalars: &[Fr]) -> Vec<Vec<u16>> {
        let wnaf_gotti_context = WnafGottiContext::new(self.t,self.b);
        let scalar=scalars[0];
        let result: Vec<Vec<u16>> = wnaf_gotti_context.gotti_window::<EdwardsProjective>(&scalar);
        result
    }

}



#[cfg(test)]
mod tests {
    use std::str::FromStr;
    use ark_ec::{CurveGroup};
    use ark_std::rand::SeedableRng;
    use ark_std::UniformRand;
    use rand_chacha::ChaCha20Rng;
    //use ark_ff::PrimeField;
    use super::*;
    use crate::{ Element};
    use crate::{msm_gotti::MSMPrecompWnafGotti, Fr};
    use crate::msm_window::MSMPrecompWnaf;

    #[test]
    fn testmain(){
        let s = "Hello world!".to_string();
        let s1 = s;
        // println!("s: {:?}", s); // 此行打开编译将报错
        println!("s1: {:?}", s1);
    }
    #[test]
    fn correctness_gotti_precompute_one_table_test() {
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
    fn correctness_gotti_msm_mul_one_test() {
        // Create a vector of 256 elements, each being a multiple of the prime subgroup generator
        // 创建一个包含 256 个元素的向量，每个元素都是素数子群生成元的倍数

        let basis_num = 1;
        let mut basic_crs = Vec::with_capacity(basis_num);
        for i in 0..basis_num {
            basic_crs.push(Element::prime_subgroup_generator() * Fr::from((i + 1) as u64));
        }
        let mut scalars = vec![];
        //q-1
        scalars.push(Fr::from_str("13108968793781547619861935127046491459309155893440570251786403306729687672800").unwrap());


        let precompute=MSMPrecompWnafGotti::new(&basic_crs, 2,4);
        let mem_byte_size=precompute.tables.len()*precompute.tables[0].len()*4*32;
        println!("precompute_size: {:?}", mem_byte_size);
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


    #[test]
    fn correctness_gotti_window() {
        // Create a vector of 256 elements, each being a multiple of the prime subgroup generator
        // 创建一个包含 256 个元素的向量，每个元素都是素数子群生成元的倍数

        let basis_num = 1;
        let mut basic_crs = Vec::with_capacity(basis_num);
        for i in 0..basis_num {
            basic_crs.push(Element::prime_subgroup_generator() * Fr::from((i + 1) as u64));
        }
        let mut scalars = vec![];
        //q-1
        scalars.push(Fr::from_str("13108968793781547619861935127046491459309155893440570251786403306729687672800").unwrap());

        let precompute=MSMPrecompWnafGotti::new(&basic_crs, 2,4);

        let _got_result = precompute.gotti_window(&scalars);

        //let mut scalars = vec![];

    }
    #[test]
    fn correctness_gotti_msm_mul_one_test_normal() {
        // Create a vector of 256 elements, each being a multiple of the prime subgroup generator
        // 创建一个包含 256 个元素的向量，每个元素都是素数子群生成元的倍数

        let basis_num = 1;
        let mut basic_crs = Vec::with_capacity(basis_num);
        for i in 0..basis_num {
            basic_crs.push(Element::prime_subgroup_generator() * Fr::from((i + 1) as u64));
        }
        let mut scalars = vec![];
        //q-1
        scalars.push(Fr::from_str("13108968793781547619861935127046491459309155893440570251786403306729687672800").unwrap());


        let precompute=MSMPrecompWnafGotti::new(&basic_crs, 2,4);
        let mem_byte_size=precompute.tables.len()*precompute.tables[0].len()*4*32;
        println!("precompute_size: {:?}", mem_byte_size);
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
    #[test]
    fn correctness_window_msm_one_g() {
        // Create a vector of 256 elements, each being a multiple of the prime subgroup generator
        // 创建一个包含 256 个元素的向量，每个元素都是素数子群生成元的倍数

        let basis_num = 1;
        let mut basic_crs = Vec::with_capacity(basis_num);
        for i in 0..basis_num {
            basic_crs.push(Element::prime_subgroup_generator() * Fr::from((i + 1) as u64));
        }
        let mut scalars = vec![];
        //q-1
        scalars.push(Fr::from_str("13108968793781547619861935127046491459309155893440570251786403306729687672800").unwrap());


        let precompute=MSMPrecompWnaf::new(&basic_crs, 4);
        let mem_byte_size=precompute.tables.len()*precompute.tables[0].len()*4*32;
        println!("precompute_size: {:?}", mem_byte_size);
        use std::time::Instant;
        let start = Instant::now();
        for _i in 0..999 {
            let _got_result = precompute.mul(&scalars);
        }
        let got_result = precompute.mul(&scalars);

        let duration = start.elapsed();
        println!("Time elapsed in mul is: {:?}", duration/1000);

        let affine_result= got_result.0.into_affine();
        let string_x="33549696307925229982445904590536874618633472405590028303463218160177641247209";
        let string_y="19188667384257783945677642223292697773471335439753913231509108946878080696678";
        let x= affine_result.x.to_string();
        let y= affine_result.y.to_string();
        assert_eq!(string_x, x);
        assert_eq!(string_y, y);
        println!("got_result: {:?}", affine_result);

    }

    #[test]

    fn correctness_window_msm_multi_scalar_one_g() {
        // Create a vector of 256 elements, each being a multiple of the prime subgroup generator
        // 创建一个包含 256 个元素的向量，每个元素都是素数子群生成元的倍数

        let scalar_num = 1;
        let basis_num = 1;
        let mut basic_crs = Vec::with_capacity(basis_num);
        for i in 0..basis_num {
            basic_crs.push(Element::prime_subgroup_generator() * Fr::from((i + 1) as u64));
        }

        let mut rng = ChaCha20Rng::from_seed([2u8; 32]);
        let scalars: Vec<Fr> = (0..scalar_num).map(|_| Fr::rand(&mut rng)).collect();

        let precompute = MSMPrecompWnaf::new(&basic_crs, 4);
        let mem_byte_size = precompute.tables.len() * precompute.tables[0].len() * 4 * 32;
        println!("precompute_size: {:?}", mem_byte_size);
        use std::time::Instant;
        let start = Instant::now();
        for i in 0..scalar_num {

            let mut scalars0 = vec![];
            scalars0.push(scalars[i]);
            let _got_result = precompute.mul(&scalars0);

        }
        let duration = start.elapsed();
        println!("Time elapsed in mul is: {:?}", duration / scalar_num as u32);

    }
}
