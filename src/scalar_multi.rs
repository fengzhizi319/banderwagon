use ark_ec::{PrimeGroup};
use ark_ff::{BigInteger, PrimeField, Zero};
use ark_std::vec::Vec;
use crate::element::Element;
use ark_ed_on_bls12_381_bandersnatch::{EdwardsProjective, Fq};
use rayon::prelude::*;
use crate::Fr;


#[derive(Clone, Debug)]
pub struct MSMExtendPrecompWnaf {
    window_size: usize,
    tables: Vec<Vec<ExtendPoint>>,
}

#[derive(Clone, Debug, Default)]
pub struct ExtendPoint {
    pub x: Fq,
    pub y: Fq,
    //pub t: Fq,
}

extern "C" {
    fn fr_add(res: &mut [u64; 4], a: &[u64; 4], b: &[u64; 4]);
    fn mul_by_5(res: &mut [u64; 4]);
    fn fr_mul(res: &mut [u64; 4], a: &[u64; 4], b: &[u64; 4]);
}

impl MSMExtendPrecompWnaf {
    pub fn new(bases: &[Element], window_size: usize) -> MSMExtendPrecompWnaf {
        let wnaf_context = WnafContext::new(window_size);

        let r = MSMExtendPrecompWnaf {
            tables: bases.par_iter().map(|base|{
                wnaf_context.extend_table(base.0)
            }).collect(),
            window_size,
        };

        r
    }

    pub fn mul_index(&self, scalar: Fr, index: usize) -> Element {
        let wnaf_context = WnafContext::new(self.window_size);
        Element(
            wnaf_context
                .mul_with_extend_table(&self.tables[index], &scalar),
        )
    }

    pub fn mul_index_inline(&self, scalar: Fr, index: usize) -> Element {
        let wnaf_context = WnafContext::new(self.window_size);
        Element(
            wnaf_context
                .mul_with_extend_table_inline(&self.tables[index], &scalar),
        )
    }

    pub fn fr_asm_mul(a: &[Fq], b: &[Fq]) -> EdwardsProjective  {
        let mut result = EdwardsProjective::zero();
        a.iter().zip(b.iter()).for_each(|(a, b)| {
            unsafe { fr_mul(&mut result.x.0.0, &a.0.0, &b.0.0) };
        });
        result
    }

    pub fn fr_asm_mul_inline(a: &[Fq], b: &[Fq]) -> EdwardsProjective  {
        let mut result = EdwardsProjective::zero();
        a.iter().zip(b.iter()).for_each(|(a, b)| {
            asm_mul(&mut result.x.0.0, &a.0.0, &b.0.0);
        });
        result
    }

    pub fn fr_mul_x(a: &[Fq], b: &[Fq]) -> EdwardsProjective  {
        let mut result = EdwardsProjective::zero();
        a.iter().zip(b.iter()).for_each(|(a, b)| {
            result.x = a * b;
        });
        result
    }

    pub fn fr_mul_iseq(a: &[Fq], b: &[Fq]) -> bool  {
        for i in 0..a.len() {
            let mut _cmp_b = Fq::default();
            let mut _cmp_c = Fq::default();
            let _cmp_a = a[i] * b[i];

            unsafe { fr_mul(&mut _cmp_b.0.0, &a[i].0.0, &b[i].0.0) };
            if _cmp_b != _cmp_a {
                println!("Error _cmp_b : {:?} != {:?}", _cmp_b.0.0, _cmp_a.0.0);
                return false;
            }
            asm_mul(&mut _cmp_c.0.0, &a[i].0.0, &b[i].0.0);
            if _cmp_c != _cmp_a {
                println!("Error _cmp_c : {:?} != {:?}", _cmp_c.0.0, _cmp_a.0.0);
                return false;
            }

        }
        true
    }

    pub fn ecc_asm_mul(&self, a: &[Fr]) -> Element {
        let mut res: Element = Element::zero();
        for i in 0..a.len() {
            res = self.mul_index(a[i], i%256);
        }
        res
    }

    pub fn fix_ecc_asm_mul(&self, _a: &[Fr]) -> Element {
        let mut res: Element = Element::zero();

        use std::str::FromStr;

        let scalars = vec![
            Fr::from_str("13108968793781547619861935127046491459309155893440570251786403306729687672800").unwrap()
        ];

        for _i in 0.._a.len() {
            res = self.mul_index(scalars[0], 0);
        }
        res
    }

    pub fn fixg_ecc_asm_mul(&self, _a: &[Fr]) -> Element {
        let mut res: Element = Element::zero();

        for _i in 0.._a.len() {
            res = self.mul_index(_a[_i], 0);
        }
        res
    }

    pub fn ecc_asm_mul_inline(&self, a: &[Fr]) -> Element {
        let mut res: Element = Element::zero();
        for i in 0..a.len() {
            res = self.mul_index_inline(a[i], i%256);
        }
        res
    }

    pub fn fix_ecc_asm_mul_inline(&self, _a: &[Fr]) -> Element {
        let mut res: Element = Element::zero();
        use std::str::FromStr;

        let scalars = vec![
            Fr::from_str("13108968793781547619861935127046491459309155893440570251786403306729687672800").unwrap()
        ];

        for _i in 0.._a.len() {
            res = self.mul_index_inline(scalars[0], 0);
        }
        res
    }

    pub fn fixg_ecc_asm_mul_inline(&self, _a: &[Fr]) -> Element {
        let mut res: Element = Element::zero();

        for _i in 0.._a.len() {
            res = self.mul_index_inline(_a[_i], 0);
        }
        res
    }

    pub fn asm_ffi_ecc_add(&self, l: usize) -> EdwardsProjective {
        let mut result = EdwardsProjective::zero();
        let vl = self.tables[0].len();
        for i in 0..l {
            extended_add_2d(&mut result, &self.tables[i%256][i%vl]);
        }
        result
    }

    pub fn asm_inline_ecc_add(&self, l: usize) -> EdwardsProjective {
        let mut result = EdwardsProjective::zero();
        let vl = self.tables[0].len();
        for i in 0..l {
            extended_add_2d_inline(&mut result, &self.tables[i%256][i%vl]);
        }
        result
    }

    /*pub fn mul(&self, scalars: &[Fr]) -> Element {
        let wnaf_context = WnafContext::new(self.window_size);
        let result: EdwardsProjective = scalars
            .iter()
            .zip(self.tables.iter())
            .filter(|(scalar, _)| !scalar.is_zero())
            .map(|(scalar, table)| wnaf_context.mul_with_table(table, scalar).unwrap())
            .sum();

        Element(result)
    }
    // TODO: This requires more benchmarking and feedback to see if we should
    // TODO put this behind a config flag
    pub fn mul_par(&self, scalars: &[Fr]) -> Element {
        let wnaf_context = WnafContext::new(self.window_size);
        let result: EdwardsProjective = scalars
            .par_iter()
            .zip(self.tables.par_iter())
            .filter(|(scalar, _)| !scalar.is_zero())
            .map(|(scalar, table)| wnaf_context.mul_with_table(table, scalar).unwrap())
            .sum();

        Element(result)
    }*/
}

/// A helper type that contains all the context required for computing
/// a window NAF multiplication of a group element by a scalar.
pub struct WnafContext {
    pub window_size: usize,
}

impl WnafContext {
    /// Constructs a new context for a window of size `window_size`.
    ///
    /// # Panics
    ///
    /// This function will panic if not `2 <= window_size < 64`
    pub fn new(window_size: usize) -> Self {
        assert!(window_size >= 2);
        assert!(window_size < 64);
        Self { window_size }
    }

    pub fn table<G: PrimeGroup>(&self, mut base: G) -> Vec<G> {
        let win_num = (256 + self.window_size - 1) / self.window_size;
        let mut table = Vec::with_capacity(win_num * (1 << (self.window_size - 1)));
        let half_size = 1 << (self.window_size - 1);

        for _ in 0..win_num {
            let dbl = base;
            table.push(base);
            for _i in 1..half_size {
                base += &dbl;
                table.push(base);
            }
            base += base;
        }

        table
    }

    pub fn extend_table(&self, base: EdwardsProjective) -> Vec<ExtendPoint> {
        Element::batch_extended_point_normalized(&self.table(base))
    }

    /// Computes scalar multiplication of a group element `g` by `scalar`.
    ///
    /// This method uses the wNAF algorithm to perform the scalar
    /// multiplication; first, it uses `Self::table` to calculate an
    /// appropriate table of multiples of `g`, and then uses the wNAF
    /// algorithm to compute the scalar multiple.
    pub fn mul<G: PrimeGroup>(&self, g: G, scalar: &G::ScalarField) -> G {
        let table = self.table(g);
        self.mul_with_table(&table, scalar).unwrap()
    }

    pub fn mul_with_extend_table(&self, base_table: &[ExtendPoint], scalar: &Fr) -> EdwardsProjective {
        let data = WnafContext::to_extend_scalar_data(scalar, self.window_size);

        let mut c = 0;
        let half_size = 1<<(self.window_size-1);
        let mut result = EdwardsProjective::zero();

        for i in 0..data.len() {
            let mut idx = (data[i] + c) as usize;
            if idx == 0 { continue; }

            c = 0;
            if idx > half_size {
                idx = (1 << self.window_size) - idx;
                if idx != 0 {
                    let neg_point = ExtendPoint {
                        x: -base_table[idx-1+i*half_size].x,
                        y: base_table[idx-1+i*half_size].y,
                        //t: -base_table[idx-1+i*half_size].t,
                    };
                    extended_add_2d(&mut result, &neg_point);
                }
                c = 1;
            } else {
                extended_add_2d(&mut result, &base_table[idx-1+i*half_size]);
            }
        }
        result
    }


    pub fn mul_with_extend_table_inline(&self, base_table: &[ExtendPoint], scalar: &Fr) -> EdwardsProjective {
        let data = WnafContext::to_extend_scalar_data(scalar, self.window_size);

        let mut c = 0;
        let half_size = 1<<(self.window_size-1);
        let mut result = EdwardsProjective::zero();

        for i in 0..data.len() {
            let mut idx = (data[i] + c) as usize;
            if idx == 0 { continue; }

            c = 0;
            if idx > half_size {
                idx = (1 << self.window_size) - idx;
                if idx != 0 {
                    let neg_point = ExtendPoint {
                        x: -base_table[idx-1+i*half_size].x,
                        y: base_table[idx-1+i*half_size].y,
                        //t: -base_table[idx-1+i*half_size].t,
                    };
                    extended_add_2d_inline(&mut result, &neg_point);
                }
                c = 1;
            } else {
                extended_add_2d_inline(&mut result, &base_table[idx-1+i*half_size]);
            }
        }
        result
    }

    /// Computes scalar multiplication of a group element by `scalar`.
    /// `base_table` holds precomputed multiples of the group element; it can be
    /// generated using `Self::table`. `scalar` is an element of
    /// `G::ScalarField`.
    ///
    /// Returns `None` if the table is too small.

    pub fn mul_with_table<G: PrimeGroup>(&self, base_table: &[G], scalar: &G::ScalarField) -> Option<G> {
        if 1 << (self.window_size - 1) > base_table.len() {
            return None;
        }
        // The modification principle comes from go-ipa(https://github.com/crate-crypto/go-ipa)
        let data = WnafContext::to_scalar_data::<G>(scalar, self.window_size);

        let mut c = 0;
        let thr = 1<<(self.window_size-1);
        let mut result = G::zero();

        for i in 0..data.len() {
            let mut idx = (data[i] + c) as usize;
            if idx == 0 { continue; }

            c = 0;
            if idx > thr {
                idx = (1 << self.window_size) - idx;
                if idx != 0 {
                    result -= &base_table[idx-1+i*thr];
                }
                c = 1;
            } else {
                result += &base_table[idx-1+i*thr];
            }
        }
        Some(result)
    }

    #[inline]
    fn to_scalar_data<G: PrimeGroup>(scalar: &G::ScalarField, w: usize) -> Vec<u64> {
        let source = WnafContext::to_u64::<G>(scalar);
        let mask = (1<<w) -1;
        let mut data = vec![];
        let mut off = w;

        for i in 0..source.len() {
            let s = if off != w {
                let mask = (1<<(w-off)) - 1;
                let j = data.len() - 1;
                data[j] += (source[i] & mask) << off;
                (source[i] >> (w - off), 64 - w + off)
            } else {
                (source[i], 64)
            };

            for j in (0..s.1).step_by(w) {
                let d = (s.0>>j)& mask;
                data.push(d);
                off = j;
            }
            off = s.1 - off;
        }


        data
    }

    #[inline]
    fn to_extend_scalar_data(scalar: &Fr, w: usize) -> Vec<u64> {
        let source = scalar.into_bigint().0;
        let mask = (1<<w) -1;
        let mut data = vec![];
        let mut off = w;

        for i in 0..source.len() {
            let s = if off != w {
                let mask = (1<<(w-off)) - 1;
                let j = data.len() - 1;
                data[j] += (source[i] & mask) << off;
                (source[i] >> (w - off), 64 - w + off)
            } else {
                (source[i], 64)
            };

            for j in (0..s.1).step_by(w) {
                let d = (s.0>>j)& mask;
                data.push(d);
                off = j;
            }
            off = s.1 - off;
        }


        data
    }

    #[inline]
    fn to_u64<G: PrimeGroup>(scalar: &G::ScalarField) -> Vec<u64> {
        let b = scalar.into_bigint();
        let mut num = b.num_bits();
        num = if num & 63 == 0{
            num >> 6
        } else {
            (num >> 6) + 1
        };

        let mut res = vec![0u64; num as usize];
        for i in 0..res.len() {
            res[i] = b.as_ref()[i];
        }
        res
    }
}

pub fn extended_neg(p: &ExtendPoint) -> ExtendPoint {
    ExtendPoint {
        x: -p.x,
        y: p.y,
        //t: -p.t,
    }
}

pub fn extended_add_2d_inline(result: &mut EdwardsProjective, p2: &ExtendPoint) {
    let mut a = Fq::default();
    let mut b = Fq::default();
    let mut c = Fq::default();
    let mut tmp = Fq::default();

    asm_mul(&mut c.0.0, &p2.x.0.0, &p2.y.0.0);

    asm_mul(&mut a.0.0, &result.x.0.0, &p2.x.0.0);
    asm_mul(&mut b.0.0, &result.y.0.0, &p2.y.0.0);
    //asm_mul(&mut tmp.0.0, &result.t.0.0, &p2.t.0.0);
    asm_mul(&mut tmp.0.0, &result.t.0.0, &c.0.0);

    asm_mul(&mut c.0.0, &tmp.0.0, &[12167860994669987632u64, 4043113551995129031u64, 6052647550941614584u64, 3904213385886034240u64]);
    asm_mul(&mut tmp.0.0, &(result.x + result.y).0.0, &(p2.x + p2.y).0.0);
    let e = tmp - a - b;
    let f = result.z - c;
    let g = result.z + c;
    unsafe { mul_by_5(&mut a.0.0) };
    let h = b + a;
    asm_mul(&mut result.x.0.0, &e.0.0, &f.0.0);
    asm_mul(&mut result.y.0.0, &g.0.0, &h.0.0);
    asm_mul(&mut result.t.0.0, &e.0.0, &h.0.0);
    asm_mul(&mut result.z.0.0, &f.0.0, &g.0.0);
}

pub fn extended_add_2d(result: &mut EdwardsProjective, p2: &ExtendPoint) {
    let mut a = Fq::default();
    let mut b = Fq::default();
    let mut c = Fq::default();
    let mut tmp = Fq::default();
    //calculate t
    unsafe { fr_mul(&mut c.0.0, &p2.x.0.0, &p2.y.0.0) };
    //let a = result.x * p2.x;
    unsafe { fr_mul(&mut a.0.0, &result.x.0.0, &p2.x.0.0) };

    //let b = result.y * p2.y;
    unsafe { fr_mul(&mut b.0.0, &result.y.0.0, &p2.y.0.0) };

    //let c = result.t * p2.t;
    //unsafe { fr_mul(&mut tmp.0.0, &result.t.0.0, &p2.t.0.0) };
    unsafe { fr_mul(&mut tmp.0.0, &result.t.0.0, &c.0.0) };

    //let mut d_c = Fq::default();
    //d_c.0.0 = [12167860994669987632u64, 4043113551995129031u64, 6052647550941614584u64, 3904213385886034240u64];
    //let c = c * d_c;

    unsafe { fr_mul(&mut c.0.0, &tmp.0.0, &[12167860994669987632u64, 4043113551995129031u64, 6052647550941614584u64, 3904213385886034240u64]) };

    //let d = result.z;

    //let tmp = result.x + result.y;
    //let mut e = p2.x + p2.y;
    //e = e * tmp - a - b;
    unsafe { fr_mul(&mut tmp.0.0, &(result.x + result.y).0.0, &(p2.x + p2.y).0.0) };
    let e = tmp - a - b;

    let f = result.z - c;

    let g = result.z + c;

    //let h = -a;
    //let h_mul_by_5 = h * Fq::from(5u64);
    //let h = b - h_mul_by_5;
    unsafe { mul_by_5(&mut a.0.0) };

    let h = b + a;

    //result.x = e * f;
    unsafe { fr_mul(&mut result.x.0.0, &e.0.0, &f.0.0) };
    //println!("ax {:?}", result.x.0.0); //ax [11263135436132592868, 16380222789396204840, 265672958677140400, 1080725393601644686]
    //result.y = g * h;
    unsafe { fr_mul(&mut result.y.0.0, &g.0.0, &h.0.0) };
    //println!("ay {:?}", result.y.0.0); //ay [16929979193753380128, 11289447262915390539, 6404279535991601865, 6514005977447446075]
    //result.t = e * h;
    unsafe { fr_mul(&mut result.t.0.0, &e.0.0, &h.0.0) };
    //println!("at {:?}", result.t.0.0); //at [8041554257280657769, 17492341825778612391, 15317427848265567591, 2085175354211983811]
    //result.z = f * g;
    unsafe { fr_mul(&mut result.z.0.0, &f.0.0, &g.0.0) };
    //println!("az {:?}", result.z.0.0); //az [2924513684221029011, 12405931109110048374, 3236999531136880018, 5137933263988462555]

}


use std::arch::asm;

#[derive(Debug, Clone, Copy)]
#[repr(C)]
struct Uint256 {
    data: [u64; 4],
}

const Q: [u64; 4] = [
    0xffffffff00000001,
    0x53bda402fffe5bfe,
    0x3339d80809a1d805,
    0x73eda753299d7d48,
];

const Q_INV: u64 = 0xfffffffeffffffff;

#[inline(always)]
pub fn asm_mul(result: &mut [u64; 4], x: &[u64; 4], y: &[u64; 4]) {

    unsafe {
        asm!(
            "push %r15",
            "push %r14",
            "push %r13",
            "push %r12",
            "push %rbp",
            "push %rcx",
            "push %rbx",

            // x -> rsi, y -> rdx, z -> rdi
            "movq %rdi, %r12",     // save res ptr
            "movq 0(%rsi), %rdi",
            "movq 8(%rsi), %r8",
            "movq 16(%rsi), %r9",
            "movq 24(%rsi), %r10",

            "movq %r12, %rsi",     // restore res ptr
            "movq %rdx, %r11",

            // clear the flags
            "xorq %rax, %rax",
            "movq 0(%r11), %rdx",

            // (A,t[0])  := x[0]*y[0] + A
            "mulxq %rdi, %r14, %r13",

            // (A,t[1])  := x[1]*y[0] + A
            "mulxq %r8, %rax, %rcx",
            "adoxq %rax, %r13",

            // (A,t[2])  := x[2]*y[0] + A
            "mulxq %r9, %rax, %rbx",
            "adoxq %rax, %rcx",

            // (A,t[3])  := x[3]*y[0] + A
            "mulxq %r10, %rax, %rbp",
            "adoxq %rax, %rbx",

            // A += carries from ADCXQ and ADOXQ
            "movq $0, %rax",
            "adoxq %rax, %rbp",

            // m := t[0]*q'[0] mod W
            "movq ${q_inv}, %rdx",
            "imulq %r14, %rdx",

            // clear the flags
            "xorq %rax, %rax",

            // C,_ := t[0] + m*q[0]
            "movq ${q0}, %r15",
            "mulxq %r15, %rax, %r12",
            "adcxq %r14, %rax",
            "movq %r12, %r14",

            // (C,t[0]) := t[1] + m*q[1] + C
            "adcxq %r13, %r14",
            "movq ${q1}, %r15",
            "mulxq %r15, %rax, %r13",
            "adoxq %rax, %r14",

            // (C,t[1]) := t[2] + m*q[2] + C
            "adcxq %rcx, %r13",
            "movq ${q2}, %r15",
            "mulxq %r15, %rax, %rcx",
            "adoxq %rax, %r13",

            // (C,t[2]) := t[3] + m*q[3] + C
            "adcxq %rbx, %rcx",
            "movq ${q3}, %r15",
            "mulxq %r15, %rax, %rbx",
            "adoxq %rax, %rcx",

            // t[3] = C + A
            "movq $0, %rax",
            "adcxq %rax, %rbx",
            "adoxq %rbp, %rbx",

            // clear the flags
            "xorq %rax, %rax",
            "movq 8(%r11), %rdx",

            // (A,t[0])  := t[0] + x[0]*y[1] + A
            "mulxq %rdi, %rax, %rbp",
            "adoxq %rax, %r14",

            // (A,t[1])  := t[1] + x[1]*y[1] + A
            "adcxq %rbp, %r13",
            "mulxq %r8, %rax, %rbp",
            "adoxq %rax, %r13",

            // (A,t[2])  := t[2] + x[2]*y[1] + A
            "adcxq %rbp, %rcx",
            "mulxq %r9, %rax, %rbp",
            "adoxq %rax, %rcx",

            // (A,t[3])  := t[3] + x[3]*y[1] + A
            "adcxq %rbp, %rbx",
            "mulxq %r10, %rax, %rbp",
            "adoxq %rax, %rbx",

            // A += carries from ADCXQ and ADOXQ
            "movq $0, %rax",
            "adcxq %rax, %rbp",
            "adoxq %rax, %rbp",

            // m := t[0]*q'[0] mod W
            "movq ${q_inv}, %rdx",
            "imulq %r14, %rdx",
            // clear the flags
            "xorq %rax, %rax",

            // C,_ := t[0] + m*q[0]
            "movq ${q0}, %r15",
            "mulxq %r15, %rax, %r12",
            "adcxq %r14, %rax",
            "movq %r12, %r14",

            // (C,t[0]) := t[1] + m*q[1] + C
            "adcxq %r13, %r14",
            "movq ${q1}, %r15",
            "mulxq %r15, %rax, %r13",
            "adoxq %rax, %r14",

            // (C,t[1]) := t[2] + m*q[2] + C
            "adcxq %rcx, %r13",
            "movq ${q2}, %r15",
            "mulxq %r15, %rax, %rcx",
            "adoxq %rax, %r13",

            // (C,t[2]) := t[3] + m*q[3] + C
            "adcxq %rbx, %rcx",
            "movq ${q3}, %r15",
            "mulxq %r15, %rax, %rbx",
            "adoxq %rax, %rcx",

            // t[3] = C + A
            "movq $0, %rax",
            "adcxq %rax, %rbx",
            "adoxq %rbp, %rbx",

            // clear the flags
            "xorq %rax, %rax",
            "movq 16(%r11), %rdx",

            // (A,t[0])  := t[0] + x[0]*y[2] + A
            "mulxq %rdi, %rax, %rbp",
            "adoxq %rax, %r14",

            // (A,t[1])  := t[1] + x[1]*y[2] + A
            "adcxq %rbp, %r13",
            "mulxq %r8, %rax, %rbp",
            "adoxq %rax, %r13",

            // (A,t[2])  := t[2] + x[2]*y[2] + A
            "adcxq %rbp, %rcx",
            "mulxq %r9, %rax, %rbp",
            "adoxq %rax, %rcx",

            // (A,t[3])  := t[3] + x[3]*y[2] + A
            "adcxq %rbp, %rbx",
            "mulxq %r10, %rax, %rbp",
            "adoxq %rax, %rbx",

            // A += carries from ADCXQ and ADOXQ
            "movq  $0, %rax",
            "adcxq %rax, %rbp",
            "adoxq %rax, %rbp",

            // m := t[0]*q'[0] mod W
            "movq ${q_inv}, %rdx",
            "imulq %r14, %rdx",

            // clear the flags
            "xorq %rax, %rax",

            // C,_ := t[0] + m*q[0]
            "movq ${q0}, %r15",
            "mulxq %r15, %rax, %r12",
            "adcxq %r14, %rax",
            "movq  %r12, %r14",

            // (C,t[0]) := t[1] + m*q[1] + C
            "adcxq %r13, %r14",
            "movq ${q1}, %r15",
            "mulxq %r15, %rax, %r13",
            "adoxq %rax, %r14",

            // (C,t[1]) := t[2] + m*q[2] + C
            "adcxq %rcx, %r13",
            "movq ${q2}, %r15",
            "mulxq %r15, %rax, %rcx",
            "adoxq %rax, %r13",

            // (C,t[2]) := t[3] + m*q[3] + C
            "adcxq %rbx, %rcx",
            "movq ${q3}, %r15",
            "mulxq %r15, %rax, %rbx",
            "adoxq %rax, %rcx",

            // t[3] = C + A
            "movq  $0, %rax",
            "adcxq %rax, %rbx",
            "adoxq %rbp, %rbx",

            // clear the flags
            "xorq %rax, %rax",
            "movq 24(%r11), %rdx",

            // (A,t[0])  := t[0] + x[0]*y[3] + A
            "mulxq %rdi, %rax, %rbp",
            "adoxq %rax, %r14",

            // (A,t[1])  := t[1] + x[1]*y[3] + A
            "adcxq %rbp, %r13",
            "mulxq %r8, %rax, %rbp",
            "adoxq %rax, %r13",

            // (A,t[2])  := t[2] + x[2]*y[3] + A
            "adcxq %rbp, %rcx",
            "mulxq %r9, %rax, %rbp",
            "adoxq %rax, %rcx",

            // (A,t[3])  := t[3] + x[3]*y[3] + A
            "adcxq %rbp, %rbx",
            "mulxq %r10, %rax, %rbp",
            "adoxq %rax, %rbx",

            // A += carries from ADCXQ and ADOXQ
            "movq  $0, %rax",
            "adcxq %rax, %rbp",
            "adoxq %rax, %rbp",

            // m := t[0]*q'[0] mod W
            "movq ${q_inv}, %rdx",
            "imulq %r14, %rdx",

            // clear the flags
            "xorq %rax, %rax",

            // C,_ := t[0] + m*q[0]
            "movq ${q0}, %r15",
            "mulxq %r15, %rax, %r12",
            "adcxq %r14, %rax",
            "movq  %r12, %r14",

            // (C,t[0]) := t[1] + m*q[1] + C
            "adcxq %r13, %r14",
            "movq ${q1}, %r15",
            "mulxq %r15, %rax, %r13",
            "adoxq %rax, %r14",

            // (C,t[1]) := t[2] + m*q[2] + C
            "adcxq %rcx, %r13",
            "movq ${q2}, %r15",
            "mulxq %r15, %rax, %rcx",
            "adoxq %rax, %r13",

            // (C,t[2]) := t[3] + m*q[3] + C
            "adcxq %rbx, %rcx",
            "movq ${q3}, %r15",
            "mulxq %r15, %rax, %rbx",
            "adoxq %rax, %rcx",

            // t[3] = C + A
            "movq  $0, %rax",
            "adcxq %rax, %rbx",
            "adoxq %rbp, %rbx",

            // reduce element(R14,R13,CX,BX) using temp registers (SI,R12,R11,DI)
            "movq    %rsi, %rax",
            "movq    %r14, %rsi",
            "movq    ${q0}, %r15",
            "subq    %r15, %r14",
            "movq    %r13, %r12",
            "movq    ${q1}, %r15",
            "sbbq    %r15, %r13",
            "movq    %rcx, %r11",
            "movq    ${q2}, %r15",
            "sbbq    %r15, %rcx",
            "movq    %rbx, %rdi",
            "movq    ${q3}, %r15",
            "sbbq    %r15, %rbx",
            "cmovc   %rsi, %r14",
            "cmovc   %r12, %r13",
            "cmovc   %r11, %rcx",
            "cmovc   %rdi, %rbx",

            "movq   %r14, 0(%rax)",
            "movq   %r13, 8(%rax)",
            "movq   %rcx, 16(%rax)",
            "movq   %rbx, 24(%rax)",


            "pop %rbx",
            "pop %rcx",
            "pop %rbp",
            "pop %r12",
            "pop %r13",
            "pop %r14",
            "pop %r15",
            inout("rsi") x.as_ptr() => _,
            inout("rdx") y.as_ptr() => _,
            inout("rdi") result.as_mut_ptr() => _,
            q0 = const Q[0],
            q1 = const Q[1],
            q2 = const Q[2],
            q3 = const Q[3],
            q_inv = const Q_INV,
            out("rax") _,
            // out("rcx") _,
            out("r8") _,
            out("r9") _,
            out("r10") _,
            out("r11") _,

            options(att_syntax),
            clobber_abi("C")
        );
    }

}

pub struct WnafGottiContext {
    pub t: usize,
    pub b: usize,
}


impl WnafGottiContext {
    pub fn new(t: usize,b: usize) -> Self {
        assert!(b >= 2);
        assert!(b < 64);
        Self { t,b }
    }
    fn fill_window<G: PrimeGroup>(&self, basis: &mut [G], table: &mut [G]) {
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
        self.fill_window(&mut basis[1..], &mut table[.. l/ 2]);
        // 遍历 table 数组的前半部分
        for i in 0..l / 2 {
            // 将 table 数组后半部分的对应元素设置为当前元素与 basis 数组第一个元素的和
            table[l / 2 + i] = table[i] + basis[0];
        }
    }



    pub fn mul<G: PrimeGroup>(&self, g: G, scalar: &G::ScalarField) -> G {
        let table = self.table(g);
        self.mul_with_table(&table, scalar).unwrap()
    }


    fn scalar_to_u64<G: PrimeGroup>(scalar: &G::ScalarField) -> Vec<u64> {
        let b = scalar.into_bigint();
        let mut num = b.num_bits();
        num = if num & 63 == 0 {
            num >> 6
        } else {
            (num >> 6) + 1
        };

        let mut res = Vec::with_capacity(num as usize);
        res.extend_from_slice(&b.as_ref()[..num as usize]);
        res
    }
    
    fn fr_to_u64(scalar: &Fr) -> Vec<u64> {
        let b = scalar.into_bigint();
        let mut num = b.num_bits();
        num = if num & 63 == 0 {
            num >> 6
        } else {
            (num >> 6) + 1
        };

        let mut res = Vec::with_capacity(num as usize);
        res.extend_from_slice(&b.as_ref()[..num as usize]);
        res
    }
    
    pub(crate) fn gotti_window<G: PrimeGroup>(&self, scalar: &G::ScalarField) -> Vec<Vec<u16>> {
        let b = scalar.into_bigint();
        let mut num = b.num_bits();
        num = if num & 63 == 0 {
            num >> 6
        } else {
            (num >> 6) + 1
        };
        let fr_bits = 253;
        let mut naf_vec: Vec<Vec<u8>> = vec![Vec::new(); self.t];

        let mut scalar_u64_4_ = Vec::with_capacity(num as usize);
        scalar_u64_4_.extend_from_slice(&b.as_ref()[..num as usize]);
        for t_i in 0..self.t {
            for k in (0..fr_bits).step_by(self.t) {
                let scalar_bit_pos = k + self.t - t_i - 1;
                // println!("k: {:?}", k);
                // println!("scalar_bit_pos: {:?}", scalar_bit_pos);
                if scalar_bit_pos < fr_bits && !scalar.is_zero() {
                    let limb = scalar_u64_4_[scalar_bit_pos >> 6];
                    let bit = ((limb >> (scalar_bit_pos & 63)) & 1) as u8;
                    naf_vec[t_i].push(bit);
                }
            }
        }
        let mut new_naf_vec: Vec<Vec<u16>> = Vec::with_capacity(self.t);

        for vec in naf_vec.iter() {
            let mut new_vec = Vec::with_capacity((vec.len() + self.b - 1) / self.b);
            let mut combined: u16 = 0;
            let mut bit_pos = 0;

            for &bit in vec.iter() {
                combined |= (bit as u16) << bit_pos;
                bit_pos += 1;

                if bit_pos == self.b {
                    new_vec.push(combined);
                    combined = 0;
                    bit_pos = 0;
                }
            }

            if bit_pos > 0 {
                new_vec.push(combined);
            }

            new_naf_vec.push(new_vec);
        }
        println!("new_naf_vec: {:?}", new_naf_vec);

        return new_naf_vec;
    }

    fn reverse_bits(&self,mut num: u32, t: usize) -> u32 {
        let mut reversed = 0;
        for _ in 0..t {
            reversed <<= 1;
            reversed |= num & 1;
            num >>= 1;
        }
        reversed
    }


    pub  fn rearrange_table<G: Clone>(&self,table: &[G]) -> Vec<G> {
        let mut table1 = vec![table[0].clone(); table.len()];
        for i in 0..table.len() {
            let reversed_index = self.reverse_bits(i as u32,self.b) as usize;
            table1[reversed_index] = table[i].clone();
        }
        table1
    }
    pub fn table<G: PrimeGroup>(&self,base: G)-> Vec<G> {

        let fr_bits =253;
        //b为窗口bit长度，t为预处理的倍数，window_size为每个窗口内的值的所有可能的值
        let window_size = 1 << self.b;
        //用t把标量分为t个部分，如果t=2，那么标量分为2部分，每部分的位数为fr_bits/2，kP=2A1+A0，
        // A1为高位，A0为低位，A1，A0即为划分出来的部分。A1、A0的所有可能的个数为points_per_column，
        // 如t=2，b=5，那么points_per_column=127
        let points_per_column = (fr_bits + self.t - 1) / self.t as usize;
        //窗口的个数
        let window_count = (points_per_column + self.b - 1) / self.b; // 等价于 ceil(256 / self.window_size)
        //let table_size = window_count * (1 << (self.b - 1));
        //let threshold = 1 << (self.b - 1);
        let mut table_basis = Vec::with_capacity(points_per_column);
        let mut dbl = base.clone();
        table_basis.push(base);
        //计算2^0*G,2^2*G,2^4*G,2^6*G,...,2^252*G
        for _i in 1..(points_per_column) {
            for _j in 0..self.t {
                dbl= dbl.double();
            }
            table_basis.push(dbl.clone());
        }
        let mut nn_table = vec![G::zero(); window_count * window_size];
        for i in 0..window_count {
            let w=i;
            let start = w * self.b as usize;
            let mut end = (w + 1) * self.b as usize;
            if end > table_basis.len() {
                end = table_basis.len();
            }
            let window_basis: &mut [G] = &mut table_basis[start..end];
            let mut table = vec![G::zero(); window_size];
            self.fill_window(window_basis, &mut *table);
            let  new_table=self.rearrange_table(&table);
            for (j, table_element) in new_table.into_iter().enumerate() {
                nn_table[w * window_size + j] = table_element;
            }
        }
        nn_table
    }
    pub fn mul_with_table<G: PrimeGroup>(&self, base_pre_table: &[G], mon_scalar: &G::ScalarField) -> Option<G> {
        if 1 << (self.b - 1) > base_pre_table.len() {
            return None;
        }

        let window_size = 1 << self.b;
        let mut accum = G::zero();
        let fr_bits = 253;
        let scalar_u64 = WnafGottiContext::scalar_to_u64::<G>(mon_scalar);

        for t_i in 0..self.t {
            if t_i > 0 {
                accum = accum.double();
            }

            let mut curr_window = 0;
            let mut window_scalar = 0;
            let mut window_bit_pos = 0;

            for k in (0..fr_bits).step_by(self.t) {
                let scalar_bit_pos = k + self.t - t_i - 1;

                if scalar_bit_pos < fr_bits && !mon_scalar.is_zero() {
                    let limb = scalar_u64[scalar_bit_pos >> 6];
                    let bit = (limb >> (scalar_bit_pos & 63)) & 1;
                    window_scalar |= (bit as usize) << (window_bit_pos);
                }

                window_bit_pos += 1;
                if window_bit_pos == self.b {
                    if window_scalar > 0 {
                        let window_precomp = &base_pre_table[curr_window * window_size..(curr_window + 1) * window_size];
                        accum += window_precomp[window_scalar];
                    }
                    curr_window += 1;
                    window_scalar = 0;
                    window_bit_pos = 0;
                }
            }

            if window_scalar > 0 {
                let window_slice = &base_pre_table[curr_window * window_size..(curr_window + 1) * window_size];
                accum += window_slice[window_scalar];
            }
        }

        Some(accum)
    }


    pub fn mul_with_table_2dasm(&self, base_pre_table: &[ExtendPoint], mon_scalar: &Fr) -> Option<EdwardsProjective> {
        if 1 << (self.b - 1) > base_pre_table.len() {
            return None;
        }

        let window_size = 1 << self.b;
        let mut accum = EdwardsProjective::zero();
        let fr_bits = 253;
        let scalar_u64 = WnafGottiContext::fr_to_u64(mon_scalar);

        for t_i in 0..self.t {
            if t_i > 0 {
                accum += accum;
            }

            let mut curr_window = 0;
            let mut window_scalar = 0;
            let mut window_bit_pos = 0;

            for k in (0..fr_bits).step_by(self.t) {
                let scalar_bit_pos = k + self.t - t_i - 1;

                if scalar_bit_pos < fr_bits && !mon_scalar.is_zero() {
                    let limb = scalar_u64[scalar_bit_pos >> 6];
                    let bit = (limb >> (scalar_bit_pos & 63)) & 1;
                    window_scalar |= (bit as usize) << (window_bit_pos);
                }

                window_bit_pos += 1;
                if window_bit_pos == self.b {
                    if window_scalar > 0 {
                        let window_precomp = &base_pre_table[curr_window * window_size..(curr_window + 1) * window_size];
                        //accum += window_precomp[window_scalar];
                        
                        extended_add_2d(&mut accum, &window_precomp[window_scalar]);
                        
                    }
                    curr_window += 1;
                    window_scalar = 0;
                    window_bit_pos = 0;
                }
            }

            if window_scalar > 0 {
                let window_slice = &base_pre_table[curr_window * window_size..(curr_window + 1) * window_size];
                //accum += window_slice[window_scalar];
                extended_add_2d(&mut accum, &window_slice[window_scalar]);
            }
        }

        Some(accum)
    }



}


#[cfg(test)]
mod tests {
    use super::*;
    use crate::{try_reduce_to_element, msm_window::MSMPrecompWnaf, Element, Fr};

    use ark_std::rand::SeedableRng;
    use ark_std::UniformRand;
    use std::time::Instant;
    use rand_chacha::ChaCha20Rng;

    fn generate_random_elements(num_required_points: usize, seed: &'static [u8]) -> Vec<Element> {
        use sha2::{Digest, Sha256};

        let _choose_largest = false;

        // Hash the seed + i to get a possible x value
        let hash_to_x = |index: u64| -> Vec<u8> {
            let mut hasher = Sha256::new();
            hasher.update(seed);
            hasher.update(index.to_be_bytes());
            let bytes: Vec<u8> = hasher.finalize().to_vec();
            bytes
        };

        (0u64..)
            .map(hash_to_x)
            .filter_map(|hash_bytes| try_reduce_to_element(&hash_bytes))
            .take(num_required_points)
            .collect()
    }

    #[test]
    fn mertic_msm() {
        let g_base: Vec<_> = generate_random_elements(256, b"eth_verkle_oct_2021").into_iter().collect();

        let precomp = MSMPrecompWnaf::new(&g_base, 10);

        for num_loads in [32] {
            //let mut rng = rand::thread_rng();
            let mut rng = ChaCha20Rng::from_seed([0u8; 32]);
            /*let ax_array: Vec<Fr> = (0..256)
                .map(|i| {
                    if i < num_loads {
                        Fr::rand(&mut rng)
                    } else {
                        Fr::from(0)
                    }
                })
                .collect();

            let start = Instant::now();
            msm_by_array(&committer, &ax_array);
            let duration = start.elapsed();
            println!("msm_by_array({}) took: {:?}", num_loads, duration);*/

            let ax_array: Vec<(Fr, usize)> = (0..num_loads).map(|i| (Fr::rand(&mut rng), i)).collect();

            msm_by_delta(&precomp, &ax_array);

        }
    }

    fn msm_by_delta(c: &MSMPrecompWnaf, frs: &[(Fr, usize)]) -> Element {
        let start = Instant::now();
        let e = commit_sparse(c, frs.to_vec());
        let duration = start.elapsed();
        println!("msm_by_delta({}) took: {:?}", frs.len(), duration);
        e
    }

    fn commit_sparse(precomp: &MSMPrecompWnaf, val_indices: Vec<(Fr, usize)>) -> Element {
        let mut result = Element::zero();

        for (value, lagrange_index) in val_indices {
            result += precomp.mul_index(value, lagrange_index) //self.scalar_mul(value, lagrange_index)
        }
        result
    }
}
#[test]
fn correctness_precompute_one_table_test() {
    // Create a vector of 256 elements, each being a multiple of the prime subgroup generator
    // 创建一个包含 256 个元素的向量，每个元素都是素数子群生成元的倍数
    let mut basic_g = Vec::with_capacity(1);

    basic_g.push(Element::prime_subgroup_generator());
    let wnaf_gotti_context = WnafGottiContext::new(2,8);

    // Perform the operation
    let precompute = wnaf_gotti_context.table(basic_g[0].0);
    println!("precompute len: {:?}", precompute.len());

}
