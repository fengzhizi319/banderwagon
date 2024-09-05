
use ark_std::rand::SeedableRng;
use ark_std::UniformRand;
use banderwagon::{Element, Fr, MSMPrecompWnaf};
use std::time::Instant;
use ffi_interface::Context;
use ipa_multipoint::committer::Committer;
use rand_chacha::ChaCha20Rng;
//use timetrace_ffi::*;

fn main() {
    mertic_msm();
}

fn mertic_msm() {
    let committer = Context::new_with_window_size(10);

    //set_output_filename("timetrace_ffi_test.log");
    //set_keepoldevents(true);

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
        
        msm_by_delta(&committer, &ax_array);
        
    }
    //tt_print();
}

fn msm_by_array(c: &Context, frs: &[Fr]) -> Element {
    c.committer.commit_lagrange(frs)
}

fn msm_by_delta(c: &Context, frs: &[(Fr, usize)]) -> Element {
    let start = Instant::now();
    let e = c.committer.commit_sparse(frs.to_vec());
    let duration = start.elapsed();
    println!("msm_by_delta({}) took: {:?}", frs.len(), duration);
    e
}