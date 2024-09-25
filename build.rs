

fn main() {
    println!("cargo:rerun-if-changed=src/ecc_ops.s");
    cc::Build::new()
        .file("src/ecc_ops.s")
        .compile("libeccops.a");
}
