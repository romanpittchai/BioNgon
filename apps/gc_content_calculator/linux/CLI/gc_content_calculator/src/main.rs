use std::env::{self};
use std::process::{self};

use gc_content_calculator::NucleotideCounter;

fn main() {
    let args: Vec<String> = env::args().collect();
    let config: NucleotideCounter = NucleotideCounter::build(&args).unwrap_or_else(|err: &str| {
        eprintln!("Problem parsing arguments: {err}");
        process::exit(1);
    });

    if let Err(e) = gc_content_calculator::run(config) {
        eprintln!("Error: {e}");
        process::exit(1);
    }
}
