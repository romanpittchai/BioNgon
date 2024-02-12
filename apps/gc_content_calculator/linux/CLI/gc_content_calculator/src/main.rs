use std::env::{self};
use std::process::{self};

use gc_content_calculator::NucleotideCounter;

fn main() {
    // Parse the command line arguments and create 
    // an instance of NucleotideCounter.
    let config: NucleotideCounter = NucleotideCounter::build(env::args())
        .unwrap_or_else(
            |err: &str| {
        eprintln!("Problem parsing arguments: {err}");
        process::exit(1);
    });

    if let Err(e) = gc_content_calculator::run(config) {
        eprintln!("Error: {e}");
        process::exit(1);
    }
}
