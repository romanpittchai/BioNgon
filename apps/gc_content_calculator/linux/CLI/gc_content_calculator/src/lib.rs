//! # GC Content Calculator
//! 
//! A library for counting nucleatides. 
//! It is used for files in the FASTA format.

use std::error::Error;
pub use self::creation_and_counting::{NucleotideCounter, NucleotideCountType};


pub mod creation_and_counting {
    //! # Creation and counting
    //! 
    //! The module that contains all the basic logic of working with calculations. 
    //! The _NucleotideCounter_ structure is defined in the module, 
    //! as well as several methods of this structure. 
    //! 1. The "_build_" - method is used to create an instance of the structure;
    //! 2. The "_read_nuacliotides_" - method is for reading a file and counting nucleatides;
    //! 3. The "_count_other_nucleatides_" - method is for calculating additional nucleotide data;
    //! 4. The "_calc_percentages_" - method is an auxiliary method for calculating percentages;
    //! 5. The "_percentage_of_nucleotides_and_output_" - method is used to calculate percentages 
    //! and output data to the console.
    use std::{
        collections::HashMap, error::Error, 
        fs::File, io::{BufRead, BufReader}, 
        path::Path, sync::{Arc, Mutex, MutexGuard}, 
        thread::{self, JoinHandle},
    };
    use chrono::{DateTime, Datelike, Local, Timelike};

    /// Enum representing different types of nucleotide count values.
    #[derive(Debug, Clone, PartialEq)]
    pub enum NucleotideCountType {
        UnInt(u64),
        Float(f64),
    }

    /// A structure representing the nucleotide occurrence counter. 
    /// It also stores the path to the file, the file name, 
    /// the file header and the file processing time.
    
    #[derive(Debug, Clone)]
    pub struct NucleotideCounter {
        /// _Path to the file being processed._
        pub file_path: String,
        /// _Name of the file._
        pub filename: String,
        /// _Header information extracted from the file._
        pub file_header: String,
        /// _Formatted date and time._
        pub formatted_datetime: String,
        /// _Counts of nucleotides in the file._
        pub nucleatides: HashMap<char, Vec<NucleotideCountType>>,
        /// _Additional nucleotide counts for specific calculations._
        pub other_calc_nucleatides: HashMap<String, Vec<NucleotideCountType>>,
    }

    impl NucleotideCounter {
        const NUCLEOTIDES_LIST: [char; 6] = ['A', 'C', 'G', 'T', 'U', 'O'];
        const OTHER_PARAMETERS_NUCLEATIDE_CONST: [&'static str; 3] = [
                "Total number of characters", "Number of ATGCU", 
                "Number of GC content in characters"
            ];
        const NUCLEATIDES_LIST_GC: [char; 2] = ['G', 'C'];
        const NUM_THREADS: usize = 2;
        const LEN_BUF: usize = 8;
        const LEN_BUF_KW: usize = 1024;

        pub fn build(mut args: impl Iterator<Item = String>) -> Result<NucleotideCounter, &'static str> {
            // Constructs a new `NucleotideCounter` instance based on command line arguments.
            args.next();

            let file_path: String = match args.next() {
                Some(arg) => arg,
                None => return Err("Didn't get a filepath string"),
            };
            let path: &Path = Path::new(&file_path);
            let filename: String = if let Some(file_name) = path.file_name() {
                file_name.to_str().unwrap_or_default().to_string()
            } else {
                String::from("Empty")
            };
            let file_path: String = file_path.to_string();
            
            let file_header: String = String::from("Empty");

            let current_datetime: DateTime<Local> = Local::now();
            let formatted_datetime: String = format!(
                "{:02}:{:02} - {:02}.{:02}.{}",
                current_datetime.hour(),
                current_datetime.minute(),
                current_datetime.day(),
                current_datetime.month(),
                current_datetime.year(),
            );

            let mut nucleatides: HashMap<char, Vec<NucleotideCountType>> = HashMap::new();

            for nucleatide in &Self::NUCLEOTIDES_LIST {
                nucleatides.insert(*nucleatide, vec![
                    NucleotideCountType::UnInt(0), 
                    NucleotideCountType::Float(0.0)
                    ]
                );
            }

            let mut other_calc_nucleatides: HashMap<String, Vec<NucleotideCountType>> = HashMap::new();

            for parameter in &Self::OTHER_PARAMETERS_NUCLEATIDE_CONST {
                other_calc_nucleatides.insert(parameter.to_string(), vec![
                    NucleotideCountType::UnInt(0), 
                    NucleotideCountType::Float(0.0)
                    ]
                );
            }

            Ok(NucleotideCounter { 
                file_path, filename, 
                file_header, formatted_datetime, 
                nucleatides, other_calc_nucleatides 
            })
        }

        pub fn read_nuacliotides(&mut self) -> Result<(), Box<dyn Error>> {
            // Reads nucleotides from a file and updates the nucleotide counts.
            let file: File = File::open(&self.file_path)?;
            let reader: BufReader<File> = BufReader::with_capacity(
                Self::LEN_BUF * Self::LEN_BUF_KW, file
            );

            let lines: Result<Vec<String>, std::io::Error> = reader.lines().collect();
            let lines: Vec<String> = lines.map_err(
                |e| {
                    eprintln!("Error reading lines: {}", e); 
                    e
                }
            )?;
            let chunk_size: usize = lines.len() / Self::NUM_THREADS;

            let nucleotides_counter: Arc<Mutex<Self>> = Arc::new(
                Mutex::new(self.clone())
            );
            let mut handles: Vec<JoinHandle<()>> = vec![];
            let nucleotide_count_int: Arc<Mutex<[u64; 6]>> = Arc::new(
                Mutex::new([0; 6])
            );

            for i in 0..Self::NUM_THREADS {
                let start: usize = i * chunk_size;
                let end: usize = if i == Self::NUM_THREADS - 1 {
                    lines.len()
                } else {
                    (i + 1) * chunk_size
                };

                let lines_chunk: Vec<String>= lines[start..end].to_vec();
                let struct_counter: Arc<Mutex<NucleotideCounter>> = Arc::clone(
                    &nucleotides_counter
                );
                let counter: Arc<Mutex<[u64; 6]>> = Arc::clone(
                    &nucleotide_count_int
                );
                
                let handle = thread::spawn(move || {
                    let mut struct_counter: MutexGuard<'_, NucleotideCounter> = struct_counter.lock().unwrap(); 
                    for line in lines_chunk {
                        let sequence: String = line.to_string();
                        if sequence.starts_with('>') {
                            struct_counter.file_header = sequence.chars().skip(1).collect();
                            continue;
                        }

                        for char in sequence.chars() {
                            let mut found_match: bool = false;
                            for i in 0..Self::NUCLEOTIDES_LIST.len() {
                                if char == Self::NUCLEOTIDES_LIST[i] {
                                    let mut counter: MutexGuard<'_, [u64; 6]> = counter.lock().unwrap();
                                    counter[i] += 1;
                                    found_match = true;
                                }
                            }
                            if !found_match {
                                let mut counter: MutexGuard<'_, [u64; 6]> = counter.lock().unwrap();
                                counter[5] += 1;
                            }
                        }
                    }
                });
                handles.push(handle);
            }
            for handle in handles {
                handle.join().unwrap();
            }

            let struct_counter: MutexGuard<'_, NucleotideCounter> = nucleotides_counter.lock().unwrap();
            self.file_header = struct_counter.file_header.clone();

            let counter: MutexGuard<'_, [u64; 6]> = nucleotide_count_int.lock().unwrap();
            for i in 0..Self::NUCLEOTIDES_LIST.len() {
                if let Some(count_vec) = self.nucleatides.get_mut(
                    &Self::NUCLEOTIDES_LIST[i]
                ) {
                    if let Some(NucleotideCountType::UnInt(value)) = count_vec.get_mut(0) {
                        *value = counter[i];
                    } else {
                        count_vec.push(NucleotideCountType::UnInt(1));
                    }
                }
            };
            Ok(())
        }

        pub fn count_other_nucleatides(&mut self) -> Result<(), Box<dyn Error>> {
            // Counts additional nucleotides and updates the counts.
            let mut total_number_of_characters: u64 = 0;
            for nucleotide in &self.nucleatides {
                if let Some(NucleotideCountType::UnInt(value)) = nucleotide.1.get(0) {
                    total_number_of_characters += value;
                }
            }

            if let Some(count_vec) = self.other_calc_nucleatides.get_mut(
                Self::OTHER_PARAMETERS_NUCLEATIDE_CONST[0]
            ) {
                if let Some(NucleotideCountType::UnInt(total_count)) = count_vec.get_mut(0) {
                    *total_count = total_number_of_characters;
                } else {
                    count_vec.push(NucleotideCountType::UnInt(total_number_of_characters));
                }
            }

            let mut other_nucleatides: u64 = 0;
            if let Some(count_vec) = self.nucleatides.get(
                &Self::NUCLEOTIDES_LIST[5]
            ) {
                if let Some(
                    NucleotideCountType::UnInt(other_nucleatides_count)
                ) = count_vec.get(0) {
                    other_nucleatides = *other_nucleatides_count;
                }
            }

            if let Some(count_vec) = self.other_calc_nucleatides.get_mut(
                Self::OTHER_PARAMETERS_NUCLEATIDE_CONST[1]
            ) {
                if let Some(NucleotideCountType::UnInt(atgcu)) = count_vec.get_mut(0) {
                    *atgcu = total_number_of_characters - other_nucleatides;
                }
            }
            
            let mut gc_sum: [u64; 2] = [0; 2];
            for i in 0..gc_sum.len() {
                if let Some(count_vec) = self.nucleatides.get(
                    &Self::NUCLEATIDES_LIST_GC[i]
                ) {
                    if let Some(NucleotideCountType::UnInt(gc)) = count_vec.get(0) {
                        gc_sum[i] = *gc;
                    } 
                }
            }

            if let Some(count_vec) = self.other_calc_nucleatides.get_mut(
                Self::OTHER_PARAMETERS_NUCLEATIDE_CONST[2]
            ) { if let Some(NucleotideCountType::UnInt(gc)) = count_vec.get_mut(0) {
                    *gc = gc_sum.iter().sum();
                }
            }
            Ok(())
        }

        fn calc_percentages(args: (u64, u64)) -> f64 {
            // Auxiliary function for calculating the percentage of nucleotides.
            let (content, total) = args;
            (((content as f64) / (total as f64)) * 100.0 * 100.0).round() / 100.0
        }

        pub fn percentage_of_nucleotides_and_output(&mut self) -> Result<(), Box<dyn Error>> {
            // Calculates the percentage of nucleotides and outputs the results.
            let border_line_snowflake: String = "*".repeat(25);
            let border_line_dash: String = "-".repeat(25);

            println!();
            println!("{}", border_line_snowflake);
            println!("Filename: {}", &self.filename);
            println!("{}", border_line_dash);
            println!("File header: {}", &self.file_header);
            println!("{}", border_line_dash);
            println!("Date of processing: {}", &self.formatted_datetime);
            println!("{}", border_line_dash);

            let mut total_number_of_characters: u64 = 0;
            let mut total_number_of_characters_perc: f64 = 0.0;
            
            if let Some(count_vec) = self.other_calc_nucleatides.get_mut(
                Self::OTHER_PARAMETERS_NUCLEATIDE_CONST[0]
            ) {
                if let Some(NucleotideCountType::UnInt(total_count)) = count_vec.get_mut(0) {
                    total_number_of_characters = *total_count;
                }
                if let Some(NucleotideCountType::Float(total_count_perc)) = count_vec.get_mut(1) {
                    *total_count_perc = 100.0;
                    total_number_of_characters_perc = *total_count_perc;
                }
                println!(
                    "Total number of characters {} - {}%", 
                    total_number_of_characters, 
                    total_number_of_characters_perc
                );
            }
            
            let mut atgcu: u64 = 0;
            let mut atgcu_perc_c: f64 = 0.0;
            if let Some(count_vec) = self.other_calc_nucleatides.get_mut(
                Self::OTHER_PARAMETERS_NUCLEATIDE_CONST[1]
            ) {
                if let Some(NucleotideCountType::UnInt(atgcu_count)) = count_vec.get_mut(0) {
                    atgcu = *atgcu_count;
                }
                if let Some(NucleotideCountType::Float(actgu_perc)) = count_vec.get_mut(1) {
                    *actgu_perc = Self::calc_percentages((atgcu, total_number_of_characters));
                    atgcu_perc_c = *actgu_perc;
                }
                println!("Number of ATGCU {} - {}%", atgcu, atgcu_perc_c);
            }

            let mut gc_content: u64 = 0;
            let mut gc_content_perc: f64 = 0.0;
            if let Some(count_vec) = self.other_calc_nucleatides.get_mut(
                Self::OTHER_PARAMETERS_NUCLEATIDE_CONST[2]
            ) {
                if let Some(NucleotideCountType::UnInt(gc_nucl_content)) = count_vec.get_mut(0) {
                    gc_content = *gc_nucl_content;
                }
                if let Some(NucleotideCountType::Float(gc_nucl_content_perc)) = count_vec.get_mut(1) {
                    *gc_nucl_content_perc = Self::calc_percentages((gc_content, total_number_of_characters)); 
                    gc_content_perc = *gc_nucl_content_perc;
                }
                println!(
                    "Number of GC content in characters {} - {}%",
                    gc_content, gc_content_perc
                );
            }

            let mut nucleatide_count: u64 = 0;
            let mut nucleatide_count_perc: f64 = 0.0;
            for nucleatide in &Self::NUCLEOTIDES_LIST {
                if let Some(count_vec) = self.nucleatides.get_mut(&nucleatide) {
                    if let Some(NucleotideCountType::UnInt(value)) = count_vec.get_mut(0) {
                        nucleatide_count = *value;
                    }
                    if let Some(NucleotideCountType::Float(value)) = count_vec.get_mut(1) {
                        *value = Self::calc_percentages((nucleatide_count, total_number_of_characters));
                        nucleatide_count_perc = *value;

                    }
                    if nucleatide == &Self::NUCLEOTIDES_LIST[5] {
                        println!(
                            "Other characters - {} - {}%", 
                            nucleatide_count, 
                            nucleatide_count_perc
                        );
                    } else {
                        println!(
                            "{} - {} - {}%", nucleatide, 
                            nucleatide_count, 
                            nucleatide_count_perc
                        );
                    }
                }
            }
            println!("{}", border_line_snowflake);
            println!();

            Ok(())
        }
    } 

}


/// Runs the nucleotide counting application.
///
/// # Arguments
///
/// * `config` - An instance of `NucleotideCounter` containing configuration parameters.
///
/// # Returns
///
/// * `Result` indicating success or an error message.
pub fn run(mut config: NucleotideCounter) -> Result<(), Box<dyn Error>> {
    config.read_nuacliotides()?;
    config.count_other_nucleatides()?;
    config.percentage_of_nucleotides_and_output()?;
    Ok(())
}

#[cfg(test)]
mod tests {
    // The test module.
    use super::*;
    use std::{
        error::Error, collections::HashMap,
    };
    use chrono::{DateTime, Datelike, Local, Timelike};
    use tempfile::{NamedTempFile, TempPath};
    use std::io::Write;

    #[test]
    fn run_nucleatides() -> Result<(), Box<dyn Error>> {
        // Checking the correct nucleotide count.
        let contents: &str = "\
>1 dna:chromosome chromosome:GRCh38:1:1:248956422:1 REF
NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNTAACCCTAACCCTAACCCTA
ACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTA
ACCCTAACCCTAACCCTAACCCTAACCCAACCCTAACCCTAACCCTAACCCTAACCCTAA
CCCTAACCCCTAACCCTAACCCTAACCCTAACCCTAACCTAACCCTAACCCTAACCCTAA
CCCTAACCCTAACCCTAACCCTAACCCTAACCCCTAACCCTAACCCTAAACCCTAAACCC
TAACCCTAACCCTAACCCTAACCCTAACCCCAACCCCAACCCCAACCCCAACCCCAACCC
CAACCCTAACCCCTAACCCTAACCCTAACCCTACCCTAACCCTAACCCTAACCCTAACCC
TAACCCTAACCCCTAACCCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACG";
    
    let mut file: NamedTempFile = NamedTempFile::new()?;
    writeln!(file, "{}", contents)?;

    let current_datetime: DateTime<Local> = Local::now();
        let formatted_datetime: String = format!(
            "{:02}:{:02} - {:02}.{:02}.{}",
            current_datetime.hour(),
            current_datetime.minute(),
            current_datetime.day(),
            current_datetime.month(),
            current_datetime.year(),
        );

    let mut config: NucleotideCounter = NucleotideCounter {
        file_path: file.path().to_string_lossy().into_owned(),
        filename: "test_filename".to_string(),
        file_header: "test_header".to_string(),
        formatted_datetime: formatted_datetime,
        nucleatides: HashMap::new(),
        other_calc_nucleatides: HashMap::new(),
    };

    for nucliotide in ['A', 'C', 'G', 'T', 'U', 'O'].iter() {
        config.nucleatides.insert(*nucliotide, vec![
            NucleotideCountType::UnInt(0), 
            NucleotideCountType::Float(0.0)
            ]
        );
    }

    let other_parameters_nucleatide_count: [&'static str; 3] = [
            "Total number of characters", "Number of ATGCU", 
            "Number of GC content in characters"
        ];

    for parameter in other_parameters_nucleatide_count {
        config.other_calc_nucleatides.insert(parameter.to_string(), vec![
            NucleotideCountType::UnInt(0), 
            NucleotideCountType::Float(0.0)
            ]
        );
    }

    config.read_nuacliotides()?;
    config.count_other_nucleatides()?;
    config.percentage_of_nucleotides_and_output()?;
    let test_file_header: &'static str = 
        "1 dna:chromosome chromosome:GRCh38:1:1:248956422:1 REF";

    assert_eq!(
        config.nucleatides.get(&'A'), 
        Some(&vec![NucleotideCountType::UnInt(147), 
        NucleotideCountType::Float(30.95)])
    );
    assert_eq!(
        config.nucleatides.get(&'C'), 
        Some(&vec![NucleotideCountType::UnInt(227), 
        NucleotideCountType::Float(47.79)])
    );
    assert_eq!(
        config.nucleatides.get(&'G'), 
        Some(&vec![NucleotideCountType::UnInt(1), 
        NucleotideCountType::Float(0.21)])
    );
    assert_eq!(
        config.nucleatides.get(&'T'), 
        Some(&vec![NucleotideCountType::UnInt(66), 
        NucleotideCountType::Float(13.89)])
    );
    assert_eq!(
        config.nucleatides.get(&'U'), 
        Some(&vec![NucleotideCountType::UnInt(0), 
        NucleotideCountType::Float(0.00)])
    );
    assert_eq!(
        config.nucleatides.get(&'O'), 
        Some(&vec![NucleotideCountType::UnInt(34), 
        NucleotideCountType::Float(7.16)])
    );

    assert_eq!(
        config.other_calc_nucleatides.get(
            "Total number of characters"), 
            Some(&vec![NucleotideCountType::UnInt(475), 
            NucleotideCountType::Float(100.0)])
        );
    assert_eq!(
        config.other_calc_nucleatides.get(
        "Number of ATGCU"), 
        Some(&vec![NucleotideCountType::UnInt(441), 
        NucleotideCountType::Float(92.84)])
    );
    assert_eq!(
        config.other_calc_nucleatides.get(
        "Number of GC content in characters"), 
        Some(&vec![NucleotideCountType::UnInt(228), 
        NucleotideCountType::Float(48.0)])
    );
    
    assert_eq!(config.file_header, test_file_header.to_string());

    let temp_path: TempPath = file.into_temp_path();
    temp_path.close()?;

    Ok(())
    }
}