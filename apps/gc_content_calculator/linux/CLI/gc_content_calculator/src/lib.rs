use std::collections::HashMap;
use std::error::Error;
use std::path::Path;
use std::fs::File;
use std::io::{BufRead, BufReader};

#[derive(Debug, Clone, PartialEq)]
pub enum NucleotideCountType {
    UnInt(u64),
    Float(f64),
}

pub struct NucleotideCounter {
    pub file_path: String,
    pub filename: String,
    pub file_header: String,
    pub nucleatides: HashMap<char, Vec<NucleotideCountType>>,
    pub other_calc_nucleatides: HashMap<String, Vec<NucleotideCountType>>,
}

impl NucleotideCounter {
    pub fn build(args: &[String]) -> Result<NucleotideCounter, &'static str> {
        if args.len() < 2 {
            return Err("Not enough arguments!");
        }

        let file_path: &str = &args[1];
        let path: &Path = Path::new(file_path);
        let filename: String = if let Some(file_name) = path.file_name() {
            file_name.to_str().unwrap_or_default().to_string()
        } else {
            String::new()
        };
        let file_path: String = file_path.to_string();
        
        let file_header: String = String::from("Empty");

        let nucleatides_list: [char; 6] = ['A', 'C', 'G', 'T', 'U', 'O'];
        let mut nucleatides: HashMap<char, Vec<NucleotideCountType>> = HashMap::new();

        for nucleatide in nucleatides_list {
            nucleatides.insert(nucleatide, vec![
                NucleotideCountType::UnInt(0), 
                NucleotideCountType::Float(0.0)
                ]
            );
        }

        let other_parameters_nucleatide_count: [&str; 4] = [
            "Total number of characters", "Number of ATGCU", 
            "Other characters", "Number of GC content in characters"
        ];
        let mut other_calc_nucleatides: HashMap<String, Vec<NucleotideCountType>> = HashMap::new();

        for parameter in other_parameters_nucleatide_count {
            other_calc_nucleatides.insert(parameter.to_string(), vec![
                NucleotideCountType::UnInt(0), 
                NucleotideCountType::Float(0.0)
                ]
            );
        }

        Ok(NucleotideCounter { file_path, filename, file_header, nucleatides, other_calc_nucleatides })
    }

    pub fn read_nuacliotides(&mut self) -> Result<(), Box<dyn Error>> {
        let file: File = File::open(&self.file_path)?;
        let reader: BufReader<File> = BufReader::with_capacity(8 * 1024, file);
    
        for line in reader.lines() {
            let sequence: String = line?;
            let header_symbol: char = '>';

            if sequence.starts_with(header_symbol) {
                self.file_header = sequence.chars().skip(1).collect();
                continue;
            }

            for c in sequence.chars() {
                if let Some(count_vec) = self.nucleatides.get_mut(&c) {
                    if let Some(NucleotideCountType::UnInt(value)) = count_vec.get_mut(0) {
                        *value += 1;
                    } else {
                        count_vec.push(NucleotideCountType::UnInt(1));
                    }
                } else {
                    let char_o: char = 'O';
                    if let Some(count_vec_o) = self.nucleatides.get_mut(&char_o) {
                        if let Some(NucleotideCountType::UnInt(value)) = count_vec_o.get_mut(0) {
                            *value += 1;
                        } else {
                            count_vec_o.push(NucleotideCountType::UnInt(1));
                        }
                    }
                }
            }
        }
        Ok(())
    }

    pub fn count_other_nucleatides(&mut self) -> Result<(), Box<dyn Error>> {
        Ok(())
    }

    pub fn percentage_of_nucleotides(&mut self) -> Result<(), Box<dyn Error>> {
        Ok(())
    }
} 


pub fn run(mut config: NucleotideCounter) -> Result<(), Box<dyn Error>> {
    config.read_nuacliotides()?;
    dbg!(config.file_path);
    dbg!(config.filename);
    dbg!(config.file_header);
    dbg!(config.nucleatides);
    dbg!(config.other_calc_nucleatides);
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use tempfile::{NamedTempFile, TempPath};
    use std::io::Write;

    #[test]
    fn run_nucleatides() -> Result<(), Box<dyn Error>>{
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

    let mut config: NucleotideCounter = NucleotideCounter {
        file_path: file.path().to_string_lossy().into_owned(),
        filename: "test_filename".to_string(),
        file_header: "test_header".to_string(),
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

    let other_parameters_nucleatide_count: [&str; 4] = [
            "Total number of characters", "Number of ATGCU", 
            "Other characters", "Number of GC content in characters"
        ];

    for parameter in other_parameters_nucleatide_count {
        config.other_calc_nucleatides.insert(parameter.to_string(), vec![
            NucleotideCountType::UnInt(0), 
            NucleotideCountType::Float(0.0)
            ]
        );
    }

    config.read_nuacliotides()?;
    let test_file_header: &str = "1 dna:chromosome chromosome:GRCh38:1:1:248956422:1 REF";

    assert_eq!(config.nucleatides.get(&'A').map(|v| v[0].clone()), Some(NucleotideCountType::UnInt(147)));
    assert_eq!(config.nucleatides.get(&'C').map(|v| v[0].clone()), Some(NucleotideCountType::UnInt(227)));
    assert_eq!(config.nucleatides.get(&'G').map(|v| v[0].clone()), Some(NucleotideCountType::UnInt(1)));
    assert_eq!(config.nucleatides.get(&'T').map(|v| v[0].clone()), Some(NucleotideCountType::UnInt(66)));
    assert_eq!(config.nucleatides.get(&'U').map(|v| v[0].clone()), Some(NucleotideCountType::UnInt(0)));
    assert_eq!(config.nucleatides.get(&'O').map(|v| v[0].clone()), Some(NucleotideCountType::UnInt(34)));

    assert_eq!(config.other_calc_nucleatides.get("Total number of characters"), Some(&vec![NucleotideCountType::UnInt(0), NucleotideCountType::Float(0.0)]));
    assert_eq!(config.other_calc_nucleatides.get("Number of ATGCU"), Some(&vec![NucleotideCountType::UnInt(0), NucleotideCountType::Float(0.0)]));
    assert_eq!(config.other_calc_nucleatides.get("Other characters"), Some(&vec![NucleotideCountType::UnInt(0), NucleotideCountType::Float(0.0)]));
    assert_eq!(config.other_calc_nucleatides.get("Number of GC content in characters"), Some(&vec![NucleotideCountType::UnInt(0), NucleotideCountType::Float(0.0)]));
    
    assert_eq!(config.file_header, test_file_header.to_string());

    let temp_path: TempPath = file.into_temp_path();
    temp_path.close()?;

    Ok(())
    }
}