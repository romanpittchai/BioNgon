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

        let other_parameters_nucleatide_count: [&str; 3] = [
            "Total number of characters", "Number of ATGCU", 
            "Number of GC content in characters"
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

    fn read_nuacliotides(&mut self) -> Result<(), Box<dyn Error>> {
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

    fn count_other_nucleatides(&mut self) -> Result<(), Box<dyn Error>> {
        let mut total_number_of_characters: u64 = 0;
        for nucleatide in &self.nucleatides {
            let count_vec: &&Vec<NucleotideCountType> = &nucleatide.1;
            if let Some(NucleotideCountType::UnInt(value)) = count_vec.get(0) {
                total_number_of_characters += value;
            }
        }

        if let Some(count_vec) = self.other_calc_nucleatides.get_mut(
            "Total number of characters"
        ) {
            if let Some(NucleotideCountType::UnInt(total_count)) = count_vec.get_mut(0) {
                *total_count = total_number_of_characters;
            } else {
                count_vec.push(NucleotideCountType::UnInt(total_number_of_characters));
            }
        }

        let mut other_nucleatides: u64 = 0;
        if let Some(count_vec) = self.nucleatides.get_mut(&'O') {
            if let Some(NucleotideCountType::UnInt(other_nucleatides_count)) = count_vec.get_mut(0) {
                other_nucleatides = *other_nucleatides_count;
            }
        }

        if let Some(count_vec) = self.other_calc_nucleatides.get_mut(
            "Number of ATGCU"
        ) {
            if let Some(NucleotideCountType::UnInt(atgcu)) = count_vec.get_mut(0) {
                *atgcu = total_number_of_characters - other_nucleatides;
            }
        }
        
        let nucleatides_list: [char; 2] = ['G', 'C'];
        for nucleotide in &nucleatides_list {
            if let Some(count_vec) = self.other_calc_nucleatides.get_mut(
                "Number of GC content in characters"
            ) {
                if let Some(NucleotideCountType::UnInt(gc_content)) = count_vec.get_mut(0) {
                    if let Some(NucleotideCountType::UnInt(value)) = self.nucleatides.get(
                        nucleotide
                    ).and_then(|vec| vec.get(0)) {
                        *gc_content += value;
                    }
                }
            }
        }
        Ok(())
    }

    fn percentage_of_nucleotides(&mut self) -> Result<(), Box<dyn Error>> {
        let mut total_number_of_characters: u64 = 0;
        if let Some(count_vec) = self.other_calc_nucleatides.get_mut(
            "Total number of characters"
        ) {
            if let Some(NucleotideCountType::UnInt(total_count)) = count_vec.get_mut(0) {
                total_number_of_characters = *total_count;
            }
            if let Some(NucleotideCountType::Float(total_count_perc)) = count_vec.get_mut(1) {
                *total_count_perc = 100.0;
            }
        }

        let mut gc_content: u64 = 0;
        if let Some(count_vec) = self.other_calc_nucleatides.get_mut(
            "Number of GC content in characters"
        ) {
            if let Some(NucleotideCountType::UnInt(gc_nucl_content)) = count_vec.get_mut(0) {
                gc_content = *gc_nucl_content;
            }
            if let Some(NucleotideCountType::Float(gc_nucl_content_perc)) = count_vec.get_mut(1) {
                *gc_nucl_content_perc = ((
                    (gc_content as f64) / (total_number_of_characters as f64)
                ) * 100.0 * 100.0).round() / 100.0;
            }
        }
        
        let mut atgcu: u64 = 0;
        if let Some(count_vec) = self.other_calc_nucleatides.get_mut(
            "Number of ATGCU"
        ) {
            if let Some(NucleotideCountType::UnInt(atgcu_count)) = count_vec.get_mut(0) {
                atgcu = *atgcu_count;
            }
            if let Some(NucleotideCountType::Float(actgu_perc)) = count_vec.get_mut(1) {
                *actgu_perc = ((
                    (atgcu as f64) / (total_number_of_characters as f64)
                ) * 100.0 * 100.0).round() / 100.0;
            }
        }

        let mut other_nucl: u64 = 0;
        if let Some(count_vec) = self.other_calc_nucleatides.get_mut(
            "Other characters"
        ) {
            if let Some(NucleotideCountType::UnInt(other_nucl_count)) = count_vec.get_mut(0) {
                other_nucl = *other_nucl_count;
            }
            if let Some(NucleotideCountType::Float(other_nucl_perc)) = count_vec.get_mut(1) {
                *other_nucl_perc = ((
                    (other_nucl as f64) / (total_number_of_characters as f64)
                ) * 100.0 * 100.0).round() / 100.0;
            }
        }
        let nucleatides_list: [char; 6] = ['A', 'C', 'G', 'T', 'U', 'O'];
        let mut nucleatide_count: u64 = 0;
        for nucleatide in nucleatides_list {
            if let Some(count_vec) = self.nucleatides.get_mut(&nucleatide) {
                if let Some(NucleotideCountType::UnInt(value)) = count_vec.get_mut(0) {
                    nucleatide_count = *value;
                }
                if let Some(NucleotideCountType::Float(value)) = count_vec.get_mut(1) {
                    *value = ((
                        (nucleatide_count as f64) / (total_number_of_characters as f64)
                    ) * 100.0 * 100.0).round() / 100.0;
                }
            }
        }

        Ok(())
    }
} 


pub fn run(mut config: NucleotideCounter) -> Result<(), Box<dyn Error>> {
    config.read_nuacliotides()?;
    config.count_other_nucleatides()?;
    config.percentage_of_nucleotides()?;
    dbg!(config.file_path);
    dbg!(config.filename);
    dbg!(config.file_header);
    dbg!(config.other_calc_nucleatides);
    dbg!(config.nucleatides);
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
    config.count_other_nucleatides()?;
    config.percentage_of_nucleotides()?;
    let test_file_header: &str = "1 dna:chromosome chromosome:GRCh38:1:1:248956422:1 REF";

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