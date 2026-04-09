#![allow(unsafe_op_in_unsafe_fn)]
use pyo3::prelude::*;
use pyo3::types::PyDict;
use numpy::{IntoPyArray, ndarray::Array1};
use noodles::bam;
use noodles::sam;
use std::fs::File;
use std::path::Path;
use anyhow::Result;
use std::collections::HashMap;
use rayon::prelude::*;

#[pyclass(from_py_object)]
#[derive(Clone)]
pub struct BamQualityInfo {
    #[pyo3(get)]
    mapped: usize,
    #[pyo3(get)]
    unmapped: usize,
    #[pyo3(get)]
    no_coordinate: usize,
    #[pyo3(get)]
    filter_rmdup: usize,
    #[pyo3(get)]
    filter_mapq: usize,
    #[pyo3(get)]
    pre_retro: usize,
    #[pyo3(get)]
    post_retro: usize,
    #[pyo3(get)]
    pair_fail: usize,
}

impl BamQualityInfo {
    fn new() -> Self {
        BamQualityInfo {
            mapped: 0,
            unmapped: 0,
            no_coordinate: 0,
            filter_rmdup: 0,
            filter_mapq: 0,
            pre_retro: 0,
            post_retro: 0,
            pair_fail: 0,
        }
    }
}

#[pyfunction]
#[pyo3(signature = (infile, reference, binsize, rmdup, threads))]
pub fn wcx_convert_core(
    py: Python,
    infile: &str,
    reference: Option<&str>,
    binsize: usize,
    rmdup: bool,
    threads: usize,
) -> PyResult<(Py<PyDict>, BamQualityInfo)> {
    
    let path = Path::new(infile);
    if !path.exists() {
        return Err(pyo3::exceptions::PyFileNotFoundError::new_err(format!("File {} not found", infile)));
    }

    let pbuilder = rayon::ThreadPoolBuilder::new().num_threads(threads);
    let pool = pbuilder.build().map_err(|e| {
        pyo3::exceptions::PyRuntimeError::new_err(format!("Rayon thread pool init error: {}", e))
    })?;

    let output = pool.install(|| {
        process_file(path, reference, binsize, rmdup)
    });

    match output {
        Ok((counts, qual)) => {
            let dict = PyDict::new(py);
            for (chr, array) in counts {
                let py_array = array.into_pyarray(py);
                dict.set_item(chr, py_array)?;
            }
            Ok((dict.into(), qual))
        },
        Err(e) => Err(pyo3::exceptions::PyRuntimeError::new_err(format!("Error: {}", e)))
    }
}

fn process_record<R: noodles::sam::alignment::Record>(
    record: &R,
    rmdup: bool,
    binsize: usize,
    thread_counts: &mut [i32],
    larp: &mut isize,
    larp2: &mut isize,
    bqi: &mut BamQualityInfo,
) -> Result<()> {
    let flags = record.flags()?;
    if flags.is_unmapped() { bqi.unmapped += 1; } else { bqi.mapped += 1; }
    let pos_opt = record.alignment_start().transpose()?.map(usize::from);
    if pos_opt.is_none() { bqi.no_coordinate += 1; }
    
    let is_paired = flags.is_segmented();
    let is_proper_pair = flags.is_properly_segmented();
    let mapq = record.mapping_quality().transpose()?.map(u8::from).unwrap_or(255);
    let pos = pos_opt.unwrap_or(1) as isize - 1; 
    let mate_pos = record.mate_alignment_start().transpose()?.map(usize::from).unwrap_or(1) as isize - 1;

    if is_paired {
        if !is_proper_pair { bqi.pair_fail += 1; return Ok(()); }
        if rmdup && *larp == pos && *larp2 == mate_pos { bqi.filter_rmdup += 1; }
        else if mapq >= 1 {
            let loc = (pos as usize) / binsize;
            if loc < thread_counts.len() { thread_counts[loc] += 1; }
        } else { bqi.filter_mapq += 1; }

        *larp2 = mate_pos; bqi.pre_retro += 1; *larp = pos;
    } else {
        if rmdup && *larp == pos { bqi.filter_rmdup += 1; }
        else if mapq >= 1 {
            let loc = (pos as usize) / binsize;
            if loc < thread_counts.len() { thread_counts[loc] += 1; }
        } else { bqi.filter_mapq += 1; }
        bqi.pre_retro += 1; *larp = pos;
    }
    Ok(())
}

fn process_file(
    path: &Path,
    _reference: Option<&str>,
    binsize: usize,
    rmdup: bool,
) -> Result<(HashMap<String, Array1<i32>>, BamQualityInfo)> {
    
    let ext = path.extension().unwrap_or_default().to_string_lossy();
    
    let define_refs = |header: &sam::Header| -> HashMap<usize, (String, usize)> {
        let mut refs = HashMap::new();
        for (i, (name, reference)) in header.reference_sequences().iter().enumerate() {
            let mut chr_name = String::from_utf8_lossy(name.as_ref()).to_string();
            if chr_name.to_lowercase().starts_with("chr") {
                chr_name = chr_name[3..].to_string();
            }
            if chr_name == "X" { chr_name = "23".to_string(); }
            if chr_name == "Y" { chr_name = "24".to_string(); }
            
            if let Ok(num) = chr_name.parse::<usize>() 
                && (1..=24).contains(&num) {
                refs.insert(i, (chr_name, reference.length().get()));
            }
        }
        refs
    };

    let valid_refs_map = if ext == "bam" {
        let mut reader = bam::io::Reader::new(File::open(path)?);
        let header = reader.read_header()?;
        define_refs(&header)
    } else if ext == "cram" {
        let mut reader = noodles::cram::io::Reader::new(File::open(path)?);
        reader.read_file_definition()?;
        let header = reader.read_file_header()?;
        define_refs(&header)
    } else {
        anyhow::bail!("Unsupported extension for parallel evaluation: {}", ext);
    };

    let regions_list: Vec<_> = valid_refs_map.iter().map(|(id, (name, len))| {
        (*id, name.clone(), *len)
    }).collect();

    let results: Result<Vec<_>> = regions_list.into_par_iter().map(|(_id, chr_name, len)| {
        let mut bqi = BamQualityInfo::new();
        
        let num_bins = len / binsize + 1;
        let mut thread_counts = vec![0; num_bins];
        
        let mut larp: isize = -1;
        let mut larp2: isize = -1;
        
        let region = match chr_name.parse::<noodles::core::Region>() {
            Ok(r) => r,
            Err(_) => return Ok((chr_name, thread_counts, bqi))
        };

        if ext == "bam" {
            if let Ok(mut reader) = bam::io::indexed_reader::Builder::default().build_from_path(path)
                && let Ok(header) = reader.read_header()
                && let Ok(query_iter) = reader.query(&header, &region) {
                    for record in query_iter.flatten() {
                        let _ = process_record(
                            &record, rmdup, binsize, &mut thread_counts, &mut larp, &mut larp2,
                            &mut bqi
                        );
                    }
                }
        } else if ext == "cram"
            && let Ok(mut reader) = noodles::cram::io::indexed_reader::Builder::default().build_from_path(path) {
                let _ = reader.read_file_definition();
                if let Ok(header) = reader.read_file_header()
                    && let Ok(query_iter) = reader.query(&header, &region) {
                        for record in query_iter.flatten() {
                            let _ = process_record(
                                &record, rmdup, binsize, &mut thread_counts, &mut larp, &mut larp2,
                                &mut bqi
                            );
                        }
                    }
                }

        for &count in &thread_counts {
            bqi.post_retro += count as usize;
        }

        Ok((chr_name.clone(), thread_counts, bqi))
    }).collect();

    let mut final_map = HashMap::new();
    let mut final_bqi = BamQualityInfo::new();
    
    for (chr, counts, qual) in results? {
        final_map.insert(chr, Array1::from_vec(counts));
        final_bqi.mapped += qual.mapped;
        final_bqi.unmapped += qual.unmapped;
        final_bqi.no_coordinate += qual.no_coordinate;
        final_bqi.filter_rmdup += qual.filter_rmdup;
        final_bqi.filter_mapq += qual.filter_mapq;
        final_bqi.pre_retro += qual.pre_retro;
        final_bqi.post_retro += qual.post_retro;
        final_bqi.pair_fail += qual.pair_fail;
    }

    Ok((final_map, final_bqi))
}

#[pymodule]
mod wisecondorx_rs {
    use super::*;

    #[pymodule_init]
    fn init(m: &Bound<'_, PyModule>) -> PyResult<()> {
        m.add_class::<BamQualityInfo>()?;
        m.add_function(wrap_pyfunction!(wcx_convert_core, m)?)?;
        Ok(())
    }
}