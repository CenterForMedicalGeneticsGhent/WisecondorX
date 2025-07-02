package wisecondorx

import (
	"fmt"
)

func WcxNewRef(infiles []string, prefix string, nipt bool, yfrac float64, refsize int, binsize int) {
	// Create a new reference using healthy reference samples

	samples := make([]map[string][]int32, len(infiles))
	// Load the input files
	for _, infile := range infiles {
		infileData, infileBinsize, err := LoadSampleNpzFile(infile)
		if err != nil {
			fmt.Println("Error: Unable to load file: ", infile)
			return
		}

		// Scale the sample to the new bin size
		scaledSample, err := ScaleSample(infileData, int(infileBinsize), binsize)
		if err != nil {
			fmt.Println("Error: Unable to scale sample: ", err)
			return
		}
		samples = append(samples, scaledSample)
	}

	sexes, _ := trainSexDeterminationModel(samples, yfrac)

	// Do some sanity checks if a NIPT reference is requested
	if nipt == true {
		femaleCount := 0
		for _, sex := range sexes {
			if sex == "F" {
				femaleCount++
			}
		}
		if femaleCount < 5 {
			fmt.Println("A NIPT reference should have at least 5 female feti samples")
			return
		}
	}
	if nipt == false {
		for i, sample := range samples {
			samples[i] = correctSex(sample, sexes[i])
		}
	}
}

func trainSexDeterminationModel(samples []map[string][]int32, yfrac float64) (sexes []string, cutoff float64) {
	// Train the sex determination model
	sexes = make([]string, len(samples))
	// Create the training data
	//sexes := make([]string, len(samples))
	yfracs := make([]float64, len(samples))

	// Calculate the Y fraction for each sample
	for _, sample := range samples {
		sliceSumY := IntSliceSum(sample["chrY"])
		sliceSumAll := 0.0
		for _, bins := range sample {
			sliceSumAll += float64(IntSliceSum(bins))
		}
		yfrac := float64(sliceSumY) / sliceSumAll
		yfracs = append(yfracs, yfrac)
	}

	// Train the Guassian Mixture Model
	return sexes, cutoff
}

func correctSex(sample map[string][]int32, sex string) map[string][]int32 {
	//Levels gonosomal reads with the one at the autosomes.
	if sex == "M" {
		correctedChrX := make([]int32, len(sample["chrX"]))
		for _, v := range sample["chrX"] {
			correctedChrX = append(correctedChrX, v*2)
		}
		sample["chrX"] = correctedChrX
		correctedChrY := make([]int32, len(sample["chrY"]))
		for _, v := range sample["chrY"] {
			correctedChrY = append(correctedChrY, v*2)
		}
		sample["chrY"] = correctedChrY
	}
	return sample
}

func calculateMask(samples []map[string][]int32) (totalMask int, binsPerChromosome map[string][]int32) {
	// Calculate the mask for the samples
	totalMask = 0
	binsPerChromosome = make(map[string][]int32)

	// Get all chromosomes
	chromosomes := make([]string, 0)
	for chr, _ := range samples[0] {
		chromosomes = append(chromosomes, chr)
	}
	// Get max length of the bin slice among all chromosomes
	maxLen := 0
	for _, sample := range samples {
		for _, bins := range sample {
			if len(bins) > maxLen {
				maxLen = len(bins)
			}
		}
	}

	return totalMask, binsPerChromosome
}
