package wisecondorx

import (
	"fmt"

	"github.com/sbinet/npyio/npz"
)

func WcxNewRef(infiles []string, prefix string, nipt bool, yfrac float64, refsize int, binsize int) {
	// Create a new reference using healthy reference samples

	samples := make([]map[string][]int32, len(infiles))
	// Load the input files
	for _, infile := range infiles {
		// Load the input file
		f, err := npz.Open(infile)
		if err != nil {
			fmt.Println("Error: Unable to open file", infile)
			return
		}
		defer f.Close()

		// Get the binsize of the input file
		var sourceBinsize int32 = 0
		err = f.Read("binsize", &binsize)
		if err != nil {
			fmt.Println("Error: Unable to read binsize from file", infile)
			return
		}

		// Load the data per chromosome
		infileData := make(map[string][]int32)
		for _, k := range f.Keys() {
			var bins []int32
			err = f.Read(k, &bins)
			if err != nil {
				fmt.Println("Error: Unable to read data from file", infile)
				return
			}
			infileData[k] = bins
		}

		// Scale the sample to the new bin size
		scaledSample, err := scaleSample(infileData, int(sourceBinsize), binsize)
		if err != nil {
			fmt.Println("Error: Unable to scale sample: ", err)
			return
		}
		samples = append(samples, scaledSample)
	}

	sexes, cutoff := trainSexDeterminationModel(samples, yfrac)

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

	totalMask, binsPerChromosome := calculateMask(samples)

}

func scaleSample(sample map[string][]int32, from int, to int) (map[string][]int32, error) {
	// Scale the sample to a new bin size

	// Check if the new bin size is equal to the old bin size
	if to == from {
		return sample, nil
	}

	// Check if the new bin size is larger than the old bin size
	if to < from {
		return nil, fmt.Errorf("new bin size must be larger than the old bin size")
	}

	// Check if the new bin size is a multiple of the old bin size
	if to%from != 0 {
		fmt.Println("Error: New bin size must be a multiple of the old bin size")
		return nil, fmt.Errorf("new bin size must be a multiple of the old bin size")
	}

	sampleScaled := make(map[string][]int32)
	scaleFactor := to / from

	for chromosome, bins := range sample {
		scaledBinsLen := len(bins) / scaleFactor
		sampleScaled[chromosome] = make([]int32, scaledBinsLen)
		for i := 0; i < scaledBinsLen; i++ {
			sampleScaled[chromosome][i] = IntSliceSum(bins[i*scaleFactor : i*scaleFactor+scaleFactor])
		}

	}

	return sampleScaled, nil
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
