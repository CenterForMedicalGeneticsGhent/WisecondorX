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

func trainSexDeterminationModel(samples []map[string][]int32, yfrac float64) {
	// Train the sex determination model

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

}
