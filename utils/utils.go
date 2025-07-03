package utils

import (
	"fmt"

	"github.com/sbinet/npyio/npz"
	"gonum.org/v1/gonum/floats"
)

func LoadSampleNpzFile(infile string) (samples map[string][]float64, binsize int32, err error) {
	// Load the input npz file
	in, err := npz.Open(infile)
	if err != nil {
		fmt.Println("Error: Unable to open file", infile)
		return nil, 0, err
	}
	defer in.Close()

	// Get the binsize of the input file
	binsize = 0
	err = in.Read("binsize", &binsize)
	if err != nil {
		fmt.Println("Error: Unable to read binsize from file", infile)
		return nil, 0, err
	}

	// Load the data per chromosome
	for _, k := range in.Keys() {
		var bins []float64
		err = in.Read(k, &bins)
		if err != nil {
			fmt.Println("Error: Unable to read data from file", infile)
			return
		}
		samples[k] = bins
	}
	return samples, binsize, nil
}

func ScaleSample(sample map[string][]float64, from int, to int) (map[string][]float64, error) {
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

	sampleScaled := make(map[string][]float64)
	scaleFactor := to / from

	for chromosome, bins := range sample {
		scaledBinsLen := len(bins) / scaleFactor
		sampleScaled[chromosome] = make([]float64, scaledBinsLen)
		for i := 0; i < scaledBinsLen; i++ {
			sampleScaled[chromosome][i] = floats.Sum(bins[i*scaleFactor : i*scaleFactor+scaleFactor])
		}

	}

	return sampleScaled, nil
}
