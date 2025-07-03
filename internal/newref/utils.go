package newref

import (
	"fmt"

	"gonum.org/v1/gonum/floats"
)

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
