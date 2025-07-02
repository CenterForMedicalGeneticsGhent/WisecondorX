package wisecondorx

import (
	"fmt"

	"github.com/sbinet/npyio/npz"
)

func IntSliceSum(intSlice []int32) int32 {
	// Calculate the sum of the slice
	var total int32 = 0
	for _, v := range intSlice {
		total += v
	}
	return total
}

func LoadSampleNpzFile(infile string) (samples map[string][]int32, binsize int32, err error) {
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
		var bins []int32
		err = in.Read(k, &bins)
		if err != nil {
			fmt.Println("Error: Unable to read data from file", infile)
			return
		}
		samples[k] = bins
	}
	return samples, binsize, nil
}

func LoadReferenceNpzFile(reference string) (samples map[string][]int32, binsize int32, err error) {
	// Load the reference npz file
	ref, err := npz.Open(reference)
	if err != nil {
		fmt.Println("Error: Unable to open file", reference)
		return nil, 0, err
	}
	defer ref.Close()

	// Get the binsize of the reference file
	binsize = 0
	err = ref.Read("binsize", &binsize)
	if err != nil {
		fmt.Println("Error: Unable to read binsize from file", reference)
		return nil, 0, err
	}

	// Load the data per chromosome
	for _, k := range ref.Keys() {
		var bins []int32
		err = ref.Read(k, &bins)
		if err != nil {
			fmt.Println("Error: Unable to read data from file", reference)
			return
		}
		samples[k] = bins
	}
	return samples, binsize, nil
}

func ScaleSample(sample map[string][]int32, from int, to int) (map[string][]int32, error) {
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
