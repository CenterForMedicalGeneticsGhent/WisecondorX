package utils

import (
	"fmt"

	"github.com/sbinet/npyio/npz"
)

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
