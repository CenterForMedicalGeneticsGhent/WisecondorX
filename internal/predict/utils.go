package predict

import (
	"fmt"

	"github.com/sbinet/npyio/npz"
)

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
