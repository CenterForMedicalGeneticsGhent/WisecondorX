package wisecondorx

import (
	"fmt"

	"github.com/sbinet/npyio/npz"
)

func WcxPredict(infile string, reference string, prefix string, minrefbins int, maskrepeats int, alpha float64, beta float64, zscore float64, blacklist string, sex string, seed int) {
	// Predict chromosomal aberrations

	// Load the input npz file
	inFileData, inFileBinsize, err := LoadSampleNpzFile(infile)
	if err != nil {
		fmt.Println("Error: Unable to load file: ", infile)
		return
	}

	// Load the reference npz file
	_, refFileBinsize, err := LoadReferenceNpzFile(reference)
	if err != nil {
		fmt.Println("Error: Unable to load file: ", reference)
		return
	}

	// Scale the sample to the reference bin size
	_, err = ScaleSample(inFileData, int(inFileBinsize), int(refFileBinsize))
	if err != nil {
		fmt.Println("Error: Unable to scale sample: ", err)
		return
	}

	// Load the reference npz file
	ref, err := npz.Open(reference)
	if err != nil {
		fmt.Println("Error: Unable to open file", reference)
		return
	}
	defer ref.Close()
}
