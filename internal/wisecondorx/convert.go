package wisecondorx

import (
	"fmt"
	"os"
	"regexp"

	"github.com/biogo/hts/bam"
	"github.com/biogo/hts/sam"
	"github.com/sbinet/npyio/npz"
)

func WcxConvert(infile string, prefix string, binsize int, RemoveDup bool, reference string) {
	// Convert and filter aligned reads to .wcx format

	// Read sam/bam/cram file
	// Open the file
	f, err := os.Open(infile)
	if err != nil {
		fmt.Printf("Error: Unable to open file %s", infile)
		return
	}
	defer f.Close()

	// if infile ends with .bam
	b, err := bam.NewReader(f, 0)
	if err != nil {
		fmt.Printf("Error: Unable to read bam file %s", infile)
		return
	}
	defer b.Close()

	// Generate a map of chromosome names and the length of the chromosome
	// and a map of the number of bins per chromosome

	canonicalChromosomeRegex := regexp.MustCompile(`^chr[0-9XY]+$`)

	binsPerChromosome := make(map[string][]int32)
	for _, ref := range b.Header().Refs() {
		// skip the chromosome if not autosomes or gonosomes
		if canonicalChromosomeRegex.Match([]byte(ref.Name())) == false {
			continue
		}
		binsPerChromosome[ref.Name()] = make([]int32, ref.Len()/binsize+1)
	}

	currentChromosome := ""
	larp := -1
	larp2 := -1

	// Loop over the reads
	for {
		samRecord, err := b.Read()
		if err != nil {
			break
		}

		// Skip the read if not on autosomes or gonosomes
		if canonicalChromosomeRegex.Match([]byte(samRecord.Ref.Name())) == false {
			continue
		}

		if currentChromosome != samRecord.Ref.Name() {
			fmt.Println("Processing chromosome: ", samRecord.Ref.Name())
			currentChromosome = samRecord.Ref.Name()
		}

		// Skip the read if duplicate && RemoveDup is true
		if samRecord.Flags&sam.Duplicate != 0 && RemoveDup == true {
			continue
		}

		// read is unpaired
		if samRecord.Flags&sam.Paired == 0 {
			if larp != samRecord.Pos {
				if samRecord.MapQ >= 1 {
					binsPerChromosome[samRecord.Ref.Name()][samRecord.Pos/binsize]++
				}
			}
			larp = samRecord.Pos
			// go to next read
			continue
		}

		// read is paired
		if samRecord.Flags&sam.Paired != 0 {
			// read is properly paired
			if samRecord.Flags&sam.ProperPair != 0 {
				if larp != samRecord.Pos && larp2 != samRecord.MatePos {
					if samRecord.MapQ >= 1 {
						binsPerChromosome[samRecord.Ref.Name()][samRecord.Pos/binsize]++
					}
				}
			}
			larp = samRecord.Pos
			larp2 = samRecord.MatePos
		}
	}

	// Write the output file
	npz, err := npz.Create(prefix + ".npz")
	if err != nil {
		fmt.Printf("Error: Unable to create file %s.npz: %s", prefix, err)
	}
	defer npz.Close()

	// Write bins per chromosome
	for chromosome, bins := range binsPerChromosome {
		err = npz.Write(chromosome, bins)
		if err != nil {
			fmt.Printf("Error: Unable to write to file %s.npz: %s", prefix, err)
		}
	}
	// Write binsize
	err = npz.Write("binsize", int32(binsize))
	if err != nil {
		fmt.Printf("Error: Unable to write to file %s.npz: %s", prefix, err)
	}

	err = npz.Close()
	if err != nil {
		fmt.Printf("Error: Unable to close file %s.npz: %s", prefix, err)
	}
}
