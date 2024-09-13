package wisecondorx

func IntSliceSum(intSlice []int32) int32 {
	// Calculate the sum of the slice
	var total int32 = 0
	for _, v := range intSlice {
		total += v
	}
	return total
}
