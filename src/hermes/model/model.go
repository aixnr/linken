package model

import (
	"log"
)

type Reads struct {
	// the input CSV (RunParam.InputCSV) needs to have these columns (no header)
	//   Sample
	//   R1
	//   R2 (for paired)
	Sample  string
	R1      string
	R2      string
	Trimmed []string
	Bwa     string
}

type RunParam struct {
	// This struct is for use with flag calls at main()
	InputCsv    *string
	ReadMode    *string
	ThreadCount *string
	SeqPlatform *string
	BwaIndex    *string
}

type Collect struct {
	// This struct is used mainly for interpolating text/template calls
	// R1 and R2 are trimmed reads from trimgalore
	// Collect struct is not shared between Process calls
	R1       string
	R2       string
	BwaIndex string
	Threads  string
	Bwa      string
	Sample   string
}

func (rp RunParam) Echo() {
	log.Printf("[INFO] Received input table '%s'.\n", *rp.InputCsv)
	log.Printf("[INFO] Configured with '%s' as the read mode (1 = single, 2 = paired).\n", *rp.ReadMode)
	log.Printf("[INFO] Sequencing platform: %s.\n", *rp.SeqPlatform)
	log.Printf("[INFO] Default executing threads n=%s.\n", *rp.ThreadCount)
}

func (read *Reads) SetTrimmed(trimmed []string) {
	read.Trimmed = trimmed
}
