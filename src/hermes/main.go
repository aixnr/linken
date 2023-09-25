package main

import (
	"encoding/csv"
	"flag"
	"linken-hermes/model"
	"linken-hermes/process"
	"log"
	"os"
)

func main() {
	inputCsv := flag.String("i", "paired.csv", "Input CSV file listing all raw fastq.gz read files.")
	readMode := flag.String("r", "2", "1 for single, 2 for paired sequencing.")
	seqPlatform := flag.String("p", "illumina", "Accepts either 'illumina' or 'nanopore'.")
	threadCount := flag.String("n", "4", "Number of threads to pass to executors.")
	bwaIndex := flag.String("ref", "./index/ref", "Reference index")

	flag.Parse()

	rp := model.RunParam{
		InputCsv:    inputCsv,
		ReadMode:    readMode,
		SeqPlatform: seqPlatform,
		ThreadCount: threadCount,
		BwaIndex:    bwaIndex,
	}

	rp.Echo()

	file, err := os.Open(*rp.InputCsv)
	if err != nil {
		log.Fatalf("Cannot open %s", *rp.InputCsv)
	}

	defer file.Close()

	pairedReads := CsvReader(file, rp)
	process.ProcessCollect(pairedReads, rp)
}

func CsvReader(file *os.File, rp model.RunParam) []model.Reads {
	reader := csv.NewReader(file)
	PairedReads := []model.Reads{}

	records, err := reader.ReadAll()
	if err != nil {
		log.Fatalf("Error reading records in %s.", *rp.InputCsv)
	}

	if *rp.ReadMode == "1" {
		// Use PairedReads, but leave empty for the second field
		for _, eachRecord := range records {
			PairedReads = append(PairedReads, model.Reads{
				Sample: eachRecord[0],
				R1:     eachRecord[1],
				R2:     ""})
		}
	}

	if *rp.ReadMode == "2" {
		// Use PaireadReads, now uses both fields
		for _, eachRecord := range records {
			PairedReads = append(PairedReads, model.Reads{
				Sample: eachRecord[0],
				R1:     eachRecord[1],
				R2:     eachRecord[2]})
		}
	}

	return PairedReads
}
