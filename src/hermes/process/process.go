package process

import (
	"fmt"
	"linken-hermes/model"
	"log"
	"os/exec"
	"strings"
)

func ProcessCollect(reads []model.Reads, rp model.RunParam) {
	// Create the required directories
	CreateDir()

	// Start looping, sample by sample
	for _, read := range reads {
		log.Printf("[INFO] ProcessCollect: Processing %s,\n", read.Sample)
		ProcessFastQc(&read, rp)
		ProcessTrim(&read, rp)
		ProcessMap(&read, rp)
		ProcessCalibrateAndCall(&read, rp)
		ProcessCoverage(&read, rp)
		ProcessExtractVariantToTable(&read, rp)
	}
}

func ProcessFastQc(read *model.Reads, rp model.RunParam) {
	if *rp.ReadMode == "1" {
		fmt.Println("Sorry, not implemented yet.")

	} else if *rp.ReadMode == "2" {
		cmd_fastqc := exec.Command("fastqc", "--threads", *rp.ThreadCount, read.R1, read.R2, "--outdir", "output_fastqc/")

		_, err := cmd_fastqc.Output()
		if err != nil {
			log.Fatalf("[ERR] FASTQC failed at %s, %s", read.Sample, err.Error())
		}
	}

	log.Println("[INFO] Completed ProcessFastQc.")
}

func ProcessTrim(read *model.Reads, rp model.RunParam) {
	if *rp.ReadMode == "1" {
		fmt.Println("Sorry, not implemented yet")

	} else if *rp.ReadMode == "2" {
		collect := model.Collect{
			Sample:  read.Sample,
			R1:      read.R1,
			R2:      read.R2,
			Threads: *rp.ThreadCount,
		}

		ShellCaller("trimgalore.sh", collect)
	}

	log.Println("[INFO] Completed ProcessTrim.")

	log.Println("[INFO] Collecting trimmed outputs...")
	trimmed := []string{}
	raw := []string{read.R1, read.R2}

	for i, r := range raw {
		str_noExt := strings.Split(r, ".fastq.gz")[0]
		str_noPrefix := strings.Split(str_noExt, "raw_reads/")[1]
		str_trimGalore := fmt.Sprintf("output_trimgalore/%s_val_%d.fq.gz", str_noPrefix, i+1)
		trimmed = append(trimmed, str_trimGalore)
	}

	read.SetTrimmed(trimmed)
}

func ProcessMap(read *model.Reads, rp model.RunParam) {
	if *rp.ReadMode == "1" {
		fmt.Println("Sorry, not implemented yet")
	} else if *rp.ReadMode == "2" {

		read.Bwa = fmt.Sprintf("%s.bwa.bam", read.Sample)

		collect := model.Collect{
			R1:       read.Trimmed[0],
			R2:       read.Trimmed[1],
			Threads:  *rp.ThreadCount,
			BwaIndex: *rp.BwaIndex,
			Sample:   read.Sample,
			Bwa:      read.Bwa,
		}

		ShellCaller("map.sh", collect)
		log.Printf("[INFO] Completed assembling %s.", read.Sample)
	}
}

func ProcessCalibrateAndCall(read *model.Reads, rp model.RunParam) {
	// This function is specific to illumina platform!
	// Check the bash script template, see templates/call.sh
	collect := model.Collect{
		Threads:  *rp.ThreadCount,
		BwaIndex: *rp.BwaIndex,
		Sample:   read.Sample,
		Bwa:      read.Bwa,
	}

	log.Printf("[INFO] Calibrating and calling variants for %s...", read.Sample)
	ShellCaller("call.sh", collect)
	log.Printf("[INFO] Completed calling %s.", read.Sample)
}

func ProcessCoverage(read *model.Reads, rp model.RunParam) {
	collect := model.Collect{
		Bwa:      read.Bwa,
		BwaIndex: *rp.BwaIndex,
		Sample:   read.Sample,
	}

	log.Printf("[INFO] Measuring coverage depth for %s...", read.Sample)
	ShellCaller("coverage.sh", collect)
	log.Printf("[INFO] Completed measuring coverage depth for %s.", read.Sample)
}

func ProcessExtractVariantToTable(read *model.Reads, rp model.RunParam) {
	collect := model.Collect{
		Sample: read.Sample,
	}

	log.Printf("[INFO] Collecting called variants for %s into tsv file...", read.Sample)
	ShellCaller("variant-tsv.sh", collect)
	log.Printf("[INFO] Completed collecting called variants for %s.", read.Sample)
	log.Printf("[INFO] Generated plots for %s, see 'output_plots'.", read.Sample)
}
