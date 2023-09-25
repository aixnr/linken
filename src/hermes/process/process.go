package process

import (
	"embed"
	"fmt"
	"linken-hermes/model"
	"log"
	"os"
	"os/exec"
	"strings"
	"text/template"
)

//go:embed templates/*
var scriptTemplate embed.FS

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
		ProcessResidueIdentifier(&read, rp)
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
		cmd_trimgalore := exec.Command("trimgalore", "--paired", read.R1, read.R2, "--cores", *rp.ThreadCount, "--output_dir", "output_trimgalore/")

		_, err := cmd_trimgalore.Output()
		if err != nil {
			log.Fatalf("[ERR] TRIMGALORE failed at %s, %s", read.Sample, err.Error())
		}
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

		path_script_map := "script_map.sh"
		script_map, err := os.Create(path_script_map)
		if err != nil {
			log.Fatalf("[ERR] Error creating %s at %s, %s", path_script_map, read.Sample, err.Error())
		}

		defer script_map.Close()

		t := template.Must(template.New("map.sh").ParseFS(scriptTemplate, "templates/map.sh"))
		if err := t.Execute(script_map, collect); err != nil {
			log.Fatalf("[ERR] Cannot write %s for %s, %s", path_script_map, read.Sample, err.Error())
		}

		log.Printf("[INFO] Assembling %s...", read.Sample)
		run_map := exec.Command("bash", path_script_map)
		_, err = run_map.Output()
		if err != nil {
			log.Fatalf("[ERR] Error assembling at %s, %s", read.Sample, err.Error())
		}

		// remove the script
		log.Printf("[INFO] Completed assembling %s.", read.Sample)
		os.Remove(path_script_map)
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

	path_script_call := "script_call.sh"
	script_call, err := os.Create(path_script_call)
	if err != nil {
		log.Fatalf("[ERR] Error creating %s at %s, %s", path_script_call, read.Sample, err.Error())
	}

	defer script_call.Close()

	t := template.Must(template.New("call.sh").ParseFS(scriptTemplate, "templates/call.sh"))
	if err := t.Execute(script_call, collect); err != nil {
		panic(err)
	}

	log.Printf("[INFO] Calibrating and calling variants for %s...", read.Sample)
	run_call := exec.Command("bash", path_script_call)
	_, err = run_call.Output()
	if err != nil {
		log.Fatalf("[ERR] Error calling at %s, %s", read.Sample, err.Error())
	}

	// remove the script
	log.Printf("[INFO] Completed calling %s.", read.Sample)
	os.Remove(path_script_call)
}

func ProcessCoverage(read *model.Reads, rp model.RunParam) {
	collect := model.Collect{
		Bwa:      read.Bwa,
		BwaIndex: *rp.BwaIndex,
		Sample:   read.Sample,
	}

	path_script_coverage := "coverage.sh"
	script_coverage, err := os.Create(path_script_coverage)
	if err != nil {
		log.Fatalf("[ERR] Error creating %s at %s, %s", path_script_coverage, read.Sample, err.Error())
	}

	defer script_coverage.Close()

	t := template.Must(template.New("coverage.sh").ParseFS(scriptTemplate, "templates/coverage.sh"))
	if err := t.Execute(script_coverage, collect); err != nil {
		panic(err)
	}

	log.Printf("[INFO] Measuring coverage depth for %s...", read.Sample)
	run_coverage := exec.Command("bash", path_script_coverage)
	_, err = run_coverage.Output()
	if err != nil {
		log.Fatalf("[ERR] Error measuring coverage at %s, %s.", read.Sample, err.Error())
	}

	// remove the script
	log.Printf("[INFO] Completed measuring coverage depth for %s.", read.Sample)
	os.Remove(path_script_coverage)
}

func ProcessExtractVariantToTable(read *model.Reads, rp model.RunParam) {
	collect := model.Collect{
		Sample: read.Sample,
	}

	path_script_vcftsv := "vcf-tsv.sh"
	script_vcftsv, err := os.Create(path_script_vcftsv)
	if err != nil {
		log.Fatalf("[ERR] Error creating %s at %s, %s", path_script_vcftsv, read.Sample, err.Error())
	}

	defer script_vcftsv.Close()

	t := template.Must(template.New("variant-tsv.sh").ParseFS(scriptTemplate, "templates/variant-tsv.sh"))
	if err := t.Execute(script_vcftsv, collect); err != nil {
		panic(err)
	}

	log.Printf("[INFO] Collecting called variants for %s into tsv file...", read.Sample)
	run_vcftsv := exec.Command("bash", path_script_vcftsv)
	_, err = run_vcftsv.Output()
	if err != nil {
		log.Fatalf("[ERR] Error collecting called variants at %s, %s.", read.Sample, err.Error())
	}

	// remove the script
	log.Printf("[INFO] Completed collecting called variants for %s.", read.Sample)
	log.Printf("[INFO] Generated plots for %s, see 'output_plots'.", read.Sample)
	os.Remove(path_script_vcftsv)
}

func ProcessResidueIdentifier(read *model.Reads, rp model.RunParam) {
	collect := model.Collect{
		Sample:   read.Sample,
		BwaIndex: *rp.BwaIndex,
	}

	path_script_residue := "residue.sh"
	script_residue, err := os.Create(path_script_residue)
	if err != nil {
		log.Fatalf("[ERR] Error creating %s at %s, %s", path_script_residue, read.Sample, err.Error())
	}

	defer script_residue.Close()

	t := template.Must(template.New("residue.sh").ParseFS(scriptTemplate, "templates/residue.sh"))
	if err := t.Execute(script_residue, collect); err != nil {
		panic(err)
	}

	log.Printf("[INFO] Identifying residue changes for %s...", read.Sample)
	run_residue := exec.Command("bash", path_script_residue)
	_, err = run_residue.Output()
	if err != nil {
		log.Fatalf("[ERR] Error identifying residue changes at %s, %s.", read.Sample, err.Error())
	}

	// remove the script
	log.Printf("[INFO] Completed collecting residue changes for %s.", read.Sample)
	os.Remove(path_script_residue)
}
