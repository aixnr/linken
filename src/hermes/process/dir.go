package process

import (
	"log"
	"os"
)

func CreateDir() {
	// fastqc does not automate folder creation
	err := os.Mkdir("output_fastqc", 0700)
	if err != nil {
		log.Fatalf("[ERR] Failed to create 'output_fastqc' directory, %s", err.Error())
	}

	// To avoid littering project root, move mapped artifacts to output_map
	// Moving is handled by templates/call.sh at ProcessCalibrateAndCall()
	err = os.Mkdir("output_map", 0700)
	if err != nil {
		log.Fatalf("[ERR] Failed to create 'output_map' directory, %s.", err.Error())
	}

	// ProcessCoverage() dumps artifacts into 'output_coverage' directory (handled by coverage.sh)
	err = os.Mkdir("output_coverage", 0700)
	if err != nil {
		log.Fatalf("[ERR] Failed to create 'output_coverage' directory, %s", err.Error())
	}

	// ProcessExtractVariantToTable dumps artifacts into 'output_variant' (handled by variant-tsv.sh)
	err = os.Mkdir("output_variant", 0700)
	if err != nil {
		log.Fatalf("[ERR] Failed to create 'output_variant' directory, %s", err.Error())
	}

	// eevee requires output_plots to dump all figures
	err = os.Mkdir("output_plots", 0700)
	if err != nil {
		log.Fatalf("[ERR] Failed to create 'output_plots' directory, %s", err.Error())
	}
}
