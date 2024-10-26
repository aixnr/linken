package process

import (
	"embed"
	"linken-hermes/model"
	"log"
	"os"
	"os/exec"
	"text/template"
)

//go:embed templates/*
var scriptTemplate embed.FS

func ShellCaller(script_file string, collect model.Collect) {
	script, err := os.Create(script_file)
	if err != nil {
		log.Fatalf("Failed to create %s at %s", script_file, collect.Sample)
	}

	defer script.Close()

	tpl_embed := "templates/" + script_file
	t := template.Must(template.New(script_file).ParseFS(scriptTemplate, tpl_embed))

	if err := t.Execute(script, collect); err != nil {
		log.Fatal(err)
	}

	script_run := exec.Command("bash", script_file)
	_, err = script_run.Output()
	if err != nil {
		log.Printf("Error processing for %s at %s", collect.Sample, script_file)
		log.Fatal(err)
	}

	os.Remove(script_file)
}
