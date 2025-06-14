---
output:
  pdf_document:
    latex_engine: xelatex
  html_document: default
---

# SUPPA2-DLERPA
DLERPA Bioinformatica project over alternatieve splicing met SUPPA2.

## Inhoud

- [Beschrijving](#beschrijving)
- [Benodigdheden](#benodigdheden)
- [Reproduceerbare omgeving](#reproduceerbare-omgeving)
- [Bestandenstructuur](#bestandenstructuur)
- [Data](#data)
- [Uitvoering](#uitvoering)
- [Resultaten](#resultaten)
- [Auteur en contact](#auteurs-en-contact)

---

## Beschrijving

Dit project richt zich op het verkennen van **SUPPA2**, een tool voor de analyse van alternatieve splicing op single-cell niveau. Het doel is om splicing events in kaart te brengen en te analyseren hoe deze verschillen tussen celtypen of experimentele omstandigheden.

Wat ik wil bereiken is het opdoen van kennis over SUPPA2 en het toepassen ervan op single-cell data, specifiek VASA-seq data.

**Hoofdvraag:**  
Hoe kan SUPPA2 worden gebruikt voor de analyse van alternatieve splicing in VASA-seq data, en in hoeverre is het een geschikte tool voor deze toepassing?

**Deelvragen:**  
- Hoe wordt VASA-seq data verwerkt met Seurat om cellen te filteren, clusteren en geschikt te maken voor analyse met SUPPA2?  
- Hoe werkt SUPPA2, welke input is vereist en hoe kan het worden toegepast op geclusterde VASA-seq data?  
- Hoe kunnen de PSI-waarden die SUPPA2 genereert worden geïnterpreteerd en gevisualiseerd om verschillen in alternatieve splicing tussen celclusters te analyseren?

---

## Benodigdheden

Voor dit project heb je zowel R als Python nodig. SUPPA2 werkt met Pythonm terwijl de pre-processing van de vasa-seq data met R/Seurat gebeurt.

### Software

- **R** (bij voorkeur versie 4.1 of hoger)
- **RStudio** (optioneel, maar aanbevolen)
- **Python** (verzie 3.6 of hoger)
- **SUPPA2** (installatie via Github: https://github.com/comprna/SUPPA)

### R-packages

De volgende R-packages worden gebruikt in dit project:

  - `Seurat` (versie 5)
  - `dplyr`
  - `ggplot2`
  - `patchwork`
  - `tidyr`
  - `colorspace`
  - `glmGamPoi` (aanbevolen voor snellere normalisatie met SCTranform)

Je kunt alle packages installeren met behulp van renv om de juiste versies te garanderen:

```r
renv::restore()
```

```r
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

if (!require("glmGamPoi")) {
  BiocManager::install("glmGamPoi")
}
```

---

## Bestandenstructuur

Het project is opgedeeld in twee hoofdmappen, elk met eigen data, scripts en output:

```
projectmap/
|-- README.md
|-- raw_data/
|   |-- e85_srr_list.txt
|   `-- fastq_files/
|-- seurat/
|   |-- data/
|   |-- output/
|   |-- scripts/
|   |-- seurat_tutorial.Rmd
|   `-- seurat.Rmd
|-- suppa2/
|   |-- data/
|   |-- output/
|   `-- scripts/
|   |   `-- run_download_fastq.sh
`-- suppa2.Rproj
```

---

## Data

De gebruikte VASA-seq data zijn opgeslagen op de HU-server en worden niet allemaal direct opgenomen in deze repository vanwege de bestandsgrootte.

Toegang tot de data kan verkregen worden via de server op het volgende pad:  

`/home/data/projecticum/splicing/data`  

Neem contact op met de projectleider of beheerder voor toegangsinformatie en rechten.

**Opmerkingen:**

Dit project is lokaal uitgevoerd, dus in de scripts wordt standaard verwezen naar paden binnen de projectmap.
Wanneer je het project op de server draait, wordt automatisch het serverpad gebruikt.

### Bestanden op de server

Serverpad:
`/home/data/projecticum/splicing/data` 

Bestanden:

- `e85_count_matrix.mtx.gz`
- `e85_feature_metadata.csv.gz`
- `e85_sample_metadata.csv`

Deze bestanden zijn gebruikt als input voor de Seurat-analyse, waarin kwaliteitscontrole, normalisatie en clustering zijn uitgevoerd.

Om dit reproduceerbaar te maken, is in de scripts de volgende code opgenomen om automatisch het juiste pad te kiezen:

```r
if (dir.exists("/home/data/projecticum/splicing/data")) {
  data_dir <- "/home/data/projecticum/splicing/data"
} else {
  data_dir <- file.path("seurat", "data")
}
```

### Extra data voor transcriptkwantificatie
Voor het analyseren van alternatieve splicing met SUPPA2 zijn andere typen data vereist, namelijk transcriptannotatie (.gtf) en kwantificaties in TPM. Deze TPM-waarden moeten eerst worden gegenereerd op basis van ruwe sequencingdata (FASTQ-bestanden)

FASTQ-bestanden zijn daarom gedownload met behulp van een eigen Bash-script: 

- Script: `suppa2/scripts/run_download_fastq.sh`

- De bijbehorende lijst met SRR-nummers: `raw_data/e85_srr_list.txt`

- De outputmap met FASTQ-bestanden: `raw_data/fastq_files/`

Het script controleert automatisch of bestanden al bestaan en slaat deze dan over. Zo wordt dubbel downloaden voorkomen.

---

## Uitvoering

Om de analyses uit te voeren, volg je de onderstaande stappen:

1. **Voorbereiding van data met Seurat**  
   - Ga naar de map `seurat/` en open het RMarkdown-bestand `seurat.Rmd`.  
   - Voer de workflow uit met Seurat versie 5 om de VASA-seq data te filteren en te clusteren.  
   - De output wordt opgeslagen in de map `seurat/output/`.

2. **Analyse met SUPPA2**  
   *Vervolgfase...*


Zorg ervoor dat alle benodigde software en packages geïnstalleerd zijn zoals beschreven in de sectie [Benodigdheden](#benodigdheden).

---

## Resultaten

De resultaten van dit project bestaan uit:

- Uitgevoerde Seurat analyses waarbij de VASA-seq data is gefilterd, geclusterd en voorbereid voor verdere analyse.  

*Vervolgfase:*

- SUPPA2 output met alternatieve splicing events en PSI-waarden per celcluster.   

- Visualisaties en interpretaties van splicing verschillen tussen celtypen, beschikbaar als figuren.

Voor een gedetailleerd overzicht van de analyses en resultaten, zie de RMarkdown-rapport (`seurat.Rmd`) en de gegenereerde outputbestanden.


--- 

## Auteur en contact

**Naam:** Jalisa van der Zeeuw  
**E-mail:** jalisavanderzeeuw@hotmail.com
**GitHub:** [https://github.com/ProjecticumDlerpDs/SUPPA2](https://github.com/ProjecticumDlerpDs/SUPPA2)

Voor vragen, suggesties of feedback kun je contact opnemen via e-mail of een issue aanmaken op de GitHub repository.

---