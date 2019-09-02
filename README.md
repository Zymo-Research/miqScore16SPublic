# 16S MIQ (Measurement Integrety Quotient) Calculator

*Terra incognita,* or land unknown, was a term used by cartographers to describe lands that were known to exist but had not yet been explored well enough to map.  While this phrase disappeard from maps over the last few centuries, it is ever-present in science. Biases in microbiome analysis have been known to exist since nearly the time the field began, but for many experiments, they remain *terra incognita.* Exploring this area provides us with a critical next step in the evolution of microbiome research and the constant drive to improve our work.

     There is always one more thing to learn.
                                 – Steve Jobs

At Zymo Research, our motto is, "*The beauty of science is to make things simple*," and it is with this principle guiding our efforts that we offer you this software package.
Like many great scientific tools, this package hides much of the complexity of its inner workings in order simplify a common analysis.  Under the hood of this system, you'll find [DADA2](https://github.com/benjjneb/dada2 "Github"), [FIGARO](https://github.com/Zymo-Research/Figaro, and the [NGS MIQ Score Calculator](https://github.com/Zymo-Research/miqScoreNGSReadCountPublic), the last two of which were initially developed for this application.  The mission of this application is identical to the [ZymoBIOMICS](https://www.zymoresearch.com/pages/zymobiomics-portfolio) mission: increasing the reproducibility and standardization of microbiome analysis.  Due to the multiple dependencies and need for reproducibility of results and environment in this software package, we have packaged it entirely within a [Docker](https://www.docker.com) image to ensure consistency between machines and runs.

#### Publication
Publication on this method is pending.  Watch this space for news and a link once it is available on bioRXiv.

## Quick Start Guide

This guide assumes you already have Docker installed and running on your system as well as the privileges to execute Docker commands.  To learn more about Docker, please check out [Katacoda](https://www.katacoda.com) for several excellent interactive tutorials on the use of Docker.

#### Downloading and building the image
```
git clone https://github.com/Zymo-Research/miqScore16SPublic.git
cd miqScore16SPublic
docker build -t miqscore16s .
```

#### Running the analysis
The simplest method to run this analysis (in this writer's opinion, at least) is to have the data ready and waiting in a file structure already set up to mount to the container.  Using this practice will ensure that the container finds the expected input files and will output data to the expected location as well as providing an easy system to script against for automation, with iteration only requiring the movement of files.

To do this, create the following directory structure somewhere on your system where you have read and write permissions:

```
.
+-- input
|  +-- sequence
|     +-- standard_submitted_R1.fastq
|     +- -standard_submitted_R2.fastq
+-- working
+-- output
```
This directory can have any name desired (**dataMountDirectory** will be used from here on as its name) so long as it conforms to good naming practices.  The files standard_submitted_R1.fastq and standard_submitted_R2.fastq are where this docker will look for read 1 and read 2 data, respectively, if no other files are specified.  These file names will work even if the fastq files are compressed using gzip.

Due to the nature of 16S sequencing and library preparation, you will need to have a little information on your library preparation method prior to running this container.  This information will consist of the maximum expected length of your amplicon and the lengths of both your forward and reverse primer.  You will also need a name that describes the sample and conforms to good file-naming practice (avoiding the use of spaces and special characters, underscores are acceptable)  A minimal command to run this container will be as follows, assuming that the container image was built under the name miqscore16s and the directory structure described above is in a directory called dataMountDirectory:

```
docker container run -v [pathTo]/dataMountDirectory:/data -e SAMPLENAME=My_Sample_Name -e FORWARDPRIMERLENGTH=16 -e REVERSEPRIMERLENGTH=24 -e AMPLICONLENGTH=510 miqscore16s
```

For more advanced usage, several parameters can be set using environment variables passed into the container at runtime. The environment variables that can be passed are as follows:

| Variable        | Type           | Default  | Description |
| --------------- |:--------------:|:--------:|-------------|
SAMPLENAME	|	string	|	**REQUIRED**	|	A name for the sample, should conform to good file-naming practices
FORWARDPRIMERLENGTH	|	integer	|	**REQUIRED**	|	Length of the forward primer
REVERSEPRIMERLENGTH	|	integer	|	**REQUIRED**	|	Length of the reverse primer
AMPLICONLENGTH	|	integer	|	**REQUIRED**	|	Maximum expected length for amplicon
MAXREADCOUNT	|	integer	|	0	|	Maximum reads to allow for analysis, will stop if a larger set of reads is provided. (Enter a value less than 1 for no limits)
FORWARDREADS	|	string	|	/data/input/sequence/standard_submitted_R1.fastq	|	Path to forward reads within the container (likely a mounted folder)
REVERSEREADS	|	string	|	/data/input/sequence/standard_submitted_R2.fastq	|	Path to forward reads within the container (likely a mounted folder)
SEQUENCEFOLDER	|	string	|	/data/input/sequence	|	Path to folder containing input sequences within the container
RSCRIPTFOLDER	|	string	|	[folderWithPackage]/rscripts	|	Path to folder containing the R scripts to run DADA2 (you will probably never use this one)
R1MAXEE	|	integer	|	3	|	Maximum expected error to allow for forward reads
R2MAXEE	|	integer	|	6	|	Maximum expected error to allow for reverse reads
TRUNCQ	|	integer	|	2	|	Truncate a read after this quality score
MAXMISMATCH	|	integer	|	2	|	Maximum mismatches allowed during read merge
MINOVERLAP	|	integer	|	10	|	Minimum desired overlap for read merging
NOCLEANUP	|	boolean	|	FALSE	|	Leave all intermediate files behind
OUTPUTFOLDER	|	string	|	/data/output	|	Folder within the container to write output data
DATABASEFILE	|	string	|	[folderWithPackage]/reference/rdp_train_set_16.fa.gz	|	File containing the rRNA reference sequences
TRIMPARAMETERDOWNSAMPLE	|	integer	|	-1	|	Downsampling for FIGARO trimming parameter prediction
TRIMPARAMETERPERCENTILE	|	integer	|	83	|	Expected error percentile for FIGARO to use when calculating trim parameters
FILENAMINGSTANDARD	|	string	|	ZYMO	|	How sequence files will be named (other option is "illumina")

## OUTPUT
All outputs will be written to the designated output folder for the container, which will unmount upon completion of the run.

There will be two primary output files, one report in HTML format and one in JSON format.  Examples of these can be seen [here for HTML](https://github.com/Zymo-Research/miqScore16SPublic/exampleReport.html) Examples of these can be seen [here for JSON](https://github.com/Zymo-Research/miqScore16SPublic/exampleReport.json).  The HTML report is designed to be viewed in a browser and gives an overview of the results for the sample.  The JSON report, while human-readable, is designed primarily to facilitate analysis using an automated script.  It also provides much more detailed information on the results than the HTML report.  
In addition to the two primary files, there will be a log file that can be used in the event of a problem with analysis for additional information on the run.  Finally, there will be several files generated by the DADA2 pipeline.  If you are familiar with DADA2, you will be familiar with these outputs.

## Contributing

We welcome and encourage contributions to this project from the microbiomics community and will happily accept and acknowledge input (and possibly provide some free kits as a thank you).  We aim to provide a positive and inclusive environment for contributors that is free of any harassment or excessively harsh criticism. Our Golden Rule: *Treat others as you would like to be treated*.

## Versioning

We use a modification of [Semantic Versioning](https://semvar.org) to identify our releases.

Release identifiers will be *major.minor.patch*

Major release: Newly required parameter or other change that is not entirely backwards compatible
Minor release: New optional parameter
Patch release: No changes to parameters

## Authors

- **Michael M. Weinstein** - *Project Lead, Programming and Design* - [michael-weinstein](https://github.com/michael-weinstein)
- **Aishani Prem** - *Testing, Design* - [AishaniPrem](https://github.com/AishaniPrem)
- **Mingda Jin** - *Testing, Code Review* - [jinmingda](https://github.com/jinmingda)
- **Shuiquan Tang** - *Design* - [shuiquantang](https://github.com/shuiquantang)
- **Jeffrey Bhasin** - *Design, Code Review* - [jeffbhasin](https://github.com/jeffbhasin)
- **Ryan Kemp**  - *Design* - [Zymo Research](https://www.zymoresearch.com)
- **Brian Janssen**  - *Design* - [Zymo Research](https://www.zymoresearch.com)

See also the list of [contributors](https://github.com/Zymo-Research/miqScore16SPublic/contributors) who participated in this project.

## License

This project is not currently licensed, but will likely be licensed under the GNU GPLv3 License - see the [LICENSE](LICENSE) file for details.
This license restricts the usage of this application for non-open sourced systems. Please contact the authors for questions related to relicensing of this software in non-open sourced systems.

## Acknowledgments

We would like to thank the following, without whom this would not have happened:
* The Python Foundation
* The staff at Zymo Research
* The microbiomics community
* Our customers

---------------------------------------------------------------------------------------------------------------------

#### If you like this software, please let us know at info@zymoresearch.com.
#### Please support our continued development of free and open-source microbiomics applications by checking out the latest microbiomics offerings from [ZymoBIOMICS](https://www.zymoresearch.com/pages/zymobiomics-portfolio)
