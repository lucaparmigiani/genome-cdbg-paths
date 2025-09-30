# genome-cdbg-paths

Print genomes as paths in a compacted de Bruijn graph (GFA format).

This tool maps input genome sequences onto a compacted de Bruijn graph and outputs their paths in standard GFA path format. 

It is compatible with graphs produced by tools such as **Bifrost**, **GGCAT**, **TwoPaCo**, **BCALM2**, **Cuttlefish2** and any other that generate GFA representations of compacted de Bruijn graphs.

`genome-cdbg-paths` can also be used as a preprocessing step for tools that take variation graphs as input, such as [panacus](https://github.com/codialab/panacus).

---

# Quick Start

## Install

### Dependencies

- **C++17 compiler**   
- **OpenMP**   
- **zlib** 

Clone the repository and build:

```bash
git clone https://github.com/yourname/genome-cdbg-paths.git
cd genome-cdbg-paths
make
```

## Run

Basic usage:

```bash
genome-cdbg-paths <k> <graph.gfa> <genome.fasta> [more.fasta ...] [options]
```

Options:

* `--break`:  split paths at non-ACGT characters
* `-t, --threads <N>`: number of threads (default: all available cores)

---

## Example

Run with three genomes:

```bash
genome-cdbg-paths 11 example/graph.gfa example/a.fa example/b.fa example/c.fa > paths.txt
```

Example output (in GFA path format):

```
P	a#0#0	10+,11+,12+,12-,13+,13-,12+,12-,13+,13-,12+,12-,13+,13-,12+,12-,13+,13-,12+,12-,13+,13-,12+,12-,11-,1+,2+,14-,15-,3+,2+,14-,15-,3+,2+,14-,15-	*
P	b#0#0	15-,3+,2+,14-,15-,4+,5+,6+,15-,3+,2+,14-,15-,4+,16+,7+,15-,3+,2+,14-,15-,4+,16+,8+,14-,15-,3+,2+,14-,15-,3+,2+,14-,15-,4+	*
P	c#0#0	9+,17-,6+,15-,3+,2+,14-,15-,4+,5+,6+,15-,3+,2+,14-,15-,4+,16+,7+,15-,3+,2+,14-,15-,4+,16+,8+,14-,15-,3+,2+,14-,15-,3+,2+,14-,15-,4+	*
```

Each `P` line describes a genome as a walk through the compacted de Bruijn graph.

You can then add the paths to the GFA (watchout to not overwrite your `graph.gfa`):

```bash
cat paths.txt >> example/graph.gfa
```

### Path naming convention

Path names follow this convention:

```
file_name#id_sequence#split
```

* `file_name`: the source FASTA file name
* `id_sequence`: 0-based index of the sequence within the FASTA file
* `split`: index of the segment if `--break` was used (splits a sequence into multiple paths at non-ACGT characters)
