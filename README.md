# Alen

Simple terminal sequence alignment viewer.

![Screenshot](/screenshots/prot.png?raw=true "Screenshot")

### What is Alen?
It's a command-line program to view DNA or protein alignments in FASTA formats. Alen is meant for having a quick view of an alignment without having to leave the shell. It's not an alignment editor.

### How to install and run.
Alen _should_ work on most Unix systems, and Windows 10. If someone asks me to, I might add Windows 7 and 8 support.

To get started with Alen, you can install a pre-compiled binary via Conda or download one from the GitHub releases page.

To install Alen using Conda:
```
$ conda install -c bioconda alen
```

Visit the GitHub [releases page](https://github.com/jakobnissen/alen/releases) to download the binary (binaries are not yet available for Windows).

Alternatively, you can download the source code from this repo and compile it yourself using Cargo:
```
$ cargo build --release
```

### How to use
Simple usage:
```
$ alen /path/to/alignment.fasta
```

Note that Alen loads in the entire alignment into memory, so don't use it for multi-gigabyte files. Alen will auto-detect whether the alignment is a nucleotide or amino acid. For more help, type `alen --help` in the terminal.

__Commands__

* `Esc / q / Ctrl-C`: Quit
* `c`: Toggle consensus sequence comparison
* `r`: Re-render the screen.
* `Ctrl-f`: Find. Searches headers, then sequences for regex, case insensitively.
* `Ctrl-j`: Jump to column.
* `Ctrl-s`: Select rows and move them around.
* `Arrow keys`: Move 1 column/row
* `Shift-Arrow keys`: Move 10 columns/rows
* `Ctrl-Arrow keys`: Move to first/last column/row

### Configuration
No.

### Why the name?
First, it's to pun on the unimaginative named [alan](https://github.com/mpdunne/alan) and [alv](https://github.com/arvestad/alv). Second, _alen_ means _cubit_ in Danish. Like using cubits, Alen is simple, crude, but usually _good enough_.
