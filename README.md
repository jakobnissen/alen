# Alen

Simple terminal sequence alignment viewer.

![Screenshot](/screenshots/prot.png?raw=true "Screenshot")

### What is Alen?
It's a command-like program to view DNA or protein alignments in FASTA formats. While more functionality will probably be added, Alen is intended to do the bare minimum of allowing you to view your alignments. It's meant for having a quick view of an alignment without having to leave the shell. It's not an alignment editor. It's not much, honestly.

### How to install and run.
Alen _should_ work on most Unix systems, and Windows 10. If someone asks me to, I might add Windows 7 and 8 support.
You can get it by downloading it from here and compiling it using Cargo. If other people begin using this tool, I might upload precompiled binaries.

### How to use
Simple usage:
```
$ alen /path/to/alignment.fasta
```

Note that Alen loads in the entire alignment in memory, so don't use it for multi-gigabyte files. Alen will auto-detect whether the alignment is protein or amino acid alignments. For more help, type `alen --help` in the terminal.

### Why the name?
First, it's to pun on the unimaginative named [alan](https://github.com/mpdunne/alan) and [alv](https://github.com/arvestad/alv). Second, _alen_ means _cubit_ in Danish. Like using cubits, Alen is simple, crude, but usually _good enough_.