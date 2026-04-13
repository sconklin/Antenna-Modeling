# dumpnec

Tools for parsing and inspecting NEC2 antenna model input files.

## Files

### dumpnec

The original NEC2 file inspector. A standalone procedural script that parses a `.nec`
file and prints a human-readable description of each card immediately as it is
encountered. No in-memory model is built. Useful for quick inspection without any
flags or options.

**Usage:**

```
python3 dumpnec <filename.nec>
```

### necmodel.py

Parses a NEC2 `.nec` input file and builds an in-memory representation of the model.
Supports all standard geometry cards (GW, GA, GH, GM, GR, GS, GX, GE) and program
control cards (FR, EX, LD, GN, GD, NT, TL, NH, NE, RP, EK, KH).

**Usage:**

```
python3 necmodel.py <filename.nec> [--verbose] [--dump]
```

- `--verbose` / `-v`: Print each card as it is parsed
- `--dump` / `-d`: Print a structured summary of the entire model after parsing,
  in NEC card order (comments → geometry → program control → EN)

### validation-notes.txt

Notes on NEC2 modeling rules and constraints extracted from the NEC2 documentation.
Covers segment length limits, wire radius guidelines, junction rules, patch geometry
requirements, and other validity checks. Intended as a reference for a future
model-validation feature.
