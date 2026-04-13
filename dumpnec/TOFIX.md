# TOFIX — Known Issues in dumpnec / necmodel.py

Derived from comparison of NEC2 Part III User's Guide (`nec2prt3.pdf`) against
`dumpnec` and `necmodel.py`.

---

## Bugs that will crash at runtime

**1. `Arc.dump()` — AttributeError on `self.AngleEnd`**
`necmodel.py` line ~89
The attribute is defined as `self.angleEnd` (lowercase) but `dump()` references
`self.AngleEnd`. Will raise `AttributeError` whenever an arc's `dump()` is called.

**2. `Helix.addTaper()` — AttributeError on `self.wireRadius`**
`necmodel.py` line ~137
The Helix class stores wire radius as `self.radius`, not `self.wireRadius`. This
check will crash if any tapered helix is encountered.

---

## Bugs that produce wrong output

**3. `Helix` — middle two radii are swapped in display**
`necmodel.py` lines ~106-107; same error in `dumpnec`
The spec's GH field order is: S, HL, A1, B1, A2, B2, RAD — where A1/B1 are x/y
radii at z=0 and A2/B2 are x/y radii at z=HL. The code parses the fields at the
correct indices but assigns them to misnamed attributes:
- `parts[6]` (B1 = radius in y at z=0)  → stored as `radiusXZhl`
- `parts[7]` (A2 = radius in x at z=HL) → stored as `radiusYZ0`

`dump()` then labels them with swapped descriptions: "Radius in x @ z=HL" and
"Radius in y @ z=0" are printed for the wrong values.

**4. `EX` type 4 (elementary current source) — current moment never parsed, wrong value printed**
`necmodel.py` line ~408; `dumpnec` line ~327
The spec defines F4=alpha, F5=beta, F6=current moment for type 4. `Excitation.__init__`
reads `f4` and `f5` for type 4 but never reads `f6` (f6 is only read for types 1/2/3).
Both `dumpnec` and `Excitation.dump()` then print `f5` (the beta angle) labelled as
"Current Moment". The actual current moment is silently dropped.

**5. `EX` type 4 — `self.f6` never initialized**
Consequence of item 4: if code ever tried to access `self.f6` for a type-4 excitation
it would raise `AttributeError`. Currently masked because the dump code wrongly uses
`self.f5` instead.

---

## Cards present in spec but not meaningfully supported

All are in the `ops` dispatch table so the parser won't crash, but no useful output
is produced:

| Card | Current behavior |
|------|-----------------|
| GF (Read Numerical Green's Function) | `necmodel.py` raises `ValueError`; `dumpnec` prints "not supported" |
| CP (Maximum Coupling) | Both print "Not supported in nec2c" |
| NX (Next Structure) | Both print "Not supported in nec2c" |
| PQ (Print Control for Charge) | Both print "Not supported in nec2c" |
| PT (Print Control for Current) | Both print "Not supported in nec2c" |
| WG (Write NGF File) | Both print "Not supported in nec2c" |
| XQ (Execute) | Both print "Ignored in nec2c" |

Note: GF, WG, NX, CP, PT, PQ are valid NEC2 cards that xnec2c handles to varying
degrees — "not supported" may be overstated for some.

---

## Missing field / option coverage

**6. `GE` — second integer field I2 not parsed**
`necmodel.py` `GeometryEnd.__init__`
The spec defines a second integer field I2 on the GE card (used with the GF
Numerical Green's Function feature). It is not parsed.

**7. `EX` — I4 flag handling incomplete for plane wave types**
`necmodel.py` `Excitation.__init__` / `Excitation.dump()`
For types 1/2/3 the spec says I4 column-19 controls admittance asymmetry printing.
`dump()` only checks `i4` for voltage source types (0/5). For plane wave types, `i4`
is decoded via a single-digit `'%01d'` format where the spec calls for the same
two-digit treatment as voltage types.

**8. `GN` — cliff case condition is ambiguous**
`necmodel.py` `Ground.__init__` / `Ground.dump()`; same in `dumpnec`
Both implementations treat `NRADL == 0` as always meaning a two-medium cliff
problem. The spec says: if NRADL=0 *and* a second medium is desired, fields F3-F6
are used for cliff parameters; otherwise they are zero. The code unconditionally
prints cliff parameters when NRADL=0 even when all values are zero (no second
medium intended), as visible in the `dipole-17m.nec` output.

**9. `SP`/`SM`/`SC` — surface patch support incomplete**
Both files handle SP with the basic 8 fields, but the spec says a second SP card
follows the first for non-arbitrary patches (providing corners 3 and 4). Neither
implementation reads that second card. SM and SC are parsed in `necmodel.py` but
their data is never stored — there are no data model classes for surface patches,
and the SM/SC handlers only print during verbose mode without retaining anything.

**10. `GC` — I1 and I2 parsed but spec says they are unused**
`dumpnec` lines ~113-114; `necmodel.py` `Taper.__init__`
The spec says I1 and I2 on the GC card are not used. Both files read them anyway,
which is harmless but unnecessary.

---

## Minor labeling / description errors

- `GN` / `dumpnec`: "Dialectric" misspelled throughout (should be "Dielectric").
  Already corrected in `necmodel.py` `Ground.dump()` but not in `dumpnec`.
- `GC` in `dumpnec` line ~118: "Ratio to preior segment" — typo ("preior" → "prior").
- `PT` handler in `necmodel.py`: backslash instead of forward slash in
  "Page Title \ Print Control".
