## Deletion

Many VCF files have the following deletion format, by using `.` for alternate allele.
```
POS  REF  ALT
100  G    .
```
This is not recommended, as one cannot use the length difference of
`REF` and `ALT` to tell whether it is a deletion or not.

VCF files from most variant callers (Mutect2, MuSE) have the following format:
```
POS  REF  ALT
100  GAA  G
```
where reference position `100` (1-based) is `G`

Convert to MAF deletion:
```
POS  REF  ALT
101  AA   -
```

## Insertion

The `REF` for VCF insertion should always have some reference base:
```
POS  REF  ALT
200  C    CTT
```
where reference position `200` (1-based) is `C`

Convert to MAF insertion:
```
POS  REF  ALT
200  -    TT
```
where it means `TT` is inserted ***after*** the `200`th base (`C`)
