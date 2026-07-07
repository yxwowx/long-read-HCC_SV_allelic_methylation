# 04_admr_definition

Recurrently remodeled constitutional aDMR loci and somatic aDMR definition.

Two distinct, non-mutually-exclusive populations built from different parent
sets — do not compare counts across them directly:
- `constitutional/`: tumor aDMR requiring bulk tumor-normal DMR overlap (>=30%
  reciprocal, >=100bp); further stratified by cross-patient recurrence depth.
- `somatic/`: tumor aDMR with <10% overlap to any matched-normal aDMR (constitutional
  signal explicitly excluded).