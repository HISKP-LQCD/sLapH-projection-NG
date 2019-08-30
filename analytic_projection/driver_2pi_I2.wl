srcDir = DirectoryName @ $InputFileName;
PrependTo[$Path, srcDir];

Needs["qct`"]
Needs["sLapHProjection`"]

Print @ $ScriptCommandLine

totalMomentum = ToExpression /@ $ScriptCommandLine[[2 ;; 4]]
irrep = $ScriptCommandLine[[5]]

wc = WickContract[
  \[Pi]\[Pi]I2Bar[s1, s2, s3, s4, c1, c2, so[1], so[2]] **
  \[Pi]\[Pi]I2[s7, s8, s9, s10, c4, c5, si[1], si[2]]];

templates = WickContractionToTemplates @ wc;

utm2 = UniqueTotalMomenta /@ Range[0, 1];
utm2Flat = Flatten[utm2, 1];
relMomenta = {#} & /@ utm2Flat

cutoffRelative = FilterRelativeMomenta[totalMomentum, relMomenta];

timing1 = AbsoluteTiming[some = StructureButSingle[totalMomentum, irrep, cutoffRelative]][[1]];
Print @ timing1

filename = PrescriptionFilename[totalMomentum, irrep];
timing2 = AbsoluteTiming[MomentaAndTemplatesToJSONFile[some, templates, filename]][[1]];
Print @ timing2
