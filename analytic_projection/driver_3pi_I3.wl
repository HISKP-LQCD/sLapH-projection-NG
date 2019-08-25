srcDir = DirectoryName @ $InputFileName;
PrependTo[$Path, srcDir];

Needs["qct`"]
Needs["sLapHProjection`"]

Print @ $ScriptCommandLine

totalMomentum = ToExpression /@ $ScriptCommandLine[[2 ;; 4]]
irrep = $ScriptCommandLine[[5]]

wc = WickContract[
  \[Pi]\[Pi]\[Pi]I3Bar[s1, s2, s3, s4, s5, s6, c1, c2, c3, so[1], so[2], so[3]] **
  \[Pi]\[Pi]\[Pi]I3[s7, s8, s9, s10, s11, s12, c4, c5, c6, si[1], si[2], si[3]]];

templates = WickContractionToTemplates @ wc;

utm1 = UniqueTotalMomenta /@ Range[0, 1];
utm1Flat = Flatten[utm1, 1]
utm2 = UniqueTotalMomenta /@ Range[0, 1];
utm2Flat = Flatten[utm2, 1]
relMomenta = Flatten[Outer[List, utm2Flat, utm2Flat, 1], 1]

cutoffRelative = FilterRelativeMomenta[totalMomentum, relMomenta];

timing1 = AbsoluteTiming[some = StructureButSingle[totalMomentum, irrep, cutoffRelative]][[1]];
Print @ timing1

filename = PrescriptionFilename[totalMomentum, irrep];
timing2 = AbsoluteTiming[MomentaAndTemplatesToJSONFile[some, templates, filename]][[1]];
Print @ timing2
