(* ::Package:: *)

(* This is the library which contains most of the functionality. The notebooks
are supposed to just to exercise this library. *)

BeginPackage["sLapHProjection`"];

Needs["qct`"]

(* Utility functions *)

basePath = DirectoryName @ DirectoryName @ $InputFileName;

(*
MonitoredMap[f_, list_, label_ : ""] := Module[{i},
  Monitor[Table[f[list[[i]]], {i, 1, Length[list]}],
   Row[{ProgressIndicator[i, {1, Length[list] + 1}],
     TemplateApply[" `` `` of ``", {label, i, Length[list]}]}, " "]]];
*)
MonitoredMap[f_, list_, label_ : ""] := Map[f, list];


(* Reading the lattice irreps *)

ReadDataframe[path_] := Module[{assocs, bulk, data, header},
  data = Import[path, "CSV"];
  header = data[[1]];
  bulk = Drop[data, 1];
  (* https://mathematica.stackexchange.com/q/88418/1507 *)
  assocs = AssociationThread[header, #] & /@ bulk;
  Dataset[assocs]];

ConvertMatrixElement[dataset_] := Module[{normal},
  normal = Normal[dataset];
  normal[[All, "matrix element"]] = ToExpression /@ Normal @ normal[[All, "matrix element"]];
  Dataset[normal]];

ReadIrreps = ConvertMatrixElement @* ReadDataframe;

MatrixElementsToAssociation[dataset_] := GroupBy[
  dataset,
  Normal @ Values @ Take[#, 6] &,
  Association[{#[["row"]], #[["col"]]} -> #[["matrix element"]] & /@ #] &]

IrrepsToAssociations[matrixElementsAssociations_] := Module[{keysIrrep, keysIrrep2, irrep3},
  keysIrrep = Normal @ GroupBy[Keys @ matrixElementsAssociations, #[[1]] &];
  keysIrrep2 = #[[2]] -> # & /@ # & /@ keysIrrep;
  irrep3 = matrixElementsAssociations[Key @ #] & /@ Association @ # & /@ keysIrrep2;
  Map[Normal, irrep3, Infinity]]

DatasetToAssocations := IrrepsToAssociations @* MatrixElementsToAssociation;

ExtractMomentumFromFilename[filename_] := First[ToExpression /@ StringCases[
  filename,
  RegularExpression["\\((-?\\d+),(-?\\d+),(-?\\d+)\\)"] -> {"$1", "$2", "$3"}]]

irrepFiles = FileNames[basePath <> "single_cover/*-*-representations.txt"];
irrepDatasets = ReadIrreps /@ irrepFiles;
irrepAssocs = DatasetToAssocations /@ irrepDatasets;
momenta = ExtractMomentumFromFilename /@ irrepFiles;
IrrepDGammaAssoc[] = AssociationThread[momenta, irrepAssocs];


(* Cartesian representation *)

ReadEulerAngles[filename_] := Module[{oh, values},
  oh = ReadDataframe[filename];
  values = Normal[Values /@ oh];
  #[[1]] -> Pi * {ToExpression @ #[[2]], ToExpression @ #[[3]], ToExpression @ #[[4]]} & /@
    values // Association];

Parities[] =
  First @* Values /@
    Association @ IrrepDGammaAssoc[][[Key @ {0, 0, 0}]][["A1u"]];

EulerAnglesAssoc[] = ReadEulerAngles[basePath <> "single_cover/Oh-elements.txt"];

EulerAnglesParityAssoc[] = MapThread[List, {EulerAnglesAssoc[], Parities[]}];


(* Momentum transformation *)

MomentumRefScalar[0] = {0, 0, 0};
MomentumRefScalar[1] = {0, 0, 1};
MomentumRefScalar[2] = {1, 1, 0};
MomentumRefScalar[3] = {1, 1, 1};
MomentumRefScalar[4] = {0, 0, 2};

MomentumRef[momentumpcm_] := MomentumRefScalar[Total[momentumpcm^2]];

ParityEulerMatrix[{angles_, parity_}] := parity * EulerMatrix[angles];

CachedParityEulerMatrix = AssociationMap[ParityEulerMatrix, Values @ EulerAnglesParityAssoc[]];

EulerGTilde[momentumpcm_] :=
  First @ Select[Values @ EulerAnglesParityAssoc[],
    momentumpcm == CachedParityEulerMatrix[#] . MomentumRef[momentumpcm] &];

MatrixRGTilde[momentumpcm_] :=
  CachedParityEulerMatrix @ EulerGTilde @ momentumpcm;

GetParity[momentumpcm_, angles_] :=
  a /. Solve[a EulerMatrix[angles] . momentumpcm == momentumpcm][[1]];

MomentumTransform[matrixRGtilde_, eulerG_] :=
  matrixRGtilde . CachedParityEulerMatrix[eulerG] . Inverse @ matrixRGtilde;

RotateMomenta[groupElementName_String, momenta_] := With[
  {rotationMatrix =
    CachedParityEulerMatrix @
     EulerAnglesParityAssoc[][[groupElementName]]},
  Map[rotationMatrix . # &, momenta]]

MomentaOrbit[momenta_] :=
  Sort @ DeleteDuplicates @ Map[Sort @ RotateMomenta[#, momenta] &,
     Keys @ IrrepDGammaAssoc[][[Key @ MomentumRef @ Total @ momenta]][[1]]];

RemoveRedundantMomenta[individualMomenta_] :=
  Keys @ DeleteDuplicates @ AssociationMap[MomentaOrbit, individualMomenta];

UniqueTotalMomenta[momentumMag_] :=
  DeleteDuplicates @
    Values[# . MomentumRefScalar[momentumMag] & /@
      (EulerMatrix[#[[1]]] * #[[2]] &) /@ EulerAnglesParityAssoc[]];

RelativeToIndividualMomenta[totalMomentum_, relMomenta_] :=
  MatrixRGTilde[totalMomentum] . # & /@
    Catenate[{{MomentumRef[totalMomentum] - Total[relMomenta]}, relMomenta}];

RelMomentaFromIndividual[momenta_] := Drop[momenta, 1];

RelMomentaRef[totalMomentum_, relMomenta_] :=
  Inverse[MatrixRGTilde[totalMomentum]] . # & /@ relMomenta;

RelMomentaRefFromIndividual[momenta_] :=
  RelMomentaRef[Total @ momenta, RelMomentaFromIndividual @ momenta];

MomentaToString[momenta_] := StringRiffle[ Map[MomentumToString, momenta], ","];

RelMomentaRefLabelFromIndividual[momenta_] :=
  MomentaToString @ RelMomentaRefFromIndividual @ momenta;

MomentaMaxNormSq[momenta_] := Max[Norm[#]^2 & /@ momenta];

MomentaSumNormSq[momenta_] := Total[Norm[#]^2 & /@ momenta];

ContractionMomentumCutoff[0] = 4;
ContractionMomentumCutoff[1] = 5;
ContractionMomentumCutoff[2] = 6;
ContractionMomentumCutoff[3] = 7;
ContractionMomentumCutoff[4] = 4;

FilterRelativeMomenta[totalMomentum_, relMomenta_] :=
  RelMomentaRefFromIndividual /@
    Select[
      RemoveRedundantMomenta @ Map[RelativeToIndividualMomenta[totalMomentum, #] &, relMomenta],
      MomentaSumNormSq @ # <= ContractionMomentumCutoff[Norm[totalMomentum]^2] &&
      MomentaMaxNormSq @ # <= 4 &];


(* Clebsch-Gordan coefficients *)

HigherClebschGordan[{j1_, j2_}, {m1_, m2_}, {j_, m_}] :=
  If[Abs[j1 - j2] <= j <= j1 + j2 && m1 + m2 == m,
    ClebschGordan[{j1, m1}, {j2, m2}, {j, m}], 0];

HigherClebschGordan[js_, ms_, {j_, m_}] := Sum[
  HigherClebschGordan[{js[[1]], jtilde}, {ms[[1]], mtilde}, {j, m}] *
  HigherClebschGordan[Drop[js, 1], Drop[ms, 1], {jtilde, mtilde}],
    {jtilde, Abs[js[[2]] - js[[3]]],
  js[[2]] + js[[3]]},
  {mtilde, -jtilde, jtilde}];


(* Spin *)

MakeSingleOperator[momentumpi_, momentumpcm_, eulerG_, spinJi_, spinMi_, i_] :=
  Module[{eulerGtilde, matrixRGtilde, momentumTransform},
    eulerGtilde = EulerGTilde[momentumpcm];
    matrixRGtilde = CachedParityEulerMatrix @ eulerGtilde;
    momentumTransform = MomentumTransform[matrixRGtilde, eulerG];
    Sum[
      eulerG[[2]] *
      WignerD[{spinJi, spinMi1, spinMi}, eulerG[[1]][[1]], eulerG[[1]][[2]], eulerG[[1]][[3]]] *
      eulerGtilde[[2]] *
      WignerD[{spinJi, spinMi2, spinMi1}, eulerGtilde[[1]][[1]], eulerGtilde[[1]][[2]], eulerGtilde[[1]][[3]]] *
      ConjugateTranspose[SingleOperator[i, spinJi, spinMi2, momentumTransform . momentumpi]],
      {spinMi1, -spinJi, spinJi},
      {spinMi2, -spinJi, spinJi}]];

MakeMultiOperator[momentapi_, eulerG_, spinsJi_, spinsMi_] := Module[{momentumpcm, parts},
  momentumpcm = Total[momentapi];
  parts = MakeSingleOperator[momentapi[[#]], momentumpcm, eulerG, spinsJi[[#]], spinsMi[[#]], #] & /@
    Range[1, Length[momentapi]];
  NonCommutativeMultiply @@ parts];

MakeGroupSum[irrep_, irrepRow_, irrepCol_, momentapi_, spinsJi_, spinsMi_] := Module[{groupSummands},
  groupSummands = Map[
    Module[{name, values, eulerG},
      name = Keys @ #;
      values = Values @ #;
      eulerG = EulerAnglesParityAssoc[][[Key @ name]];
      Conjugate[values[[Key @ {irrepRow, irrepCol}]]] *
      MakeMultiOperator[momentapi, eulerG, spinsJi, spinsMi]] &,
    IrrepDGammaAssoc[][[Key @ MomentumRef @ Total @ momentapi]][[Key @ irrep]]];
  Plus @@ groupSummands / Length[groupSummands]];

MakeMagneticSum[irrep_, irrepRow_, irrepCol_, momentapi_, spinJ_, spinsJi_, phasePhiM_] :=
  Sum[
    phasePhiM[irrepCol] *
    Sum[
      ClebschGordan[{spinsJi[[1]], spinsMi1}, {spinsJi[[2]], spinsMi2}, {spinJ, spinM}] *
      MakeGroupSum[irrep, irrepRow, irrepCol, momentapi,
      spinsJi, {spinsMi1, spinsMi2}],
    {spinsMi1, -spinsJi[[1]], spinsJi[[1]]},
    {spinsMi2, -spinsJi[[2]], spinsJi[[2]]}],
  {spinM, -spinJ, spinJ}];

ExtractMomenta[expr_] :=
  expr /. ConjugateTranspose[SingleOperator[_, _, _, p_]] ->
    DTMomentum[p]

MomentumProductToMomenta[factors_] := Module[
  {momenta, scalar},
  momenta = factors /. a_. DTMomentum[p_] -> p;
  scalar = Times @@ (factors /. a_. DTMomentum[p_] -> a);
  scalar * DTMomenta @@ momenta];

mm[factors_] :=
  Apply[Times, #[[All, 1]]] DTMomenta[#[[All, 2]]] &[
    factors /. a_. DTMomentum[v_?VectorQ] :> {a, v}]l

ExtractMultiMomenta[assoc_] := Block[{x1, x2},
  x1 = Map[ExtractMomenta, assoc, {5}];
  x2 = x1 /. NonCommutativeMultiply[a__] :> MomentumProductToMomenta[{a}];
  Map[Evaluate, x2, {5}]];

MomentumToString[p_] := StringJoin[ToString /@ p];

momentaToRules[momenta_, location_] :=
  ReplaceAll[momenta,
    DTMomenta[p__] :>
    AssociationThread[
      Table["p" <> location <> ToString @ i, {i, 1, Length[{p}]}],
      MomentumToString /@ {p}]];

MomentaToAssoc[expr_, location_, sign_] :=
  expr /. DTMomenta[p__] :> DTMomentaAssoc[momentaToRules[DTMomenta @@ (sign * # & /@ {p}), location]];

MomentaToAssocSourceSink[expr1_, expr2_] := Module[
  FullSimplify @ ReplaceAll[
  ExpandAll @ ReplaceAll[
    ExpandAll[Conjugate @ MomentaToAssoc[expr1, "so", +1] * MomentaToAssoc[expr2, "si", -1]],
    Conjugate[Plus[a__]] :> Plus @@ Conjugate /@ a],
  Conjugate[DTMomentaAssoc[<|a__|>]] * DTMomentaAssoc[<|b__|>] :> DTMomentaAssoc[<|a, b|>]]];

MomentaToAssocSourceSink[expr1_, expr2_] := Module[{
  assocSo = MomentaToAssoc[expr1, "so", +1],
  assocSi = MomentaToAssoc[expr2, "si", -1],
  expandedProduct,
  expandedProductReplaced},
  expandedProduct = ExpandAll[Conjugate@assocSo*assocSi];
  (* https://mathematica.stackexchange.com/a/109735/1507 *)

expandedProductReplaced = ExpandAll @
    ReplaceAll[expandedProduct, e : Conjugate[Plus[__]] :> Thread[e, Plus]];
  FullSimplify @ ReplaceAll[expandedProductReplaced,
    Conjugate[DTMomentaAssoc[<|a__|>]]*DTMomentaAssoc[<|b__|>] :> DTMomentaAssoc[<|a, b|>]]];

MakeSourceSinkMomenta[assoc_] := Module[
  {keys = Keys @ assoc,
    values = Values @ assoc,
    labels,
    outer},
  labels = StringRiffle[#, ","] & /@ keys;
  outer = Outer[MomentaToAssocSourceSink, values, values];
  AssociationThread[labels, AssociationThread[labels, #] & /@ outer]];


(* Isospin *)

\[Pi]Plus[s1_, s2_, c1_, x1_] :=
  -qct`FieldB["up", c1, s1, x1] **
  (qct`Gamma^5)[qct`SI[{s1, s2}]] **
  qct`Field["dn", c1, s2, x1];

\[Pi]Minus[s1_, s2_, c1_, x1_] :=
  qct`FieldB["dn", c1, s1, x1] **
  (qct`Gamma^5)[qct`SI[{s1, s2}]] **
  qct`Field["up", c1, s2, x1];

\[Pi]Zero[s1_, s2_, c1_, x1_] :=
  (qct`FieldB[up, c1, s1, x1]
  ** (qct`Gamma^5)[qct`SI[{s1, s2}]]
  ** qct`Field[ up, c1, s2, x1]
  - qct`FieldB[dn, c1, s1, x1]
  ** (qct`Gamma^5)[qct`SI[{s1, s2}]]
  ** qct`Field[dn, c1, s2, x1]) / Sqrt[2];

\[Pi]\[Pi]I1[s1_, s2_, s3_, s4_, c1_, c2_, x1_, x2_] :=
  (\[Pi]Plus[s1, s2, c1, x1] ** \[Pi]Minus[s3, s4, c2, x2] -
  \[Pi]Plus[s3, s4, c2, x2] ** \[Pi]Minus[s1, s2, c1, x1]) / Sqrt[2];

\[Pi]\[Pi]I1Bar[s1_, s2_, s3_, s4_, c1_, c2_, x1_, x2_] :=
  (\[Pi]Minus[s3, s4, c2, x2] ** \[Pi]Plus[s1, s2, c1, x1] -
  \[Pi]Minus[s1, s2, c1, x1] ** \[Pi]Plus[s3, s4, c2, x2]) / Sqrt[2];

\[Pi]\[Pi]I2[s1_, s2_, s3_, s4_, c1_, c2_, x1_, x2_] :=
  \[Pi]Plus[s1, s2, c1, x1] ** \[Pi]Plus[s3, s4, c2, x2];

\[Pi]\[Pi]I2Bar[s1_, s2_, s3_, s4_, c1_, c2_, x1_, x2_] :=
  \[Pi]Minus[s3, s4, c2, x2] ** \[Pi]Minus[s1, s2, c1, x1];

\[Pi]\[Pi]\[Pi]I3[s1_, s2_, s3_, s4_, s5_, s6_, c1_, c2_, c3_, x1_, x2_, x3_] :=
  \[Pi]Plus[s1, s2, c1, x1] ** \[Pi]Plus[s3, s4, c2, x2] ** \[Pi]Plus[s5, s6, c3, x3]; 

\[Pi]\[Pi]\[Pi]I3Bar[s1_, s2_, s3_, s4_, s5_, s6_, c1_, c2_, c3_, x1_, x2_, x3_] :=
  \[Pi]Minus[s5, s6, c3, x3] ** \[Pi]Minus[s3, s4, c2, x2] ** \[Pi]Minus[s1, s2, c1, x1]; 


(* Trace normalization *)

RotateGammaToFront[expr_] := If[MatchQ[expr[[1]], qct`Gamma^_], expr, RotateLeft[expr]];

Starts[traceContent_] :=
  ReplaceAll[traceContent[[2 ;; ;; 2]], {Dot -> List}];

StartScore[propagator_] := propagator /.
  {"up" -> 1000, "dn" -> 2000, so[i_] -> i, si[i_] -> i + 100} /.
  prop[f_, t_] -> f + t;

IndexOfFirst[traceContent_] :=
  Ordering[Starts[traceContent], 1,
    StartScore[#1] < StartScore[#2] &][[1]];

NormalizeTrace[traceContent_] := With[{tr2 = RotateGammaToFront @ traceContent},
  RotateLeft[tr2, 2 * (IndexOfFirst[tr2] - 1)]];

NormalizeTraceRecursive[expr_] := If[AtomQ @ expr,
  expr,
  If[MatchQ[expr, qct`trace[_]],
    qct`trace[NormalizeTrace[expr[[1]]]],
    NormalizeTraceRecursive /@ expr]];

ReplacePropagators[expr_] := expr /. qct`DE[{f_, f_}, {_, t_}] :> prop[f, t];

FlowReversalRules[] = {
  qct`trace[qct`Gamma^5 . prop["up", si[si2_]].qct`Gamma^5 . prop["dn", so[so1_]]] :>
    qct`trace[qct`Gamma^5 . prop["up", so[so1]] . qct`Gamma^5 . prop["dn", si[si2]]]
  };

FlavorSwitchRules[] = {
  qct`trace[qct`Gamma^5 . prop["up", so[so1_]] . qct`Gamma^5 . prop["dn", so[so2_]] . qct`Gamma^5 . prop["up", si[si2_]] . qct`Gamma^5 . prop["dn", si[si1_]]] :>
    Conjugate @ qct`trace[qct`Gamma^5 . prop["up", so[so1]] . qct`Gamma^5 . prop["dn", si[si1]] . qct`Gamma^5 . prop["up", si[si2]] . qct`Gamma^5 . prop["dn", so[so2]]]
  };

WickContractionToTemplates[wc_] := With[
  {normalized = NormalizeTraceRecursive @ ReplacePropagators @ QuarkContract @ wc},
  normalized /. FlavorSwitchRules[] /. FlowReversalRules[] /. DatasetNameRules[]];


(* Wick contractions *)

MakeTemplate[n_] := StringRiffle[
 Table["p`x" <> ToString @ i <> "`.d000.g`g" <> ToString @ i <> "`", {i, 1, n}], "_"]

DatasetNameRules[] = {
  (* C6cD *)
  qct`trace[qct`Gamma^g1_ . prop["up", so[so1_]] .
    qct`Gamma^g2_ . prop["dn", si[si1_]]] *
  qct`trace[qct`Gamma^g3_ . prop["up", so[so2_]] .
    qct`Gamma^g4_ . prop["dn", si[si2_]]] *
  qct`trace[qct`Gamma^g5_ . prop["up", so[so3_]] .
    qct`Gamma^g6_ . prop["dn", si[si3_]]] :>
  TemplateApply[
    "C6cD_uuuuuu_" <> MakeTemplate[6],
    <|"g1" -> g1, "g2" -> g2, "g3" -> g3, "g4" -> g4, "g5" -> g5, "g6" -> g6,
      "x1" -> "`pso" <> ToString @ so1 <> "`",
      "x2" -> "`psi" <> ToString @ si1 <> "`",
      "x3" -> "`pso" <> ToString @ so2 <> "`",
      "x4" -> "`psi" <> ToString @ si2 <> "`",
      "x5" -> "`pso" <> ToString @ so3 <> "`",
      "x6" -> "`psi" <> ToString @ si3 <> "`"|>],

  (* C6cC *)
  qct`trace[qct`Gamma^g1_ . prop["up", so[so1_]] .
    qct`Gamma^g2_ . prop["dn", si[si1_]] .
    qct`Gamma^g3_ . prop["up", so[so2_]] .
    qct`Gamma^g4_ . prop["dn", si[si2_]] .
    qct`Gamma^g5_ . prop["up", so[so3_]] .
    qct`Gamma^g6_ . prop["dn", si[si3_]]] :>
  TemplateApply[
    "C6cC_uuuuuu_" <> MakeTemplate[6],
    <|"g1" -> g1, "g2" -> g2, "g3" -> g3, "g4" -> g4, "g5" -> g5, "g6" -> g6,
      "x1" -> "`pso" <> ToString @ so1 <> "`",
      "x2" -> "`psi" <> ToString @ si1 <> "`",
      "x3" -> "`pso" <> ToString @ so2 <> "`",
      "x4" -> "`psi" <> ToString @ si2 <> "`",
      "x5" -> "`pso" <> ToString @ so3 <> "`",
      "x6" -> "`psi" <> ToString @ si3 <> "`"|>],

  (* C6cCD *)
  qct`trace[qct`Gamma^g5_ . prop["up", so[so3_]].
    qct`Gamma^g6_ . prop["dn", si[si3_]]] *
  qct`trace[qct`Gamma^g1_ . prop["up", so[so1_]].
    qct`Gamma^g2_ . prop["dn", si[si1_]].
    qct`Gamma^g3_ . prop["up", so[so2_]].
    qct`Gamma^g4_ . prop["dn", si[si2_]] ] :>
  TemplateApply[
    "C4cCD_uuuuuu_" <> MakeTemplate[4],
    <|"g1" -> g1, "g2" -> g2, "g3" -> g3, "g4" -> g4, "g5" -> g5, "g6" -> g6,
      "x1" -> "`pso" <> ToString @ so1 <> "`",
      "x2" -> "`psi" <> ToString @ si1 <> "`",
      "x3" -> "`pso" <> ToString @ so2 <> "`",
      "x4" -> "`psi" <> ToString @ si2 <> "`",
      "x4" -> "`pso" <> ToString @ so3 <> "`",
      "x5" -> "`psi" <> ToString @ si3 <> "`"|>],

  (* C4cB *)
  qct`trace[qct`Gamma^g1_ . prop["up", so[so1_]].
    qct`Gamma^g2_ . prop["dn", si[si2_]].
    qct`Gamma^g3_ . prop["up", si[si3_]].
    qct`Gamma^g4_ . prop["dn", so[so4_]]] :>
  TemplateApply[
    "C4cB_uuuu_" <> MakeTemplate[4],
    <|"g1" -> g1, "g2" -> g2, "g3" -> g3, "g4" -> g4,
      "x1" -> "`pso" <> ToString @ so1 <> "`",
      "x2" -> "`psi" <> ToString @ si2 <> "`",
      "x3" -> "`psi" <> ToString @ si3 <> "`",
      "x4" -> "`pso" <> ToString @ so4 <> "`"|>],

  (* C4cC *)
  qct`trace[qct`Gamma^g1_ . prop["up", so[so1_]].
    qct`Gamma^g2_ . prop["dn", si[si1_]].
    qct`Gamma^g3_ . prop["up", so[so2_]].
    qct`Gamma^g4_ . prop["dn", si[si2_]] ] :>
  TemplateApply[
    "C4cC_uuuu_" <> MakeTemplate[4],
    <|"g1" -> g1, "g2" -> g2, "g3" -> g3, "g4" -> g4,
      "x1" -> "`pso" <> ToString @ so1 <> "`",
      "x2" -> "`psi" <> ToString @ si1 <> "`",
      "x3" -> "`pso" <> ToString @ so2 <> "`",
      "x4" -> "`psi" <> ToString @ si2 <> "`"|>],

  (* C4cD *)
  qct`trace[qct`Gamma^g1_ . prop["up", so[so1_]].
    qct`Gamma^g2_ . prop["dn", si[si2_]]]
  qct`trace[qct`Gamma^g3_ . prop["up", so[so3_]].
    qct`Gamma^g4_ . prop["dn", si[si4_]]] :>
  TemplateApply[
    "C4cD_uuuu_" <> MakeTemplate[4],
    <|"g1" -> g1, "g2" -> g2, "g3" -> g3, "g4" -> g4,
      "x1" -> "`pso" <> ToString @ so1 <> "`",
      "x2" -> "`psi" <> ToString @ si2 <> "`",
      "x3" -> "`pso" <> ToString @ so3 <> "`",
      "x4" -> "`psi" <> ToString @ si4 <> "`"|>],

  (* C4cV *)
  qct`trace[qct`Gamma^g3_ . prop["up", si[si3_]].
    qct`Gamma^g4_ . prop["dn", si[si4_]]]
  qct`trace[qct`Gamma^g1_ . prop["up", so[so1_]].
    qct`Gamma^g2_ .prop["dn", so[so2_]]] :>
  TemplateApply[
    "C4cV_uuuu_" <> MakeTemplate[4],
    <|"g1" -> g1, "g2" -> g2, "g3" -> g3, "g4" -> g4,
      "x1" -> "`pso" <> ToString @ so1 <> "`",
      "x2" -> "`pso" <> ToString @ so2 <> "`",
      "x3" -> "`psi" <> ToString @ si3 <> "`",
      "x4" -> "`psi" <> ToString @ si4 <> "`"|>]
}


(* Isospin and Spin *)

CombineIsospinAndSpin[corrTemplates_, momentaAssoc_] :=
  momentaAssoc /. DTMomentaAssoc[rules_] :>
    ReplaceAll[corrTemplates, str_String :> TemplateApply[str, rules]];


(* Export to data frame *)

(* https://mathematica.stackexchange.com/a/191718/1507 *)
StringExpressionToAssociation[expr_] := Module[
  {exprConj, keys},
  exprConj = expr /. Conjugate[str_String] :> "conj:" <> str;
  keys = Union @ Cases[exprConj, _String, Infinity];
  AssociationMap[Coefficient[exprConj, #, 1] &, keys]];

NeedsConjugation[name_] := With[{prefix = StringTake[name, 5]},
  If[prefix == "conj:",
    {StringDrop[name, 5], True},
    {name, False}]];

(* https://mathematica.stackexchange.com/a/191718/1507 *)
DatasetnameAssocToCSV[assoc_, filename_String] := With[
  {export = KeyValueMap[Flatten @ {NeedsConjugation @ #1, N @ ReIm @ #2} &, assoc] //
    ExportString[#, "CSV"] &},
  DeleteFile[filename];
  WriteString[filename, "datasetname,conjugate,re,im\n" <> export]];


(* GEVP building *)

IrrepSize[totalMomentum_, irrep_] :=
  Last @ Last @ Keys @ First @ Association @
    IrrepDGammaAssoc[][totalMomentum][[irrep]];

MultiGroupSum[irrep_, momentapi_, irrepRow_, irrepCol_, hold_ : Identity] :=
  hold @ MakeGroupSum[irrep, irrepRow, irrepCol, momentapi, {0, 0, 0}, {0, 0, 0}];

GroupSumIrrepRowCol[totalMomentum_, irrep_, irrepRow_, irrepCol_, relMomenta_, hold_ : Identity] :=
  With[{selectedIndividualMomenta = RelativeToIndividualMomenta[totalMomentum, #] & /@ relMomenta},
    AssociationThread[Map[MomentumToString, relMomenta, {2}],
      MonitoredMap[MultiGroupSum[irrep, #, irrepRow, irrepCol, hold] &,
        selectedIndividualMomenta, "Momentum"]]];

GroupSumIrrepRow[totalMomentum_, irrep_, irrepCol_, relMomenta_, hold_ : Identity] :=
Module[{rows = Range[1, IrrepSize[totalMomentum, irrep]]},
  AssociationThread[
    ToString /@ rows,
    MonitoredMap[GroupSumIrrepRowCol[totalMomentum, irrep, #, irrepCol, relMomenta, hold] &,
      rows, "Irrep row"]]];

GroupSumWholeIrrep[totalMomentum_, irrep_, relMomenta_, hold_ : Identity] :=
Module[{cols = Range[1, IrrepSize[totalMomentum, irrep]]},
  AssociationThread[
    ToString /@ cols,
    MonitoredMap[GroupSumIrrepRow[totalMomentum, irrep, #, relMomenta, hold] &,
      cols, "Irrep col"]]];

GroupSumWholeTotalMomentum[totalMomentum_, relMomenta_, hold_ : Identity] := Module[
  {irreps = Keys @ IrrepDGammaAssoc[][totalMomentum]},
  AssociationThread[irreps,
    MonitoredMap[GroupSumWholeIrrep[totalMomentum, #, relMomenta, hold] &,
      irreps, "Irrep"]]];

DatasetnameToObject[value_, key_] := Module[
  {datasetname = key[[1]] /. Key[x_] -> x, nc},
  nc = NeedsConjugation[datasetname];
  <|"datasetname" -> nc[[1]], "re" -> Re @ value, "im" -> Im @ value, "conj" -> nc[[2]]|>];

DatasetnameAssocToObject[value_] :=
  Values @ MapIndexed[DatasetnameToObject, value];

DropEmpty[container_] := Select[container, Length @ # > 0 &];

PrescriptionToNumeric[prescription_] := <|
  "datasetname" -> prescription["datasetname"],
  "re" -> N @ prescription["re"],
  "im" -> N @ prescription["im"],
  "conj" -> prescription["conj"]|>;

PrescriptionFilename[totalMomentum_, irrep_] :=
  "prescriptions/prescription_" <> MomentumToString[totalMomentum] <> "_" <> irrep <> ".js";

MomentaAndTemplatesToJSONFile[momentaAssoc_, templates_, filename_] := Module[
  {someMomenta, someSourceSinkMomenta, gevp1, gevp2, gevp222, gevp3, gevp3C1, gevp3C2, gevpNumeric, json},
  someMomenta = ExtractMultiMomenta[momentaAssoc];
  someSourceSinkMomenta = Map[MakeSourceSinkMomenta, someMomenta, {4}];
  gevp1 = Map[CombineIsospinAndSpin[templates, #] &, someSourceSinkMomenta, {6}];
  gevp2 = Map[StringExpressionToAssociation, gevp1, {6}];
  gevp222 = Map[DatasetnameAssocToObject, gevp2, {6}];
  gevp3 = AssociationThread[MomentumToString /@ Keys[gevp222], Values[gevp222]];
  gevp3C1 = Map[DropEmpty, gevp3, {5}];
  gevp3C2 = Map[DropEmpty, gevp3C1, {4}];
  gevpNumeric = Map[PrescriptionToNumeric, gevp3C2, {7}];
  json = ExportString[gevpNumeric, "JSON"];
  If[FileExistsQ[filename], DeleteFile[filename], Null];
  WriteString[filename, json];
  json];

StructureButSingle[totalMomentum_, irrep_, relMomenta_] :=
  <|totalMomentum -> <|irrep ->
    GroupSumWholeIrrep[totalMomentum, irrep, relMomenta]|>|>;


EndPackage[];

(* vim: set cc=100 ft=mma sts=2 sw=2 :*)
