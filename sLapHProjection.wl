(* ::Package:: *)

(* This is the library which contains most of the functionality. The notebooks
are supposed to just to exercise this library. *)

BeginPackage["sLapHProjection`"];


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

irrepFiles = FileNames["Single-cover/*-*-representations.txt"];
irrepDatasets = ReadIrreps /@ irrepFiles;
irrepAssocs = DatasetToAssocations /@ irrepDatasets;
momenta = ExtractMomentumFromFilename /@ irrepFiles;
IrrepDGammaAssoc[] = AssociationThread[momenta, irrepAssocs];


(* Irrep reading utilities that are not used *)

MatrixFromAssocList[assocs_] := Module[{dataset, rows, cols, values, size},
  dataset = Dataset[assocs];
  rows = Normal @ dataset[[All, "row"]];
  cols = Normal @ dataset[[All, "col"]];
  values = Normal @ dataset[[All, "matrix element"]];
  size = Length @ values;
  Partition[ToExpression /@ values, Sqrt[size]]];

MatrixLongToActual[dataset_] := GroupBy[
  dataset,
  Normal @ Values @ Take[#, 6] &,
  MatrixFromAssocList];


(* Cartesian representation *)

ReadEulerAngles[filename_] := Module[{oh, values},
  oh = ReadDataframe[filename];
  values = Normal[Values /@ oh];
  #[[1]] -> Pi * {ToExpression @ #[[2]], ToExpression @ #[[3]], ToExpression @ #[[4]]} & /@
    values // Association];
  
EulerAnglesAssoc[] = ReadEulerAngles["Single-cover/Oh-elements.txt"];

(* Momentum transformation *)

MomentumRefScalar[0] = {0,0,0};
MomentumRefScalar[1] = {0,0,1};
MomentumRefScalar[2] = {1,1,0};
MomentumRefScalar[3] = {1,1,1};
MomentumRefScalar[4] = {0,0,2};

MomentumRef[momentumpcm_] := MomentumRefScalar[Total[momentumpcm^2]];

EulerGTilde[momentumpcm_] := 
  First @ Select[Values @ EulerAnglesAssoc[],
    momentumpcm == EulerMatrix[#] . MomentumRef[momentumpcm] &];

MatrixRGTilde = EulerMatrix @* EulerGTilde;

MomentumTransform[momentumd_,  momentumpcm_, eulerG_] :=
  Inverse[MatrixRGTilde[MomentumRef[momentumpcm]]] .
    EulerMatrix[eulerG] .
    MatrixRGTilde[MomentumRef[momentumpcm]];


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
  Module[{eulerGtilde},
    eulerGtilde = EulerAngles[MatrixRGTilde[momentumpcm]];
    Sum[
      WignerD[{spinJi, spinMi1, spinMi}, eulerG[[1]], eulerG[[2]], eulerG[[3]]] 
      WignerD[{spinJi, spinMi2, spinMi1}, eulerGtilde[[1]], eulerGtilde[[2]], eulerGtilde[[3]]]
      ConjugateTranspose[SingleOperator[i, spinJi, spinMi2, MomentumTransform[momentumpi, momentumpcm, eulerG] . momentumpi]],
      {spinMi1, -spinJi, spinJi},
      {spinMi2, -spinJi, spinJi}]];

MakeMultiOperator[momentapi_, eulerG_, spinsJi_, spinsMi_] := Module[{momentumpcm, parts},
  momentumpcm = Total[momentapi];
  parts = MakeSingleOperator[momentapi[[#]], momentumpcm, eulerG, spinsJi[[#]], spinsMi[[#]], #] & /@
    Range[1, Length[momentapi]];
  NonCommutativeMultiply @@ parts];

MakeGroupSum[irrep_, irrepRow_, irrepCol_, momentapi_, spinsJi_, spinsMi_] := Module[{groupSummands},
  groupSummands = Module[{name, values, eulerG},
    name = Keys @ #;
    values = Values @ #;
    eulerG = EulerAnglesAssoc[][[Key @ name]];
    Conjugate[values[[Key @ {irrepRow, irrepCol}]]] *
    MakeMultiOperator[momentapi, eulerG, spinsJi, spinsMi]
  ] & /@ IrrepDGammaAssoc[][[Key @ Total @ momentapi]][[Key @ irrep]];
  Plus @@ groupSummands];

MakeMagneticSum2[irrep_, irrepRow_, irrepCol_, momentapi_, spinJ_, spinsJi_, phasePhiM_] :=
  Sum[
    phasePhiM[irrepCol] *
    Sum[
      ClebschGordan[{spinsJi[[1]], spinsMi1}, {spinsJi[[2]], spinsMi2}, {spinJ, spinM}] *
      MakeGroupSum[irrep, irrepRow, irrepCol, momentapi, 
      spinsJi, {spinsMi1, spinsMi2}],
    {spinsMi1, -spinsJi[[1]], spinsJi[[1]]},
    {spinsMi2, -spinsJi[[2]], spinsJi[[2]]}],
  {spinM, -spinJ, spinJ}];


(* Isospin and Spin *)

ReplaceSingleOperatorScalar[expr_, repl_] := 
  expr /. ConjugateTranspose[SingleOperator[1, 0, 0, p1_]] ** 
    ConjugateTranspose[SingleOperator[2, 0, 0, p2_]] -> ConjugateTranspose[repl[p2, p1]]

ReplaceSingleOperatorScalarWick[expr_, repl_] := 
  expr /.
    ConjugateTranspose[SingleOperator[1, 0, 0, p1_]] ** 
    ConjugateTranspose[SingleOperator[2, 0, 0, p2_]] ** 
    ConjugateTranspose[
      ConjugateTranspose[SingleOperator[1, 0, 0, p3_]] **
      ConjugateTranspose[SingleOperator[2, 0, 0, p4_]]] ->
    repl[p1, p2, p3, p4]


(* Contractions *)

MakeTemplate[n_] := StringRiffle[
 Table["g`g" <> ToString @ i <> "`.p`x" <> ToString @ i <> "`.d000", {i, 1, n}], "_"]

GammaRules[n_] := Association[Table[
  With[{v = "g" <> ToString @ i}, v :> ToExpresion @ v],
  {i, 1, n}]]

DatasetNameRules[] = {
  (* C4cD *)
  qct`trace[qct`Gamma^g1_ . qct`DE[{"up", "up"}, {si[si1_], so[so2_]}].
    qct`Gamma^g2_ . qct`DE[{"dn", "dn"}, {so[so2_], si[si1_]}]]
  qct`trace[qct`Gamma^g3_ .qct`DE[{"up", "up"}, {si[si3_], so[so4_]}].
    qct`Gamma^g4_ . qct`DE[{"dn", "dn"}, {so[so4_], si[si3_]}]] :> 
  TemplateApply[
    "C4cD_uuuu_" <> MakeTemplate[4],
    <|"g1" -> g1, "g2" -> g2, "g3" -> g3, "g4" -> g4,
      "x1" -> "`pso" <> ToString @ so2 <> "`",
      "x2" -> "`psi" <> ToString @ si1 <> "`",
      "x3" -> "`pso" <> ToString @ so4 <> "`",
      "x4" -> "`psi" <> ToString @ si3 <> "`"|>],

  (* C4cC *)
  qct`trace[qct`Gamma^g2_ . qct`DE[{"dn", "dn"}, {so[so2_], si[si1_]}].
    qct`Gamma^g3_ . qct`DE[{"up", "up"}, {si[si2_], so[so2_]}].
    qct`Gamma^g4_ . qct`DE[{"dn", "dn"}, {so[so1_], si[si2_]}].
    qct`Gamma^g1_ . qct`DE[{"up", "up"}, {si[si1_], so[so1_]}]] :>
  TemplateApply[
    "C4cC_uuuu_" <> MakeTemplate[4],
    <|"g1" -> g1, "g2" -> g2, "g3" -> g3, "g4" -> g4,
      "x1" -> "`pso" <> ToString @ so1 <> "`",
      "x2" -> "`psi" <> ToString @ si1 <> "`",
      "x3" -> "`pso" <> ToString @ so2 <> "`",
      "x4" -> "`psi" <> ToString @ si2 <> "`"|>]
}


EndPackage[];

(* vim: set cc=100 ft=mma sts=2 sw=2 :*)
