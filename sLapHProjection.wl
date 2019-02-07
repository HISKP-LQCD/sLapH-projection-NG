(* ::Package:: *)

(* This is the library which contains most of the functionality. The notebooks
are supposed to just to exercise this library. *)


BeginPackage["sLapHProjection`"];


(* Text IO helpers *)


ReadDataframe[path_] := Module[{assocs, bulk, data, header},
	data = Import[path, "CSV"];
	header = data[[1]];
	bulk = Drop[data, 1];
	(* https://mathematica.stackexchange.com/q/88418/1507 *)
	assocs = AssociationThread[header, #] & /@ bulk;
	Dataset[assocs]];


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
IrrepDGammaAssoc = AssociationThread[momenta, irrepAssocs];


ReadEulerAngles[filename_] := Module[{oh, values},
	oh = ReadDataframe[filename];
	values = Normal[Values /@ oh];
	#[[1]] -> Pi * {#[[2]], #[[3]], #[[4]]} & /@ values // Association];
	
EulerAnglesAssoc = ReadEulerAngles["Single-cover/Oh-elements.txt"];


(* Momentum transformation *)


MomentumRefScalar[0] = {0,0,0};
MomentumRefScalar[1] = {0,0,1};
MomentumRefScalar[2] = {1,1,0};
MomentumRefScalar[3] = {1,1,1};
MomentumRefScalar[4] = {0,0,2};

MomentumRef[momentumpcm_] := MomentumRefScalar[Total[momentumpcm^2]];


EulerGTilde[momentumpcm_] := 
	First @ Select[Values @ EulerAnglesAssoc,
		momentumpcm == EulerMatrix[#] . MomentumRef[momentumpcm] &];

MatrixRGTilde = EulerMatrix @* EulerGTilde;

MomentumTransform[momentumd_,  momentumpcm_, eulerG_] :=
	Inverse[MatrixRGTilde[MomentumRef[momentumpcm]]] .
		EulerMatrix[eulerG] .
		MatrixRGTilde[MomentumRef[momentumpcm]];


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


MakeGroupSum[irrep_, momentapi_, spinsJi_, spinsMi_] := Module[{groupSummands},
	groupSummands = Module[{name, values, eulerG},
		name = Keys @ #;
		values = Values @ #;
		eulerG = EulerAnglesAssoc[[Key @ name]];
		MakeMultiOperator[momentapi, eulerG, spinsJi, spinsMi]
	] & /@ IrrepDGammaAssoc[[Key @ Total @ momentapi]][[Key @ irrep]];
	Plus @@ groupSummands];


EndPackage[];
