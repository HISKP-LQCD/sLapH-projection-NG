(* ::Package:: *)

(* This is the library which contains most of the functionality. The notebooks
are supposed to just to exercise this library. *)


BeginPackage["sLapHProjection`"];


ReadDataframe::usage = "Reads a CSV encoded data frame";

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


ReadEulerAngles[filename_] := Module[{oh, values},
	oh = ReadDataframe[filename];
	values = Normal[Values /@ oh];
	#[[1]]->{#[[2]], #[[3]], #[[4]]} & /@ values];


EndPackage[];
