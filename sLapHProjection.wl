(* ::Package:: *)

(* This is the library which contains most of the functionality. The notebooks
are supposed to just to exercise this library. *)


BeginPackage["sLapHProjection`"];


ReadDataframe[path_] := Module[{assocs, bulk, data, header},
	data := Import[path, "CSV"];
	header := data[[1]];
	bulk := Drop[data, 1];
	(* https://mathematica.stackexchange.com/q/88418/1507 *)
	assocs := AssociationThread[header, #] & /@ bulk;
	Dataset[assocs]
]


MatrixFromAssocList[assocs_] := Module[{dataset, rows, cols, values, size},
	dataset = Dataset[assocs];
	rows = Normal @ dataset[[All, "row"]];
	cols = Normal @ dataset[[All, "col"]];
	values = Normal @ dataset[[All, "matrix element"]];
	size = Length @ values;
	Partition[ToExpression /@ values, Sqrt[size]]]


MatrixLongToActual[dataset_] := GroupBy[
	dataset,
	Normal @ Values @ Take[#, 6] &,
	MatrixFromAssocList]


EndPackage[];
