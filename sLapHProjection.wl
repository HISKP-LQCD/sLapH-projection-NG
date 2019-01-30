(* ::Package:: *)

BeginPackage["sLapHProjection`"];


ReadDataframe[path_] := Module[{assocs, bulk, data, header},
	data := Import[path, "CSV"];
	header := data[[1]];
	bulk := Drop[data, 1];
	(* https://mathematica.stackexchange.com/q/88418/1507 *)
	assocs := AssociationThread[header, #] & /@ bulk;
	Dataset[assocs]
]


EndPackage[];
