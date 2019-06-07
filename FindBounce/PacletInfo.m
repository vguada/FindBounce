(* ::Package:: *)

(* Paclet Info File *)
(* BuildNumber and Internal values are inserted during build procedure. *)
Paclet[
	Name -> "FindBounce",
	Version -> "0.1.0",
	WolframVersion -> "10.+",
	Description -> "Computes tunneling transition with multiple scalar fields.",
	Creator -> "Victor Guada, Miha Nemevsek and Matevz Pintar",
	URL->"https://github.com/vguada/FindBounce",
	Thumbnail->"FrontEnd/Icon.png",
	Extensions -> {
		{"Kernel",
			Root -> ".",
			Context ->{"FindBounce`"}
		},
		{"Documentation",
			Language -> "English",
			MainPage -> "Tutorials/Examples_FindBounce"
		}
	}
]
