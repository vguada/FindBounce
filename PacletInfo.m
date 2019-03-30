(* ::Package:: *)

(* Paclet Info File *)
Paclet[
	Name -> "FindBounces",
	Version -> "0.0.1",
	WolframVersion -> "8.+",
	Description -> "Computes tunneling transition with multiple scalar fields.",
	Creator -> "Victor Guada",
	Publisher -> "  ",
	URL -> "https://github.com/",
	Thumbnail->"FrontEnd/Icon.png",
	Extensions -> {
		{"Kernel",
			Root -> ".",
			Context ->{"FindBounces`"}
		},
		{"Documentation",
			Language -> "English",
			MainPage -> "Guides/FindBounces"
		}(*,
		(* Metadata for PacletServer (https://paclets.github.io/PacletServer) *)
		{"PacletServer",
			"Tags" -> {"finite-elements","mesh","FEM"},
			"Categories" -> {"FEM"},
			"Description" -> "A package with utilities for creating and manipulating ElementMesh objects.",
			"License" -> "MIT"
		}*)
	}
]
