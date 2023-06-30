{
  description = "A basic flake for research project";
  nixConfig = {
    bash-prompt = "λ ";
  };
  inputs = {
    nixpkgs.url = "github:NixOS/nixpkgs/nixos-23.05";
    flake-utils.url = "github:numtide/flake-utils";
  };

  outputs = { self, nixpkgs, flake-utils }:
    flake-utils.lib.eachDefaultSystem (system: let
      pkgs = nixpkgs.legacyPackages.${system};
    in {
      devShells.default =  pkgs.mkShell {
        buildInputs = [

          # Testing
          pkgs.hyperfine

          # R and packages
          pkgs.R
          pkgs.rPackages.devtools
          pkgs.rPackages.languageserver

          # Documentation/writing tools
          pkgs.pandoc
        ];
      }; 
    });
}
