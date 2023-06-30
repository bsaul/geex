{
  description = "geex: m-estimation in R";
  nixConfig = {
    bash-prompt = "ψ ";
  };
  inputs = {
    nixpkgs.url = "github:NixOS/nixpkgs/nixos-23.05";
    flake-utils.url = "github:numtide/flake-utils";
    gitignore = {
      url = "github:hercules-ci/gitignore.nix";
      # Use the same nixpkgs
      inputs.nixpkgs.follows = "nixpkgs";
    };
  };

  outputs = { self, nixpkgs, flake-utils, gitignore }:
    flake-utils.lib.eachDefaultSystem (system: let
      pkgs = nixpkgs.legacyPackages.${system};

      geexRDeps = with pkgs.rPackages; [Matrix rootSolve numDeriv lme4];

      inherit (gitignore.lib) gitignoreSource;

    in {

      packages.default = pkgs.rPackages.buildRPackage {
        name = "geex";
        src  = gitignoreSource ./.;
        propagatedBuildInputs = geexRDeps;
      };

      devShells.default =  pkgs.mkShell {
        buildInputs = [

          # Testing
          pkgs.hyperfine

          # R and packages
          pkgs.R
          pkgs.rPackages.devtools
          pkgs.rPackages.languageserver
          pkgs.rPackages.lintr
          pkgs.rPackages.styler
          
          # Documentation/writing tools
          pkgs.pandoc
        ] ++ geexRDeps;
      }; 
    });
}
